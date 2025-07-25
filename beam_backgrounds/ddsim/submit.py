import sys, os, glob, shutil
import time
import argparse
import logging
import subprocess
import random
import socket

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action='store_true', help="Submit to batch system")
parser.add_argument("--gun_file", type=str, help="Input gun file")
parser.add_argument("--gp_dir", type=str, help="Input directory containing .pairs files from Guinea pig")
parser.add_argument("--geometry_file", type=str, help="Geometry file (locally)")
parser.add_argument("--geometry_tag", type=str, help="Geometry tag (from Key4hep)")
parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="longlunch")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/user/j/jaeyserm/fccee/beam_backgrounds")
parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--nJobs", help="Maximum number of events", type=int, default=9e99)
parser.add_argument("--maxMemory", help="Maximum job memory", type=float, default=2000)
parser.add_argument("-x", "--xrootd", action='store_true', help="Use XrootD transfer")
args = parser.parse_args()

# /cvmfs/sw.hsf.org/key4hep/releases/2024-03-10/x86_64-almalinux9-gcc11.3.1-opt/k4geo/0.20-hapqru/share/k4geo/FCCee/CLD/compact/CLD_o2_v05/
# $K4GEO/FCCee/CLD/compact/CLD_o2_v05/
# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml --nJobs 10 --cms_pool
# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/guineaPig_andrea_June2024_v23/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --osg_pool --nJobs 10 
# python submit.py --gun_file cards/gun/mu_theta_0-180_p_50.input --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml --osg_pool --nJobs 10 
# python submit.py --gun_file cards/gun/mu_theta_0-180_p_50.input --geometry $K4GEO/FCCee/IDEA/compact/CLD_o2_v07/CLD_o2_v07.xml --osg_pool --nJobs 10 
# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig_tracking//FCCee_Z_4IP_04may23_FCCee_Z32_grids1/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --nJobs 10 --osg_pool
# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig_tracking//FCCee_Z_4IP_04may23_FCCee_Z32_grids1/ --geometry_tag $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --nJobs 10 --osg_pool

GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10" # CLD_o2_v05
GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29" # CLD_o2_v07
#GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-01-28" # IDEA_o1_v03
#SINGULARITY = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/key4hep/k4-deploy/alma9:latest"
SINGULARITY = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"
HOSTNAME = socket.gethostname()


def get_voms_proxy_path():
    try:
        output = subprocess.check_output(['voms-proxy-info'], text=True)
        for line in output.splitlines():
            if line.strip().startswith('path'):
                return line.split(':', 1)[1].strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running voms-proxy-info: {e}")
    return None

pdg_map = {
    13: 'mu-',
    -13: 'mu+',
    11: 'e-',
    -11: 'e+',
    22: 'gamma',
    211: 'pi+',
    -211: 'pi-',
}

class DDSimProducer:

    def __init__(self, storagedir, args):
        self.geometry_tag = None
        self.geometry_file = None
        if args.geometry_tag:
            self.geometry_tag = args.geometry_tag # automatically expanded with env vars
            self.geometry_name = os.path.splitext(os.path.basename(self.geometry_tag))[0]
            if not os.path.exists(self.geometry_tag):
                logger.error(f"Geometry file {self.geometry_tag} does not exist, exiting")
                sys.exit(1)
        elif args.geometry_file:
            self.geometry_file = os.path.abspath(os.path.expandvars(args.geometry_file))
            self.geometry_dir = os.path.dirname(self.geometry_file)
            self.geometry_name = os.path.splitext(os.path.basename(self.geometry_file))[0]

            if not os.path.exists(self.geometry_file):
                logger.error(f"Geometry file {self.geometry_file} does not exist, exiting")
                sys.exit(1)
        self.storagedir = storagedir
        self.cwd = os.getcwd()

        logger.info(f"Found geometry {self.geometry_name}")

        self.local_dir = f"{self.cwd}/local_ddsim/{self.geometry_name}/"

        self.condor_queue = args.condor_queue
        self.condor_priority = args.condor_priority
        self.args = args
        self.ddsim_args = ""

    def generate_gp_submit(self):

        gun, gp = False, False
        if args.gun_file and os.path.exists(args.gun_file):
            gun = True
            self.gun_file = os.path.abspath(args.gun_file)
            self.name = os.path.basename(os.path.normpath(self.gun_file)).split('.')[0]
            logger.info(f"Found gun file {self.name}")

            multiplicity, particle, nevents = None, None, None
            thetaMin, thetaMax, momentumMin, momentumMax = None, None, None, None
            with open(self.gun_file , 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    key, value = line.split(None, 1) # split at first whitespace
                    if key == 'npart':
                        multiplicity = value
                    if key == 'nevents':
                        nevents = value
                    if key == 'theta_range':
                        thetaMin, thetaMax = value.split(',')
                    if key == 'mom_range':
                        momentumMin, momentumMax = value.split(',')
                    if key == 'pid_list':
                        if ',' in value:
                            logger.error(f"ddsim does not support a list of particles {value}")
                            quit()
                        particle = pdg_map[int(value)]

            self.ddsim_args = f"--enableGun --gun.particle {particle} --gun.multiplicity {multiplicity} --gun.thetaMin '{thetaMin}*deg' --gun.thetaMax '{thetaMax}*deg' --gun.momentumMin '{momentumMin}*GeV' --gun.momentumMax '{momentumMax}*GeV' --random.seed ${{seed}} --numberOfEvents {nevents} --gun.distribution uniform"

            seeds = [random.randint(10000, 99999) for _ in range(int(args.nJobs))]

        elif args.gp_dir and os.path.exists(args.gp_dir):
            gp = True
            self.gp_dir = os.path.abspath(args.gp_dir)
            self.name = os.path.basename(os.path.normpath(self.gp_dir))
            logger.info(f"Found guineapig directory {self.name}")

            logger.info(f"Get .pairs files from {self.gp_dir}")
            input_files = glob.glob(f"{self.gp_dir}/*.pairs")
            if len(input_files) == 0:
                logger.error(f"Input directory {self.gp_dir} does not contain .pair files, exiting")
                sys.exit(1)
            logger.info(f"Found {len(input_files)} pairs files")


            # get seeds
            seeds = []
            for n, f in enumerate(input_files):
                if n >= args.nJobs:
                    break
                seed = os.path.basename(f).replace("output_", "").replace("output0_", "").replace(".pairs", "")
                seeds.append(seed)
            self.ddsim_args = f"--inputFiles output_${{seed}}.pairs --numberOfEvents -1"

        else:
            logger.error(f"Please provide input guinea pig directory or a gun file")
            quit()



        # preparing submission dirs
        self.out_dir = f"{self.storagedir}/{self.geometry_name}/{self.name}/"
        self.log_dir = f"{self.storagedir}/{self.geometry_name}/{self.name}/logs/"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)


        # pack the geometry files to be shipped
        if self.geometry_file != None:
            logger.info("Packing geometry files")
            geometry_sandbox = f"{self.out_dir}/geometry.tar"
            os.system(f"cd {self.geometry_dir} && tar -cf {geometry_sandbox} .")
            compactFile = f"geometry/{self.geometry_name}.xml"
            geometry_action = "mkdir geometry && tar -xf geometry.tar -C geometry"
        else:
            compactFile = self.geometry_tag
            geometry_action = ""


        out_dir_xrd = self.out_dir.replace("/ceph/submit", "")

        if args.xrootd:
           xrootd_cfg = f"""
#echo "Setup certificates directories"
#ls -lrt /etc/grid-security/

if [ -d /cvmfs/grid.cern.ch/etc/grid-security/certificates ]; then
    export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
else
    echo "ERROR: Grid certs not available on CVMFS" >&2
    exit 1
fi


#echo $X509_USER_PROXY

#source /cvmfs/grid.cern.ch/alma9-ui-current/etc/profile.d/setup-alma9-test.sh
#echo "Copy output file"
#echo $X509_USER_PROXY
#export X509_USER_PROXY=x509
#echo $X509_USER_PROXY

voms-proxy-info -all
#voms-proxy-info -all -file x509

xrdcp --version
xrdcp -d 3 output_${{seed}}_sim.root root://submit50.mit.edu/{out_dir_xrd}/output_${{seed}}_sim.root
rc=$?
if [ $rc -ne 0 ]; then
    echo "xrdcp failed with exit code $rc" >&2
    #exit $rc
fi
echo "Copy output file via XrootD done"
           """
        else:
            xrootd_cfg = ""

        job = f"""
#!/bin/bash

set -e


unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

echo "Release:"
cat /proc/version

echo "Hostname:"
hostname

echo "List current working dir"
ls -lrt


echo "Set seed"
export seed=$1
echo $seed



# unpack geometry
{geometry_action}


echo "Checking CVMFS..."
if ! timeout 15 ls /cvmfs/sw.hsf.org/key4hep/setup.sh >/dev/null 2>&1; then
    echo "ERROR: CVMFS or key4hep setup.sh not found!" >&2
    exit 1
fi


echo "CVMFS looks OK. Sourcing Key4Hep"
if ! source {GP_STACK} >/dev/null 2>&1; then
    echo "ERROR: Failed to source Key4Hep setup" >&2
    exit 1
fi
echo "Key4Hep sourcing successful"


echo "Start ddsim"
SECONDS=0
#ddsim --compactFile {compactFile} --outputFile output_${{seed}}_sim.root --inputFiles output_${{seed}}.pairs --numberOfEvents -1
ddsim --compactFile {compactFile} --outputFile output_${{seed}}_sim.root {self.ddsim_args} --crossingAngleBoost 0.015
#ddsim --compactFile {compactFile} --outputFile output_${{seed}}_sim.root {self.ddsim_args}
duration=$SECONDS

#echo "lol" > output_${{seed}}_sim.root

# check ddsim exit code
# fail job in case ddsim isn't properly executed --> try again via max_retries
# on e.g. OSG, sometimes segfaults because of shared object libs not found
rc=$?
if [ $rc -ne 0 ]; then
    echo "ddsim failed with exit code $rc" >&2
    exit $rc
fi

# check if output file exists
if [ ! -f output_${{seed}}_sim.root ]; then
    echo "ddsim did not produce output_${{seed}}_sim.root" >&2
    exit 1
fi

echo "Done ddsim"
ls -lrt


{xrootd_cfg}


echo "Duration: $(($duration))"

        """



        # make executable script
        submitFn = f"{self.out_dir}/run_ddsim.sh"
        fOut = open(submitFn, "w")
        fOut.write(job)

        subprocess.getstatusoutput(f"chmod 777 {submitFn}")

        # make condor submission script
        condorFn = f'{self.out_dir}/condor.cfg'
        fOut = open(condorFn, 'w')

        fOut.write(f'universe       = vanilla\n')
        fOut.write(f'initialdir     = {self.out_dir}\n')
        fOut.write(f'executable     = {submitFn}\n')
        fOut.write(f'arguments      = $(SEED)\n')

        fOut.write(f'Log            = logs/condor_job.$(ClusterId).$(ProcId).log\n')
        fOut.write(f'Output         = logs/condor_job.$(ClusterId).$(ProcId).out\n')
        fOut.write(f'Error          = logs/condor_job.$(ClusterId).$(ProcId).error\n')

        fOut.write(f'should_transfer_files = YES\n')
        fOut.write(f'when_to_transfer_output = ON_EXIT\n')
        if gp:
            if self.geometry_file:
                fOut.write(f'transfer_input_files = {self.gp_dir}/output_$(SEED).pairs,{geometry_sandbox}\n')
            else:
                fOut.write(f'transfer_input_files = {self.gp_dir}/output_$(SEED).pairs\n')
        else:
            if self.geometry_file:
                fOut.write(f'transfer_input_files = {geometry_sandbox}\n')

        if args.xrootd:
            fOut.write(f'transfer_output_files = ""\n') # explicit no transfer files back
        else:
            fOut.write(f'transfer_output_files = output_$(SEED)_sim.root\n') # done by xrdcp

        fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        fOut.write(f'max_retries    = 150\n')
        
        # retry the job if it failed due to output transfer failure (HoldReasonCode == 12)
        # e.g. when job fails, it doesn't produce the output file, going to Hold (need to prevent it)
        fOut.write(f'periodic_release = (HoldReasonCode == 12 && NumJobStarts < 150)\n')
        fOut.write(f'on_exit_hold = False\n')

        fOut.write(f'RequestMemory  = {args.maxMemory}\n')

        if args.xrootd:
            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.log_dir}/")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = logs/{os.path.basename(proxy_path)}\n')

        # global pool
        if args.cms_pool:
            fOut.write(f'+DESIRED_Sites = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"\n')

            fOut.write(f'+AccountingGroup      = "analysis.jaeyserm"\n')
            fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')

        # OSG pool
        elif args.osg_pool:
            # https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/requirements/#additional-feature-specific-attributes
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            ##fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_SINGULARITY == TRUE )\n')




        elif 'mit.edu' in HOSTNAME:
            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.log_dir}/x509")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = logs/x509\n')
            
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')
            fOut.write(f'+DESIRED_Sites = "mit_tier2,mit_tier3"\n')
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            #fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')

        elif 'cern.ch' in HOSTNAME:
            fOut.write(f'+JobFlavour    = "{self.condor_queue}"\n')
            fOut.write(f'+AccountingGroup = "{self.condor_priority}"\n')

        seedsStr = ' \n '.join([str(s) for s in seeds])
        fOut.write(f'queue SEED in ( \n {seedsStr} \n)\n')

        fOut.close()

        subprocess.getstatusoutput(f'chmod 777 {condorFn}')
        os.system(f"condor_submit {condorFn}")

        logger.info(f"Written to {self.out_dir}")




def main():
    producer = DDSimProducer(args.storagedir, args)
    producer.generate_gp_submit()


if __name__ == "__main__":
    main()
