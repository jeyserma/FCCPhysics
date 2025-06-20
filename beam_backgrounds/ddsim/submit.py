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
parser.add_argument("--geometry_file", type=str, help="Geometry file", required=True)
parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="longlunch")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/ee/detector/ddsim/")
parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--maxEvents", help="Maximum number of events", type=int, default=9e99)
parser.add_argument("--maxMemory", help="Maximum job memory", type=float, default=2000)
args = parser.parse_args()

# /cvmfs/sw.hsf.org/key4hep/releases/2024-03-10/x86_64-almalinux9-gcc11.3.1-opt/k4geo/0.20-hapqru/share/k4geo/FCCee/CLD/compact/CLD_o2_v05/
# $K4GEO/FCCee/CLD/compact/CLD_o2_v05/
# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml --maxEvents 10 --cms_pool

GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10"
SINGULARITY = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/key4hep/k4-deploy/alma9\:latest"
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

class DDSimProducer:

    def __init__(self, geometry_file, storagedir, args):
        self.geometry_file = os.path.abspath(os.path.expandvars(geometry_file))
        self.geometry_dir = os.path.dirname(self.geometry_file)
        self.geometry_name = os.path.splitext(os.path.basename(self.geometry_file))[0]
        self.storagedir = storagedir
        self.cwd = os.getcwd()

        if not os.path.exists(self.geometry_file):
            logger.error(f"Geometry file {self.geometry_file} does not exist, exiting")
            sys.exit(1)
        logger.info(f"Found geometry {self.geometry_name}")

        self.local_dir = f"{self.cwd}/local_ddsim/{self.geometry_name}/"

        self.condor_queue = args.condor_queue
        self.condor_priority = args.condor_priority
        self.args = args

    def generate_gp_submit(self, input_dir):
        self.input_dir = os.path.abspath(input_dir)
        if not os.path.exists(self.input_dir):
            logger.error(f"Input directory {self.input_dir} does not exist, exiting")
            sys.exit(1)

        self.gp_name = os.path.basename(os.path.normpath(self.input_dir))
        logger.info(f"Found guineapig directory {self.gp_name}")

        logger.info(f"Get .pairs files from {self.input_dir}")
        input_files = glob.glob(f"{self.input_dir}/*.pairs")
        if len(input_files) == 0:
            logger.error(f"Input directory {self.input_dir} does not contain .pair files, exiting")
            sys.exit(1)
        logger.info(f"Found {len(input_files)} pairs files")


        # get seeds
        seeds = []
        for n, f in enumerate(input_files):
            seed = os.path.basename(f).replace("output_", "").replace(".pairs", "")
            seeds.append(seed)
            if n > args.maxEvents:
                break

        # preparing submission dirs
        self.out_dir = f"{self.storagedir}/{self.geometry_name}/{self.gp_name}/"
        self.log_dir = f"{self.storagedir}/{self.geometry_name}/{self.gp_name}/logs/"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)


        # pack the geometry files to be shipped
        logger.info("Packing geometry files")
        geometry_sandbox = f"{self.out_dir}/geometry.tar"
        os.system(f"cd {self.geometry_dir} && tar -cf {geometry_sandbox} .")

        out_dir_xrd = self.out_dir.replace("/ceph/submit", "")

        job = f"""
#!/bin/bash
SECONDS=0
unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

ls -lrt
ls -lrt /cvmfs

export seed=$1
echo $seed
hostname

# unpack geometry
mkdir geometry
tar -xf geometry.tar -C geometry

source {GP_STACK}

echo "START DDSIM"
ddsim --compactFile geometry/{self.geometry_name}.xml --outputFile output_${{seed}}_sim.root --inputFiles output_${{seed}}.pairs --numberOfEvents -1

# check ddsim exit code
# fail job in case ddsim isn't properly executed --> try again via max_retries
# on e.g. OSG, sometimes /cvmfs/sw.hsf.org/key4hep/setup.sh not found, or segfaults because of so libs not found
rc=$?
if [ $rc -ne 0 ]; then
    echo "ddsim failed with exit code $rc"
    exit $rc
fi

# check if output file exists
if [ ! -f output_${{seed}}_sim.root ]; then
    echo "ddsim did not produce output_${{seed}}_sim.root"
    exit 1
fi

echo "DONE DDSIM"
ls -lrt

echo "COPY OUTPUT FILE"
xrdcp output_${{seed}}_sim.root root://submit50.mit.edu/{out_dir_xrd}

rc=$?
if [ $rc -ne 0 ]; then
    echo "xrdcp failed with exit code $rc"
    #exit $rc
fi
echo "COPY OUTPUT DONE"


duration=$SECONDS
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
        fOut.write(f'transfer_input_files = {self.input_dir}/output_$(SEED).pairs,{geometry_sandbox}\n')
        #fOut.write(f'transfer_output_files = output_$(SEED)_sim.root\n') # done by xrdcp

        fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        fOut.write(f'max_retries    = 5\n')
        fOut.write(f'RequestMemory  = {args.maxMemory}\n')


        # global pool
        if args.cms_pool:
            fOut.write(f'+DESIRED_Sites = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"\n')

            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.log_dir}.")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = logs/x509up_u204569\n')
            fOut.write(f'+AccountingGroup      = "analysis.jaeyserm"\n')
            fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')

        # OSG pool
        elif args.osg_pool:
            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.log_dir}.")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = logs/x509up_u204569\n')
            
            # https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/requirements/#additional-feature-specific-attributes
            fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == True && HAS_SINGULARITY == TRUE )\n')



        elif 'mit.edu' in HOSTNAME:
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')
            fOut.write(f'+DESIRED_Sites = "mit_tier2,mit_tier3"\n')

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
    producer = DDSimProducer(args.geometry_file, args.storagedir, args)
    producer.generate_gp_submit(args.gp_dir)


if __name__ == "__main__":
    main()
