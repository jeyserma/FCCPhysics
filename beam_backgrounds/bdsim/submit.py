import sys, os, glob, shutil
import time
import argparse
import logging
import subprocess
import random
import socket
import json

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action='store_true', help="Submit to batch system")
parser.add_argument("-d", "--dryRun", action='store_true', help="dry run")
parser.add_argument("-c", "--clean", action='store_true', help="Clean directory")
parser.add_argument("-l", "--local", action='store_true', help="Execute locally")
parser.add_argument("-n", "--njobs", type=int, help="number of jobs", default=10)
parser.add_argument("-p", "--config", type=str, help="Config file", default="test_halo_exp")
parser.add_argument("--gmaddir", type=str, help="GMAD directory", default="GMAD_ORIG")
parser.add_argument("--suffix", type=str, help="Suffix", default="")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/work/submit/jaeyserm/fccee/FCCAnalyses/FCCPhysics/beam_backgrounds/bdsim/data/")
parser.add_argument("--logdir", type=str, help="Base directory to save the log files", default="condor")
parser.add_argument("--maxMemory", help="Maximum job memory", type=float, default=1500)
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
args = parser.parse_args()

# /work/submit/jaeyserm/fccee/FCCAnalyses/FCCPhysics/beam_backgrounds/bdsim/data/
# /ceph/submit/data/group/fcc/ee/detector/bdsim/

# python submit.py --cms_pool --submit --njobs 2600

EXE = "bdsim.py"
STACK = "/cvmfs/beam-physics.cern.ch/bdsim/x86_64-el9-gcc13-opt/bdsim-env-v1.7.7-g4v10.7.2.3-ftfp-boost.sh"
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

class BDSIMProducer:

    def __init__(self, args):
        self.args = args
        self.cwd = os.getcwd()

        self.njobs = args.njobs
        self.config = args.config
        
        self.storagedir = args.storagedir
        self.gmaddir = args.gmaddir
        
        self.suffix = f"_{args.suffix}" if args.suffix else ""
        self.logdir = f"{self.cwd}/{args.logdir}/{self.config}{self.suffix}/" 
        self.outdir = f"{self.storagedir}/{self.config}{self.suffix}/"
        
        self.configFile = f"{self.cwd}/config/{args.config}.json"
        if not os.path.isfile(self.configFile):
            logger.error(f"Config file not found: {self.configFile}")
            quit()
        
        # files to transfer 
        self.tffiles = ["GMAD_ORIG/FCC_30mm.dat", "GMAD_ORIG/fcc_ee_z_aperture.tfs", "GMAD_ORIG/fcc_ee_z_b1_twiss_end.tfs"]

    def clean(self):
        return
        # remove empty files
        os.system(f"cd {self.outdir} && find . -type f -size 0b -print")
        ans = input("Remove above zero-size files? (y/n)")
        if ans == "y":
            os.system(f"cd {self.outdir} && find . -type f -size 0b -delete")

    def generate_submit(self):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(self.logdir):
            os.makedirs(self.logdir)

        # pack GMAD input files
        logger.info("Packing GMAD input files")
        self.gmad_sandbox = f"{self.outdir}/gmad.tar"
        os.system(f"cd {self.gmaddir} && tar -cf {self.gmad_sandbox} .")
        gmad_action = "mkdir GMAD && tar -xf gmad.tar -C GMAD" # extraxt to GMAD dir

        # copy input file and executable to output dir
        os.system(f"cp {self.configFile} {self.outdir}/input.json")
        os.system(f"cp {EXE} {self.outdir}")
        #outdir_xrd = self.outdir.replace("/ceph/submit", "")

        njob = 0
        seeds = []
        while njob < self.njobs:
            seed = f"{random.randint(100000,999999)}"
            outputFile = f"{self.outdir}/output_{seed}.root"
            if os.path.exists(outputFile):
                logger.warning(f"Output file with seed {seed} already exists, skipping")
                continue
            seeds.append(seed)
            njob += 1



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



export seed=$1
export ngenerate=$2

echo "Seed"
echo $seed

echo "Ngenerate"
echo $ngenerate


echo "Checking CVMFS..."
if ! timeout 15 ls {STACK} >/dev/null 2>&1; then
    echo "ERROR: CVMFS or BDSIM stack setup.sh not found!" >&2
    exit 1
fi


echo "CVMFS looks OK. Sourcing stack"
if ! source {STACK} >/dev/null 2>&1; then
    echo "ERROR: Failed to source BDSIM stack" >&2 
    exit 1
fi
echo "BDSIM stack sourcing successful"


# unpack gmad
echo "Unpack GMAD files"
{gmad_action}



echo "Start BDSIM"
SECONDS=0
python bdsim.py --cfg input.json --seed=$seed

# check exit code
rc=$?
if [ $rc -ne 0 ]; then
    echo "bdsim failed with exit code $rc" >&2
    exit $rc
fi

duration0=$SECONDS

# check if output file exists
if [ ! -f output_$seed.root ]; then
    echo "BDSIM did not produce output_$seed.root" >&2
    exit 1
fi

# check if output file exists
if [ ! -f output_$seed.hepevt ]; then
    echo "BDSIM did not produce output_$seed.hepevt" >&2
    exit 1
fi











echo "Done BDSIM"
ls -lrt





echo "Duration BDSIM: $(($duration0))"

        """




        # make executable script
        submitFn = f"{self.outdir}/run_bdsim.sh"
        fOut = open(submitFn, "w")
        fOut.write(job)
        subprocess.getstatusoutput(f"chmod 777 {submitFn}")


        # make condor submission script
        condorFn = f'{self.outdir}/condor.cfg'
        fOut = open(condorFn, 'w')

        fOut.write(f'universe       = vanilla\n')
        fOut.write(f'initialdir     = {self.outdir}\n')
        fOut.write(f'executable     = {submitFn}\n')
        fOut.write(f'arguments      = $(SEED)\n')

        fOut.write(f'Log            = {self.logdir}/condor_job.$(ClusterId).$(ProcId).log\n')
        fOut.write(f'Output         = {self.logdir}/condor_job.$(ClusterId).$(ProcId).out\n')
        fOut.write(f'Error          = {self.logdir}/condor_job.$(ClusterId).$(ProcId).error\n')

        fOut.write(f'should_transfer_files      = YES\n')
        fOut.write(f'when_to_transfer_output    = ON_EXIT\n')
        fOut.write(f'transfer_input_files       = {os.path.basename(EXE)},input.json,gmad.tar\n')
        #fOut.write(f'transfer_output_files      = output_$(SEED).root,output_$(SEED).hepevt\n')
        fOut.write(f'transfer_output_files      = output_$(SEED).hepevt\n')

        fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        fOut.write(f'max_retries    = 10\n')
        # retry the job if it failed due to output transfer failure (HoldReasonCode == 12)
        # e.g. when job fails, it doesn't produce the output file, going to Hold (need to prevent it)
        ###fOut.write(f'periodic_release = (HoldReasonCode == 12 && NumJobStarts < 150)\n')
        ###fOut.write(f'on_exit_hold = False\n')
        fOut.write(f'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')
        #fOut.write(f'periodic_hold = (MemoryUsage > 1.15 * RequestMemory)\n') # hold when mem-requirements are not met

        #fOut.write(f'periodic_hold          = (MemoryUsage =!= UNDEFINED) && (MemoryUsage > 1.15 * RequestMemory)\n') 
        #fOut.write(f'periodic_hold_reason   = "Held: memory exceeded 1.15x RequestMemory"\n') 
        #fOut.write(f'periodic_hold_subcode  = 252\n') 
        


        
        fOut.write(f'RequestMemory  = {args.maxMemory}\n')

        proxy_path = get_voms_proxy_path()
        os.system(f"cp {proxy_path} {self.outdir}/")
        fOut.write(f'use_x509userproxy     = True\n')
        fOut.write(f'x509userproxy         = {self.outdir}/{os.path.basename(proxy_path)}\n')


        #elif 'cern.ch' in HOSTNAME:
        #    fOut.write(f'+JobFlavour    = "{self.condor_queue}"\n')
        #    fOut.write(f'+AccountingGroup = "{self.condor_priority}"\n')

        # OSG pool
        if args.osg_pool:
            # https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/requirements/#additional-feature-specific-attributes
            #fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            ##fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_CVMFS_sft_cern_ch == TRUE && HAS_CVMFS_sw_hsf_org == TRUE && HAS_SINGULARITY == TRUE )\n')
            #fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_SINGULARITY == TRUE &&  (GLIDEIN_Site == "Wisconsin" || GLIDEIN_Site == "UChicago")  )\n')
        elif args.cms_pool:
            # https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/requirements/#additional-feature-specific-attributes
            #fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            ##fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
            fOut.write(f'+DESIRED_Sites = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'+AccountingGroup      = "analysis.jaeyserm"\n')
            
            fOut.write(f'Requirements          = (  HAS_SINGULARITY == TRUE )\n')
            #fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_SINGULARITY == TRUE &&  (GLIDEIN_Site == "Wisconsin" || GLIDEIN_Site == "UChicago")  )\n')
        else:
            fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')
            fOut.write(f'+DESIRED_Sites = "mit_tier2,mit_tier3"\n')
            #fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            #fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            #fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            #fOut.write(f'+SingularityBindCVMFS   = True\n')
            #fOut.write(f'Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')



        seedsStr = ' \n '.join([str(s) for s in seeds])
        fOut.write(f'queue SEED in ( \n {seedsStr} \n)\n')

        fOut.close()



        subprocess.getstatusoutput(f'chmod 777 {condorFn}')
        os.system(f"condor_submit {condorFn}")

        logger.info(f"Written to {self.outdir}")



    def dryRun(self):
        # check environment ?
        rundir = f"/tmp/bdsim/{self.config}{self.suffix}/"
        subprocess.run(["rm", "-rf", rundir], check=True)
        subprocess.run(["mkdir", "-p", rundir], check=True)
        subprocess.run(["mkdir", "GMAD"], cwd=rundir, check=True)
        for f in self.tffiles:
            subprocess.run(["cp", f"{self.cwd}/{f}", "GMAD"], cwd=rundir, check=True)
        subprocess.run(["cp", self.configFile, "input.json"], cwd=rundir, check=True)
        subprocess.run(["python", f"{self.cwd}/{EXE}", "--cfg", "input.json"], cwd=rundir, check=True)

def main():
    producer = BDSIMProducer(args)
    if args.submit:
        producer.generate_submit()
    if args.clean:
        producer.clean()
    if args.dryRun:
        producer.dryRun()


if __name__ == "__main__":
    main()
