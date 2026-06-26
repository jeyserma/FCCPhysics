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
parser.add_argument("--submit", action='store_true', help="Submit to batch system")
parser.add_argument("--dryrun", action='store_true', help="dry run")

parser.add_argument("--njobs", type=int, help="number of jobs", default=10)
parser.add_argument("--ngenerate", type=int, help="number of events", default=200000)
parser.add_argument("--lattice", type=str, help="Lattice", default="LCC_v1")
parser.add_argument("--runconfig", type=str, help="Config run file", default="halo")

parser.add_argument("--suffix", type=str, help="Suffix", default="")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim")
parser.add_argument("--logdir", type=str, help="Base directory to save the log files", default="logdir")
parser.add_argument("--max_memory", help="Maximum job memory", type=float, default=2000)
parser.add_argument("--njobs_per_sub", help="Maximum number of jobs per submission", type=int, default=5000)
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
args = parser.parse_args()

# /work/submit/jaeyserm/fccee/FCCAnalyses/FCCPhysics/beam_backgrounds/bdsim/data/
# /ceph/submit/data/group/fcc/ee/detector/bdsim/

# python submit.py --cms_pool --submit --njobs 2600


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

def chunk_list(lst, chunk_size):
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]



class BDSIMProducer:

    def __init__(self, args):
        self.args = args
        self.cwd = os.getcwd()

        self.njobs = args.njobs
        self.ngenerate = args.ngenerate

        self.lattice = args.lattice
        self.runconfig = args.runconfig
        self.configdir = f"{self.cwd}/config/{self.lattice}/"
        
        self.storagedir = args.storagedir
        self.max_memory = args.max_memory  
        
        self.suffix = f"_{args.suffix}" if args.suffix else ""
        self.logdir = f"{self.cwd}/{args.logdir}/{self.lattice}_{self.runconfig}{self.suffix}/" 
        self.outdir = f"{self.storagedir}/{self.lattice}_{self.runconfig}{self.suffix}/"

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # pack sandbox
        logger.info("Creating sandbox")
        self.sandbox = f"{self.outdir}/sandbox.tar"
        os.system(f"tar -cvf {self.sandbox} -C {self.configdir} .")

        self.transfer_input_files = ['sandbox.tar']
        self.transfer_output_files  = ['output_$(SEED).hepevt']

        njob = 0
        self.seeds = []
        while njob < self.njobs:
            seed = f"{random.randint(100000,999999)}"
            outputFile = f"{self.outdir}/output_{seed}.hepevt"
            if os.path.exists(outputFile):
                logger.warning(f"Output file with seed {seed} already exists, skipping")
                continue
            self.seeds.append(seed)
            njob += 1



    def make_script(self, runconfig):
        return f"""
            #!/bin/bash
            set -euo pipefail

            if [ $# -lt 2 ]; then
                echo "Usage: $0 SEED NGENERATE" >&2
                exit 1
            fi

            seed="$1"
            ngenerate="$2"   # total number of events
            echo $seed
            echo $ngenerate

            echo "Release:"
            cat /proc/version

            echo "Hostname:"
            hostname

            echo "Current working directory:"
            pwd

            echo "List current working dir"
            ls -lrt

            echo "Checking CVMFS"
            set +u
            if ! timeout 15 ls "{STACK}" >/dev/null 2>&1; then
                echo "ERROR: Stack not found or not accessible: {STACK}" >&2
                exit 1
            fi

            echo "Sourcing stack"
            if ! source "{STACK}" >/dev/null 2>&1; then
                echo "ERROR: Failed to source stack: {STACK}" >&2
                exit 1
            fi
            set -u

            echo "Unpack sandbox"
            if ! tar -xf sandbox.tar; then
                echo "ERROR: Failed to unpack sandbox.tar" >&2
                exit 1
            fi
            ls -lrt

            SECONDS=0

            echo "Running generation"
            python "run_{runconfig}.py" --seed="$seed" --ngenerate="$ngenerate"

            if [ ! -f "output_${{seed}}.root" ]; then
                echo "ERROR: Generation did not produce output_${{seed}}.root" >&2
                exit 1
            fi

            echo "Generation step done"
            gen_duration=$SECONDS

            echo "Running conversion"
            python convert.py --input "output_${{seed}}.root"

            if [ ! -f "output_${{seed}}.hepevt" ]; then
                echo "ERROR: Conversion did not produce output_${{seed}}.hepevt" >&2
                exit 1
            fi

            ls -lrt

            total_duration=$SECONDS
            echo "Done script, generation duration ${{gen_duration}} seconds, total duration ${{total_duration}} seconds"

        """



    def generate_submit(self):

        # make executable script
        script_sandbox = self.make_script(self.runconfig)
        submitFn = f"{self.outdir}/run_bdsim.sh"
        fOut = open(submitFn, "w")
        fOut.write(script_sandbox)
        subprocess.getstatusoutput(f"chmod 777 {submitFn}")

        seeds_chunked = chunk_list(self.seeds, args.njobs_per_sub)
        for i,seeds in enumerate(seeds_chunked):
            subv = i+1
            logger.info(f"Submit {subv}/{len(seeds_chunked)} ")

            logdir = f"{self.logdir}/v{subv}/"
            if not os.path.exists(logdir):
                os.makedirs(logdir)

            # make condor submission script
            condorFn = f'{self.outdir}/condor_v{subv}.cfg'
            fOut = open(condorFn, 'w')

            fOut.write(f'universe       = vanilla\n')
            fOut.write(f'initialdir     = {self.outdir}\n')
            fOut.write(f'executable     = {submitFn}\n')
            fOut.write(f'arguments      = $(SEED) {self.ngenerate}\n')

            fOut.write(f'Log            = {logdir}/condor_job.$(ClusterId).$(ProcId).log\n')
            fOut.write(f'Output         = {logdir}/condor_job.$(ClusterId).$(ProcId).out\n')
            fOut.write(f'Error          = {logdir}/condor_job.$(ClusterId).$(ProcId).error\n')

            fOut.write(f'should_transfer_files = YES\n')
            fOut.write(f'when_to_transfer_output = ON_EXIT\n')

            fOut.write(f'transfer_input_files = {",".join(self.transfer_input_files)}\n')
            fOut.write(f'transfer_output_files = {",".join(self.transfer_output_files)}\n') # done by xrdcp

            fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
            fOut.write(f'max_retries    = 3\n')
            fOut.write(f'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')

            # Intercept memory growth before the site removes the job at 1.2 * RequestMemory
            fOut.write(f'periodic_hold = (JobStatus == 2) && (MemoryUsage > 1.10 * RequestMemory)\n')
            fOut.write(f'periodic_hold_reason = "Retrying after high memory usage"\n')
            fOut.write(f'periodic_hold_subcode = 9001\n')
            fOut.write(f'periodic_release = (((HoldReasonCode == 12) && (HoldReasonSubCode == 2) && (NumHolds < 3)) || ((HoldReasonCode == 3) && (HoldReasonSubCode == 9001) && (NumHolds < 3)) || ((HoldReasonCode == 3) && (HoldReasonSubCode == 0) && (NumHolds < 3)) )\n')
            

        
            fOut.write(f'+JobBatchName = "BDSIM_{self.lattice}_{self.runconfig}{self.suffix}_v{subv}"\n')
            fOut.write(f'RequestMemory  = {self.max_memory}\n')

            

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



    def dryrun(self):
        rundir = f"/tmp/bdsim/{self.lattice}/{self.runconfig}{self.suffix}/"

        script_init = f"""
        set -e

        rm -rf {rundir}
        mkdir -p {rundir}
        cd {rundir}
        pwd

        cp {self.sandbox} .
        ls -lrt

        """
        subprocess.run(["/bin/bash", "-c", script_init])

        script_sandbox = self.make_script(self.runconfig)
        with open(f"{rundir}/run.sh", "w") as tf:
            tf.write(script_sandbox)
        subprocess.run(["bash", "run.sh", "12345", f"{self.args.ngenerate}"], cwd=rundir)


def main():
    producer = BDSIMProducer(args)
    if args.submit:
        producer.generate_submit()
    if args.dryrun:
        producer.dryrun()


if __name__ == "__main__":
    main()
