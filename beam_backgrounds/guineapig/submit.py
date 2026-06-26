import sys, os, glob, shutil, re
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
parser.add_argument("--submit", action='store_true', help="Submit to batch system")
parser.add_argument("--merge", action='store_true', help="Merge ROOT files")
parser.add_argument("--dryrun", action='store_true', help="dry run")
parser.add_argument("--suffix", type=str, help="Suffix", default="")

parser.add_argument("-n", "--njobs", type=int, help="number of jobs", default=10)
parser.add_argument("-i", "--input_file", type=str, help="Input file", default="cards/studies.dat")
parser.add_argument("-p", "--accelerator", type=str, help="Accelerator config", default="FCCee_Z_GHC_V25p1")
parser.add_argument("-d", "--parameter_set", type=str, help="Parameter set", default="Z256_2T_grids8")

parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="longlunch")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/ee/detector/guineapig_studies/")
parser.add_argument("--logdir", type=str, help="Base directory to save the log files", default="logdir")

parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--max_memory", help="Maximum job memory", type=float, default=500)
parser.add_argument("--njobs_per_sub", help="Maximum number of jobs per submission", type=int, default=5000)
parser.add_argument("--xrootd", action='store_true', help="Use XrootD transfer")

parser.add_argument("--seeddir", type=str, help="Base directory to extract seeds", default=None)

args = parser.parse_args()



SINGULARITY = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"
SINGULARITY = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/key4hep/k4-deploy/alma9:latest"

GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29"
#GP_EXEC = f"{os.getcwd()}/guinea-pig-15122025/gp/bin/guinea"
GP_EXEC = f"{os.getcwd()}/guinea-pig-15122025_dev/build/src/guinea"
GP_EXEC = f"{os.getcwd()}/guinea-pig-15122025_dev_time/build/src/guinea"
#GP_EXEC = f"{os.getcwd()}/guinea-pig-15122025_dev/build_theta0/src/guinea"

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


pdg_map = {
    13: 'mu-',
    -13: 'mu+',
    11: 'e-',
    -11: 'e+',
    22: 'gamma',
    211: 'pi+',
    -211: 'pi-',
}




class GPProducer:

    def __init__(self, args):
        self.args = args
        self.cwd = os.getcwd()
        self.stack = GP_STACK
        self.gp_exec = GP_EXEC
        self.storagedir = args.storagedir


        self.accelerator = args.accelerator
        self.parameter_set = args.parameter_set
        self.input_file = args.input_file
        if self.input_file[0] != "/":
            self.input_file = f"{self.cwd}/{self.input_file}"
        #self.input_file_name = os.path.basename(self.input_file)

        self.suffix = f"_{args.suffix}" if args.suffix else ""
        self.outdir = f"{self.storagedir}/{self.accelerator}/{self.parameter_set}{self.suffix}/" # main output dir
        self.cfgdir = f"{self.outdir}/cfg/" # storage of executable, input file, seeds, ...
        self.stodir = f"{self.outdir}/unmerged/" # output files
        self.logdir = f"{self.cwd}/{args.logdir}/{self.accelerator}/{self.parameter_set}{self.suffix}/"


    def make_script(self):
        return f"""#!/bin/bash
                set -euo pipefail
                if [ $# -lt 1 ]; then
                    echo "Usage: $0 SEED" >&2
                    exit 1
                fi
                seed="$1"
                echo "Release:"
                cat /proc/version
                echo "Hostname:"
                hostname
                echo "Current working directory:"
                pwd

                echo "List current working dir"
                ls -lrt
                set +u
                echo "Checking CVMFS"
                if ! timeout 15 ls "{self.stack.split(' ')[0]}" >/dev/null 2>&1; then
                    echo "ERROR: Stack not found or not accessible: {self.stack.split(' ')[0]}" >&2
                    exit 1
                fi
                echo "Sourcing stack"
                if ! source {self.stack} >/dev/null 2>&1; then
                    echo "ERROR: Failed to source stack: {self.stack}" >&2
                    exit 1
                fi
                set -u

                echo "Set seed $seed"
                sed -i -e \"s/rndm_seed=100000/rndm_seed=$seed/g\" input.dat

                ls -lrt
                SECONDS=0
                echo "Running guinea-pig"
                chmod 777 run_guinea
                export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/fftw3/3.3.10-33229/x86_64-el9-gcc11-opt/lib:$LD_LIBRARY_PATH
                ./run_guinea --acc_file input.dat {self.accelerator} {self.parameter_set} output
                ls -lrt
                if [ ! -s "output.root" ]; then
                    echo "ERROR: expected output file output.root not found or empty" >&2
                    exit 100
                fi
                
                mv pairs.dat output_$seed.pairs
                mv pairs0.dat output0_$seed.pairs
                mv output output_$seed.log
                mv output.root output_$seed.root
                ls -lrt

                echo "guinea-pig step done"
                gen_duration=$SECONDS

                ls -lrt
                echo "Done script, total duration ${{gen_duration}} seconds"
                echo "${{gen_duration}}"

                """

    def generate(self):

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            os.makedirs(self.cfgdir)
            os.makedirs(self.stodir)

            # copy if does not exist!!
            if not os.path.isfile(f"{self.cfgdir}/input.dat"):
                os.system(f"cp {self.input_file} {self.cfgdir}/input.dat")
            if not os.path.isfile(f"{self.cfgdir}/run_guinea"):
                os.system(f"cp {self.gp_exec} {self.cfgdir}/run_guinea") # don't use guinea, can clash with default stack version
            os.system(f"touch {self.cfgdir}/seeds.txt")

        self.transfer_input_files = ['input.dat', 'run_guinea']
        self.transfer_output_files  = []
        self.transfer_output_files.append(f"output_$(SEED).pairs")
        self.transfer_output_files.append(f"output0_$(SEED).pairs")
        self.transfer_output_files.append(f"output_$(SEED).log")
        self.transfer_output_files.append(f"output_$(SEED).root")

        self.max_memory = args.max_memory

        self.used_seeds = set()
        with open(f"{self.cfgdir}/seeds.txt", "r") as f:
            for line in f:
                seed = line.strip()
                self.used_seeds.add(seed)

        njobs = 0
        self.seeds = []
        if self.args.seeddir == None:
            while njobs < args.njobs:
                seed = f"{random.randint(100000,999999)}"
                if seed in self.used_seeds:
                    logger.warning(f"Seed {seed} already exists, skipping")
                    continue
                self.seeds.append(seed)
                njobs += 1
        else:
            logger.info(f"Extract seeds from directory {self.args.seeddir}")

            for f in glob.glob(f"{self.args.seeddir}/output_*.root"):
                seed = re.search(r"output_(\d+)\.root$", f).group(1)
                self.seeds.append(seed)
                
            logger.info(f"Found {len(self.seeds)} seeds") 

        for seed in self.seeds:
            os.system(f"echo {seed} >> {self.cfgdir}/seeds.txt")

        # make executable script
        script_sandbox = self.make_script()
        submitFn = f"{self.cfgdir}/run_guinea.sh"
        fOut = open(submitFn, "w")
        fOut.write(script_sandbox)

        subprocess.getstatusoutput(f"chmod 777 {submitFn}")

        # get latest submission
        submf = glob.glob(f'{self.cfgdir}/condor_v*.cfg')
        submv = versions = [int(re.search(r"_v(\d+)\.cfg$", f).group(1)) for f in submf if re.search(r"_v(\d+)\.cfg$", f)]
        max_submv = max(submv) if submv else 0


        seeds_chunked = chunk_list(self.seeds, args.njobs_per_sub)
        for i,seeds in enumerate(seeds_chunked):
            subv = i+1+max_submv
            logger.info(f"Submit {i+1}/{len(seeds_chunked)} with version {subv} ")

            logdir = f"{self.logdir}/v{subv}/"
            if not os.path.exists(logdir):
                os.makedirs(logdir)


            # make condor submission script
            condorFn = f'{self.cfgdir}/condor_v{subv}.cfg'
            fOut = open(condorFn, 'w')

            fOut.write(f'universe           = vanilla\n')
            fOut.write(f'initialdir         = {self.cfgdir}\n')
            fOut.write(f'output_directory   = {self.stodir}\n')
            
            fOut.write(f'executable         = {submitFn}\n')
            fOut.write(f'arguments          = $(SEED)\n')

            fOut.write(f'Log                = {logdir}/condor_job.$(ClusterId).$(ProcId).log\n')
            fOut.write(f'Output             = {logdir}/condor_job.$(ClusterId).$(ProcId).out\n')
            fOut.write(f'Error              = {logdir}/condor_job.$(ClusterId).$(ProcId).error\n')

            fOut.write(f'should_transfer_files = YES\n')
            fOut.write(f'when_to_transfer_output = ON_EXIT\n')

            fOut.write(f'transfer_input_files = {",".join(self.transfer_input_files)}\n')
            fOut.write(f'transfer_output_files = {",".join(self.transfer_output_files)}\n') # done by xrdcp
            #if args.xrootd:
            #    fOut.write(f'transfer_output_files = ""\n') # explicit no transfer files back

                

            fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
            fOut.write(f'max_retries    = 3\n')
            fOut.write(f'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')
            #fOut.write(f'periodic_release = (NumJobStarts < 3)\n')



            # Intercept memory growth before the site removes the job at 1.2 * RequestMemory
            fOut.write(f'periodic_hold = (JobStatus == 2) && (MemoryUsage > 1.10 * RequestMemory)\n')
            fOut.write(f'periodic_hold_reason = "Retrying after high memory usage"\n')
            fOut.write(f'periodic_hold_subcode = 9001\n')
            fOut.write(f'periodic_release = (((HoldReasonCode == 12) && (HoldReasonSubCode == 2) && (NumHolds < 3)) || ((HoldReasonCode == 3) && (HoldReasonSubCode == 9001) && (NumHolds < 3)) || ((HoldReasonCode == 3) && (HoldReasonSubCode == 0) && (NumHolds < 3)) )\n')
            


            # Only auto-release the holds we created ourselves
            
            fOut.write(f'+JobBatchName = "GUINEA_{self.accelerator}_{self.parameter_set}{self.suffix}_v{subv}"\n')

            fOut.write(f'RequestMemory  = {self.max_memory}\n')

            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.cfgdir}/")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = {self.cfgdir}/{os.path.basename(proxy_path)}\n')

            # global pool
            if args.cms_pool:
                #proxy_path = get_voms_proxy_path()
                #os.system(f"cp {proxy_path} {self.log_dir}/")
                #fOut.write(f'use_x509userproxy     = True\n')
                #fOut.write(f'x509userproxy         = logs/{os.path.basename(proxy_path)}\n')

                fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
                fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
                fOut.write(f'+DESIRED_Sites = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"\n')
                fOut.write(f'+SingularityBindCVMFS   = True\n')
                fOut.write(f'+AccountingGroup      = "analysis.jaeyserm"\n')
                fOut.write(f'Requirements          = (  HAS_SINGULARITY == TRUE )\n')

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

            logger.info(f"Written to {self.outdir}")

    def dryrun(self):
        rundir = f"/tmp/guineapig/{self.accelerator}/{self.parameter_set}{self.suffix}/"

        script_init = f"""
        set -e

        rm -rf {rundir}
        mkdir -p {rundir}
        cd {rundir}
        pwd
        cp {self.gp_exec} run_guinea
        cp {self.input_file} input.dat
        ls -lrt

        """
        subprocess.run(["/bin/bash", "-c", script_init])

        script_sandbox = self.make_script()
        with open(f"{rundir}/run.sh", "w") as tf:
            tf.write(script_sandbox)
        subprocess.run(["/bin/bash", "run.sh", "123456"], cwd=rundir, env={})

    def merge(self):
        batch_size = 200

        #mergedir = os.path.join(self.outdir, "merged")
        #os.makedirs(mergedir, exist_ok=True)

        files = sorted(glob.glob(f"{self.stodir}/*.root", recursive=True))

        print(f"Found {len(files)} ROOT files")

        # split into chunks of 1000
        for i in range(0, len(files), batch_size):

            batch = files[i:i + batch_size]

            outfile = os.path.join(
                self.outdir,
                f"output_{i//batch_size:04d}.root"
            )

            cmd = [
                "hadd",
                "-f",
                "-j", "16",
                outfile
            ] + batch

            print(f"Merging {len(batch)} files -> {outfile}")
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

def main():
    producer = GPProducer(args)
    if args.submit:
        producer.generate()
    #if args.clean:
    #    producer.clean()
    if args.dryrun:
        producer.dryrun()
    if args.merge:
        producer.merge()

if __name__ == "__main__":
    main()
