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
parser.add_argument("-c", "--clean", action='store_true', help="Clean directory")
parser.add_argument("-l", "--local", action='store_true', help="Execute locally")
parser.add_argument("-n", "--njobs", type=int, help="number of jobs", default=10)
parser.add_argument("-i", "--inputfile", type=str, help="Input file", default="cards/fccee_v23.dat")
parser.add_argument("-p", "--accelerator", type=str, help="Accelerator config", default="FCCee_Z_4IP_04may23")
parser.add_argument("-d", "--parameter_set", type=str, help="Parameter set", default="FCCee_Z256")
parser.add_argument("--suffix", type=str, help="Suffix", default="")
parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="workday")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/ee/detector/guineapig_fixWindow/") # guineapig_bz_2T  guineapig
parser.add_argument("--maxMemory", help="Maximum job memory", type=float, default=2000)
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
args = parser.parse_args()



GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10"
#GP_EXEC = f"{os.getcwd()}/guinea-pig/gp/bin/guinea"
#GP_EXEC = f"{os.getcwd()}/guinea-pig_bz/gp/bin/guinea"
#GP_EXEC = f"{os.getcwd()}/guinea-pig_tracking/gp/bin/guinea"
GP_EXEC = f"{os.getcwd()}/guinea-pig_tracking/gp/bin/guinea"
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

class GPProducer:

    def __init__(self, accelerator, parameter_set, inputFile, storagedir, args):
        self.accelerator = accelerator
        self.parameter_set = parameter_set
        self.inputFile = inputFile
        self.storagedir = storagedir

        self.cwd = os.getcwd()
        
        if self.inputFile[0] != "/":
            self.inputFile = f"{self.cwd}/{self.inputFile}"
        self.baseInputFile = os.path.basename(self.inputFile)

        self.suffix = f"_{args.suffix}" if args.suffix else ""

        self.out_dir = f"{self.storagedir}/{self.accelerator}_{self.parameter_set}{self.suffix}/"
        self.log_dir = f"{self.out_dir}/logs/"
        self.local_dir = f"{self.cwd}/local_gp/{self.accelerator}_{self.parameter_set}/"

        self.condor_queue = args.condor_queue
        self.condor_priority = args.condor_priority
        self.args = args

    def clean(self):
        # remove empty files
        os.system(f"cd {self.out_dir} && find . -type f -size 0b -print")
        ans = input("Remove above zero-size files? (y/n)")
        if ans == "y":
            os.system(f"cd {self.out_dir} && find . -type f -size 0b -delete")

    def generate_submit(self, njobs):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        # copy input file and executable to output dir
        os.system(f"cp {self.inputFile} {self.out_dir}")
        os.system(f"cp {GP_EXEC} {self.out_dir}")
        out_dir_xrd = self.out_dir.replace("/ceph/submit", "")

        njob = 0
        seeds = []
        while njob < njobs:
            seed = f"{random.randint(100000,999999)}"
            outputFile = f"{self.out_dir}/output_{seed}.pairs"
            if os.path.exists(outputFile):
                logger.warning(f"Output file with seed {seed} already exists, skipping")
                continue
            seeds.append(seed)
            njob += 1



        job = f"""
#!/bin/bash

set -e

SECONDS=0
unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

echo "Release:"
cat /proc/version

echo "Hostname:"
hostname

echo "List current working dir:"
ls -lrt

echo "List cvmfs dir"
ls -lrt /cvmfs


echo "Set seed"
export seed=$1
echo $seed
sed -i -e \"s/rndm_seed=100000/rndm_seed=$seed/g\" {self.baseInputFile}




echo "Checking CVMFS..." >&2

if ! ls /cvmfs/sw.hsf.org/key4hep/setup.sh >/dev/null 2>&1; then
    echo "ERROR: CVMFS or key4hep setup.sh not found!" >&2
    exit 1
fi

echo "CVMFS looks OK. Sourcing..." >&2

if ! source {GP_STACK} >/dev/null 2>&1; then
    echo "ERROR: Failed to source Key4Hep setup" >&2
    exit 1
fi

echo "Key4Hep sourcing successful" >&2


echo "Start guinea-pig"
mv {os.path.basename(GP_EXEC)} run
chmod 777 run
./run --acc_file {self.baseInputFile} {self.accelerator} {self.parameter_set} output

# check exit code
# fail job in case guinea-pig isn't properly executed --> try again via max_retries
# on e.g. OSG, sometimes segfaults because of shared object libs not found
rc=$?
if [ $rc -ne 0 ]; then
    echo "guinea-pig failed with exit code $rc" >&2
    exit $rc
fi

# check if output file exists
if [ ! -f pairs.dat ]; then
    echo "guinea-pig did not produce pairs.dat" >&2
    exit 1
fi

# check if output file exists
if [ ! -f pairs0.dat ]; then
    echo "guinea-pig did not produce pairs0.dat" >&2
    exit 1
fi

echo "Done guinea-pig"

mv pairs.dat output_$1.pairs
mv pairs0.dat output0_$1.pairs

ls -lrt


#echo "Setup certificates directories"
#ls -lrt /etc/grid-security/

#if [ -d /cvmfs/grid.cern.ch/etc/grid-security/certificates ]; then
#    export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
#else
#    echo "ERROR: Grid certs not available on CVMFS"
#    exit 1
#fi



#echo $X509_USER_PROXY

#source /cvmfs/grid.cern.ch/alma9-ui-current/etc/profile.d/setup-alma9-test.sh
#echo "Copy output file"
#echo $X509_USER_PROXY
#export X509_USER_PROXY=x509
#echo $X509_USER_PROXY

#voms-proxy-info -all
#voms-proxy-info -all -file x509

#xrdcp --version
#xrdcp -d 3 output_${{seed}}_sim.root root://submit50.mit.edu/{out_dir_xrd}/output_${{seed}}_sim.root
#rc=$?
#if [ $rc -ne 0 ]; then
#    echo "xrdcp failed with exit code $rc" >&2
#    #exit $rc
#fi
#echo "Copy output file done"


duration=$SECONDS
echo "Duration: $(($duration))"

        """

        # make executable script
        submitFn = f"{self.out_dir}/run_gp.sh"
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
        fOut.write(f'transfer_input_files = {os.path.basename(GP_EXEC)},{self.inputFile}\n')
        fOut.write(f'transfer_output_files = output_$(SEED).pairs,output0_$(SEED).pairs\n')

        fOut.write(f'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        fOut.write(f'max_retries    = 50\n')
        # retry the job if it failed due to output transfer failure (HoldReasonCode == 12)
        # e.g. when job fails, it doesn't produce the output file, going to Hold (need to prevent it)
        fOut.write(f'periodic_release = (HoldReasonCode == 12 && NumJobStarts < 50)\n')
        fOut.write(f'on_exit_hold = False\n')
        fOut.write(f'RequestMemory  = {args.maxMemory}\n')


        if 'mit.edu' in HOSTNAME:
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

        # OSG pool
        elif args.osg_pool:
            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.log_dir}/x509")
            fOut.write(f'use_x509userproxy     = True\n')
            #fOut.write(f'x509userproxy         = /tmp/x509up_u$(id -u)\n')
            fOut.write(f'x509userproxy         = logs/x509\n')
            
            # https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/requirements/#additional-feature-specific-attributes
            #fOut.write(f'+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            ##fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            fOut.write(f'+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
            fOut.write(f'+ProjectName            = "MIT_submit"\n')
            fOut.write(f'+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_SINGULARITY == TRUE )\n')
            # HAS_CVMFS_sw_hsf_org
            # HAS_CVMFS_singularity_opensciencegrid_org
            # && HAS_CVMFS_sw_hsf_org == True





        seedsStr = ' \n '.join([str(s) for s in seeds])
        fOut.write(f'queue SEED in ( \n {seedsStr} \n)\n')

        fOut.close()


        subprocess.getstatusoutput(f'chmod 777 {condorFn}')
        os.system(f"condor_submit {condorFn}")

        logger.info(f"Written to {self.out_dir}")

    def generate_local(self):
        if os.path.exists(self.local_dir):
            logger.error(f"Please remove local directory {self.local_dir}")
            sys.exit(1)
        os.makedirs(self.local_dir)
        logger.info(f"Generate GP event locally")

        submitFn = f"{self.local_dir}/run_gp.sh"

        fOut = open(submitFn, "w")
        fOut.write("#!/bin/bash\n")
        fOut.write("SECONDS=0\n")
        fOut.write("unset LD_LIBRARY_PATH\n")
        fOut.write("unset PYTHONHOME\n")
        fOut.write("unset PYTHONPATH\n")
        fOut.write(f"source {GP_STACK}\n")
        fOut.write(f"cp {self.inputFile} {self.baseInputFile}\n")
        fOut.write(f"ls -lrt\n")

        fOut.write('echo "START GUINEA-PIG"\n')
        fOut.write(f"{GP_EXEC} --acc_file {self.baseInputFile} {self.accelerator} {self.parameter_set} output \n")
        fOut.write("echo \"DONE GUINEA-PIG\"\n")
        fOut.write("duration=$SECONDS\n")
        fOut.write("echo \"Duration: $(($duration)) seconds\"\n")
        fOut.close()

        subprocess.getstatusoutput(f"chmod 777 {submitFn}")

        os.system(f"cd {self.local_dir} && ./run_gp.sh | tee output.txt")
        logger.info(f"Done, saved to {self.local_dir}")




def main():
    producer = GPProducer(args.accelerator, args.parameter_set, args.inputfile, args.storagedir, args)
    if args.submit:
        producer.generate_submit(njobs=args.njobs)
    if args.clean:
        producer.clean()
    if args.local:
        producer.generate_local()



if __name__ == "__main__":
    main()
