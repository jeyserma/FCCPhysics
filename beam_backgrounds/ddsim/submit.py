import sys, os, glob, shutil
import time
import argparse
import logging
import subprocess
import random
import socket
import geometry

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--submit", action='store_true', help="Submit to batch system")
parser.add_argument("--dryrun", action='store_true', help="dry run")
parser.add_argument("--suffix", type=str, help="Suffix", default="")

parser.add_argument("--geometry", type=str, help="Geometry", default="IDEA_o1_v03_MDI") # IDEA_o1_v03_MDI CLD_o2_v07

parser.add_argument("--input_source", type=str, help="Input source (depending on the type)")
parser.add_argument("--input_name", type=str, help="Input name, will define output storage dir", default="test6")
parser.add_argument("--input_type", type=str, help="Source of geometry", choices=["pairs", "gun", "hepevt", "root"], default="pairs")

parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="longlunch")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/ee/detector/ddsim/samples")
parser.add_argument("--logdir", type=str, help="Base directory to save the log files", default="logdir")

parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--max_jobs", help="Maximum number of jobs", type=int, default=-1)
parser.add_argument("--njobs_per_sub", help="Maximum number of jobs per submission", type=int, default=10000)
parser.add_argument("--max_memory", help="Maximum job memory", type=float, default=-1)
parser.add_argument("--xrootd", action='store_true', help="Use XrootD transfer")
args = parser.parse_args()



SINGULARITY = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"
SINGULARITY = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/key4hep/k4-deploy/alma9:latest"

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




class DDSimProducer:

    def __init__(self, args):
        self.args = args
        self.cwd = os.getcwd()
        self.input_name = args.input_name
        self.geometry_name = args.geometry

        self.storagedir = args.storagedir
        
        self.suffix = f"_{args.suffix}" if args.suffix else ""
        self.logdir = f"{self.cwd}/{args.logdir}/{self.geometry_name}/{self.input_name}/{self.suffix}/" 
        self.outdir = f"{self.storagedir}/{self.geometry_name}/{self.input_name}{self.suffix}/"

        self.ddsim_xargs = ""

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        

        # parse geometry
        self.geodict = geometry.geodict
        if not self.geometry_name in self.geodict:
            logging.error(f"Geometry {self.geometry_name} not found")
            quit()
        self.geocfg = self.geodict[self.geometry_name]
        self.stack = self.geocfg['stack']
        self.steering_file = f"{self.cwd}/{self.geocfg['steering_file']}"
        self.geometry_sandbox = ""
        self.geometry_action = ""
        self.compact_file = self.geocfg['compact_file']
        

        self.transfer_input_files = ['steering.py']
        self.transfer_output_files  = []

        self.max_memory = self.geocfg['max_memory']
        if args.max_memory > 0:
            self.max_memory = args.max_memory           

        if self.geocfg['source'] == "stack":
            pass
            #self.compactfile = self.geocfg['compactfile']
        elif self.geocfg['source'] == "local":
            
            self.geometry_source = self.geocfg['geometry_source']

            logger.info("Packing local geometry files")
            self.geometry_sandbox = f"{self.outdir}/geometry.tar"
            self.transfer_input_files.append('geometry.tar')
            os.system(f"tar -cf {self.geometry_sandbox} -C {self.cwd}/{self.geometry_source}/ .")
            
            if args.dryrun:
                self.geometry_action = f"cp {self.geometry_sandbox} . \n"
            self.geometry_action += f"""
                mkdir k4geo
                tar -xf geometry.tar -C k4geo
                cd k4geo
                mkdir build install
                cd build
                cmake .. -DCMAKE_INSTALL_PREFIX=../install -D INSTALL_BEAMPIPE_STL_FILES=ON
                make install -j {8 if args.dryrun else 1}
                cd ..
                set +u
                k4_local_repo
                set -u
                cd ..
            """    
        else:
            logging.error("Not implemented")
            quit()
        logger.info(f"Found geometry {self.geometry_name} from stack")

        # parse input source 
        self.input_type = args.input_type
        self.input_source = args.input_source
        
        os.system(f"cp {self.steering_file} {self.outdir}/steering.py")

        
        self.seeds = []
        if self.input_type == "pairs":
            self.ddsim_xargs += f" --inputFiles ${{seed}}.pairs --outputFile ${{seed}}.root --numberOfEvents -1"

            input_files = glob.glob(f"{self.input_source}/output_*.pairs")
            if len(input_files) == 0:
                logger.error(f"Input source {self.input_source} does not contain pair files, exiting")
                sys.exit(1)
            logger.info(f"Found {len(input_files)} pairs files")

            # get files and seeds
            if args.dryrun:
                self.dryrun_file = input_files[0]
                self.dryrun_seed = os.path.basename(self.dryrun_file).replace("output_", "").replace(".pairs", "")
            else:
                njobs = 0
                for n, f in enumerate(input_files):
                    if args.max_jobs > 0 and njobs >= args.max_jobs:
                        break
                    njobs += 1
                    seed = os.path.basename(f).replace(".pairs", "")
                    if os.path.exists(f"{self.outdir}/{seed}.root"):
                        logger.info(f"File exists, skip")  
                    else:
                        logger.info(f"Add seed {seed}")
                        self.seeds.append(seed)
                self.transfer_input_files.append(f"{self.input_source}/$(SEED).pairs")
                self.transfer_output_files.append(f"$(SEED).root")

        elif self.input_type == "hepevt":
            self.ddsim_xargs += f" --inputFiles ${{seed}}.hepevt --outputFile ${{seed}}.root --numberOfEvents -1"

            input_files = glob.glob(f"{self.input_source}/*.hepevt")
            if len(input_files) == 0:
                logger.error(f"Input source {self.input_source} does not contain hepevt files, exiting")
                sys.exit(1)
            logger.info(f"Found {len(input_files)} hepevt files")

            # get files and seeds
            if args.dryrun:
                self.dryrun_file = input_files[0]
                #self.dryrun_file = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/collimation/LCC_V106p2_Z/tcrh2.hepevt"
                self.dryrun_seed = os.path.basename(self.dryrun_file).replace(".hepevt", "")
            else:
                njobs = 0
                for n, f in enumerate(input_files):
                    if args.max_jobs > 0 and njobs >= args.max_jobs:
                        break
                    njobs += 1
                    seed = os.path.basename(f).replace(".hepevt", "")
                    if os.path.exists(f"{self.outdir}/{seed}.root"):
                        logger.info(f"File exists, skip")  
                    else:
                        logger.info(f"Add seed {seed}")
                        self.seeds.append(seed)
                self.transfer_input_files.append(f"{self.input_source}/$(SEED).hepevt")
                self.transfer_output_files.append(f"$(SEED).root")
            

        
        elif self.input_type == "gun":
            self.input_source = os.path.abspath(args.input_source) # assume local
            if not os.path.exists(self.input_source): 
                logger.info(f"Cannot find gun file {self.input_source}")    
                quit()

            logger.info(f"Found gun file {self.input_source}")

            multiplicity, particle, nevents = None, None, None
            thetaMin, thetaMax, momentumMin, momentumMax = None, None, None, None
            with open(self.input_source , 'r') as f:
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

            self.seeds = [random.randint(10000, 99999) for _ in range(int(args.max_jobs))]

            self.transfer_output_files.append(f"output_$(SEED).root")
            self.ddsim_xargs += f" --outputFile output_${{seed}}.root"
            self.ddsim_xargs += f" --enableGun --gun.particle {particle} --gun.multiplicity {multiplicity} --gun.thetaMin '{thetaMin}*deg' --gun.thetaMax '{thetaMax}*deg' --gun.momentumMin '{momentumMin}*GeV' --gun.momentumMax '{momentumMax}*GeV' --random.seed ${{seed}} --numberOfEvents {nevents} --gun.distribution uniform"
            self.dryrun_file = self.input_source # not necessary, but need sth to copy 
            self.dryrun_seed = "10001"

        if self.input_type == "root":
            self.ddsim_xargs += f" --inputFiles ${{seed}}.root --outputFile ddsim_${{seed}}.root --edm4hep.mcParticleCollectionName Pairs --numberOfEvents 200"

            input_files = glob.glob(f"{self.input_source}/*.root")
            if len(input_files) == 0:
                logger.error(f"Input source {self.input_source} does not contain pair files, exiting")
                sys.exit(1)
            logger.info(f"Found {len(input_files)} ROOT files")

            # get files and seeds
            if args.dryrun:
                self.dryrun_file = input_files[0]
                self.dryrun_seed = os.path.basename(self.dryrun_file).replace(".root", "")
            else:
                njobs = 0
                for n, f in enumerate(input_files):
                    if args.max_jobs > 0 and njobs >= args.max_jobs:
                        break
                    njobs += 1
                    seed = os.path.basename(f).replace(".root", "")
                    if os.path.exists(f"{self.outdir}/ddsim_{seed}.root"):
                        logger.info(f"File exists, skip")  
                    else:
                        logger.info(f"Add seed {seed}")
                        self.seeds.append(seed)
                self.transfer_input_files.append(f"{self.input_source}/$(SEED).root")
                self.transfer_output_files.append(f"ddsim_$(SEED).root")

        else:
            logging.error(f"Source {self.input_source} not implemented")
            quit()   

    def make_script(self):
        suppress_output = "" if self.args.dryrun else " &> /dev/null "
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

                
                echo "Apply geometry actions"
                {self.geometry_action}

                ls -lrt
                SECONDS=0
                echo "Running ddsim"
                ddsim --steeringFile steering.py --compactFile {self.compact_file} {self.ddsim_xargs} {suppress_output}
                if [ ! -s "${{seed}}.root" ]; then
                    echo "ERROR: expected output file ${{seed}}.root not found or empty" >&2
                    exit 100
                fi
                
                echo "ddsim step done"
                gen_duration=$SECONDS

                ls -lrt
                echo "Done script, total duration ${{gen_duration}} seconds"




                """

    def generate(self):

        # make executable script
        script_sandbox = self.make_script()
        submitFn = f"{self.outdir}/run_ddsim.sh"
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
            fOut.write(f'arguments      = $(SEED)\n')

            fOut.write(f'Log            = {logdir}/condor_job.$(ClusterId).$(ProcId).log\n')
            fOut.write(f'Output         = {logdir}/condor_job.$(ClusterId).$(ProcId).out\n')
            fOut.write(f'Error          = {logdir}/condor_job.$(ClusterId).$(ProcId).error\n')

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
            

            fOut.write(f'+JobBatchName = "DDSIM_{self.geometry_name}_{self.input_name}{self.suffix}_v{subv}"\n')
            fOut.write(f'RequestMemory  = {self.max_memory}\n')

            proxy_path = get_voms_proxy_path()
            os.system(f"cp {proxy_path} {self.outdir}/")
            fOut.write(f'use_x509userproxy     = True\n')
            fOut.write(f'x509userproxy         = {self.outdir}/{os.path.basename(proxy_path)}\n')

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
        rundir = f"/tmp/ddsim/{self.geometry_name}/{self.input_name}{self.suffix}/"

        script_init = f"""
        set -e

        rm -rf {rundir}
        mkdir -p {rundir}
        cd {rundir}
        pwd
        cp {self.steering_file} steering.py
        cp {self.dryrun_file} .
        ls -lrt

        """
        subprocess.run(["/bin/bash", "-c", script_init])

        script_sandbox = self.make_script()
        with open(f"{rundir}/run.sh", "w") as tf:
            tf.write(script_sandbox)
        subprocess.run(["/bin/bash", "run.sh", self.dryrun_seed], cwd=rundir, env={})



def main():
    producer = DDSimProducer(args)
    if args.submit:
        producer.generate()
    #if args.clean:
    #    producer.clean()
    if args.dryrun:
        producer.dryrun()


if __name__ == "__main__":
    main()
