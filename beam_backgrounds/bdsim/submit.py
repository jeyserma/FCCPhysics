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
parser.add_argument("--analysis", action='store_true', help="report on the sample for the given runconfig (files, size, weight, what is still needed), then exit")

parser.add_argument("--njobs", type=int, help="number of jobs", default=10)
parser.add_argument("--ngenerate", type=int, help="override the ngenerate from config/<lattice>/sample_config.py", default=None)
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


# Where BDSIM comes from. Pick one with BDSIM_SOURCE.
#   /cvmfs/...  -> an environment script, sourced on the worker node
#   local path  -> a runtime tarball (see standalone/*/make_tarball.sh), shipped
#                  with each job and run on top of key4hep
BDSIM_SOURCES = {
    "bdsim_def":     "/cvmfs/beam-physics.cern.ch/bdsim/x86_64-el9-gcc13-opt/bdsim-env-v1.7.7-g4v10.7.2.3-ftfp-boost.sh",
    "bdsim_def_G4SR":  "standalone/G4SynchrotronRadiation/bdsim_g4sr.tgz",   # BDSIM 1.7.8 + Geant4 11.4.2 SR patch
}
BDSIM_SOURCE = "bdsim_def_G4SR"


def sample_name(lattice, runconfig, suffix=""):
    """
    <lattice>_<runconfig>_<source>[_<suffix>], e.g. LCC_v2short_halo_standalone_G4SR.
    python/sample_name.py parses it back.
    """
    return f"{lattice}_{runconfig}_{BDSIM_SOURCE}{suffix}"

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


def analyse_sample(args):
    """
    Report on a sample: what is on disk, how much of a bunch it represents, and
    what completing it would cost.

        weight = chargeFraction * bunchIntensity / (nfiles * ngenerate)

    Filesystem only -- it deliberately does not open the output files, so it
    needs no key4hep and stays fast on a directory of any size. The primary count
    therefore assumes every file used the ngenerate reported below.
    """
    configdir = os.path.join(os.getcwd(), "config", args.lattice)
    sys.path.insert(0, configdir)
    try:
        import sample_config
    except ImportError:
        logger.error(f"No sample_config.py in {configdir}")
        return

    if args.runconfig not in sample_config.CHARGE_FRACTION:
        logger.error(f"Unknown runconfig '{args.runconfig}'; sample_config.py knows "
                     f"{sorted(sample_config.CHARGE_FRACTION)}")
        return

    fraction = sample_config.CHARGE_FRACTION[args.runconfig]
    intensity = sample_config.BUNCH_INTENSITY
    ngenerate = args.ngenerate or sample_config.NGENERATE[args.runconfig]

    particles = fraction * intensity          # positrons in one bunch for this sample
    nfiles_w1 = particles / ngenerate         # files needed for weight 1

    suffix = f"_{args.suffix}" if args.suffix else ""
    outdir = f"{args.storagedir}/{sample_name(args.lattice, args.runconfig, suffix)}/"
    files = glob.glob(os.path.join(outdir, "output_*_edm4hep.root"))
    if not files:   # fall back to the v1 naming
        files = glob.glob(os.path.join(outdir, "output_*.hepevt"))
    nfiles = len(files)
    total_bytes = sum(os.path.getsize(f) for f in files) if files else 0
    mean_bytes = total_bytes / nfiles if nfiles else 0
    GB = 1024.0**3

    print(f"\n{args.lattice} / {args.runconfig}")
    print(f"  bunch intensity        {intensity:.4g}")
    print(f"  charge fraction        {fraction:g}")
    print(f"  particles in one bunch {particles:.4g}")
    print(f"  ngenerate per file     {ngenerate}")

    print(f"\n  {outdir}")
    print(f"  files                  {nfiles}")
    print(f"  total size             {total_bytes / GB:.2f} GB")
    if nfiles:
        print(f"  mean file size         {mean_bytes / 1e6:.2f} MB")
        print(f"  primaries simulated    {nfiles * ngenerate:.4g}")
        print(f"  current weight         {particles / (nfiles * ngenerate):.6g}")
    else:
        print(f"  current weight         - (no files yet)")

    print(f"\n  files for weight 1     {nfiles_w1:.0f}")
    print(f"  still needed           {max(0, nfiles_w1 - nfiles):.0f}")
    if mean_bytes:
        print(f"  size at weight 1       {nfiles_w1 * mean_bytes / GB:.1f} GB")
    print()



class BDSIMProducer:

    def __init__(self, args):
        self.args = args
        self.cwd = os.getcwd()

        self.njobs = args.njobs
        # ngenerate lives in config/<lattice>/sample_config.py; only passed on the
        # command line when explicitly overridden (handy for short test runs)
        self.ngenerate_arg = f"--ngenerate={args.ngenerate}" if args.ngenerate else ""

        self.lattice = args.lattice
        self.runconfig = args.runconfig
        self.configdir = f"{self.cwd}/config/{self.lattice}/"
        
        self.storagedir = args.storagedir
        self.max_memory = args.max_memory  
        
        self.suffix = f"_{args.suffix}" if args.suffix else ""
        self.name = sample_name(self.lattice, self.runconfig, self.suffix)
        self.logdir = f"{self.cwd}/{args.logdir}/{self.name}/"
        self.outdir = f"{self.storagedir}/{self.name}/"

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # pack sandbox
        logger.info("Creating sandbox")
        self.sandbox = f"{self.outdir}/sandbox.tar"
        os.system(f"tar -cvf {self.sandbox} -C {self.configdir} .")

        self.transfer_input_files = ['sandbox.tar']

        # BDSIM source: a CVMFS environment script, or a local tarball to ship
        self.bdsim_path = BDSIM_SOURCES[BDSIM_SOURCE]
        self.bdsim_mode = "cvmfs" if self.bdsim_path.startswith("/cvmfs/") else "standalone"
        if self.bdsim_mode == "standalone":
            src = self.bdsim_path
            if not os.path.isabs(src):
                src = os.path.join(os.path.dirname(os.path.abspath(__file__)), src)
            if not os.path.exists(src):
                logger.error(f"{src} not found -- run make_tarball.sh in {os.path.dirname(src)}")
                sys.exit(1)
            self.tarball_name = os.path.basename(src)
            self.standalone_tarball = f"{self.outdir}/{self.tarball_name}"
            shutil.copy(src, self.standalone_tarball)
            self.transfer_input_files.append(self.tarball_name)
        logger.info(f"BDSIM source '{BDSIM_SOURCE}' ({self.bdsim_mode}): {self.bdsim_path}")

        self.output_pattern = 'output_$(SEED)_edm4hep.root'   # was output_$(SEED).hepevt in v1
        self.transfer_output_files  = [self.output_pattern]

        njob = 0
        self.seeds = []
        while njob < self.njobs:
            seed = f"{random.randint(100000,999999)}"
            outputFile = os.path.join(self.outdir, self.output_pattern.replace("$(SEED)", seed))
            if os.path.exists(outputFile):
                logger.warning(f"Output file with seed {seed} already exists, skipping")
                continue
            self.seeds.append(seed)
            njob += 1



    def env_block(self):
        """The only part of the job script that differs between the two modes."""
        if self.bdsim_mode == "standalone":
            return f"""
            echo "Unpack standalone BDSIM"
            mkdir -p bdsim_standalone
            if ! tar -xzf {self.tarball_name} -C bdsim_standalone; then
                echo "ERROR: Failed to unpack {self.tarball_name}" >&2
                exit 1
            fi

            # setup.sh checks CVMFS and sources the pinned key4hep release itself
            echo "Sourcing standalone BDSIM (key4hep + patched Geant4)"
            set +u
            if ! source bdsim_standalone/setup.sh; then
                echo "ERROR: Failed to source bdsim_standalone/setup.sh" >&2
                exit 1
            fi
            set -u
            echo "bdsim: $(which bdsim)"
            # guard: the patched Geant4 must be the one loaded, not key4hep's 11.4.0
            if ldd "$(which bdsim)" | grep -q "key4hep.*geant4"; then
                echo "ERROR: bdsim resolves Geant4 from key4hep, not the shipped build" >&2
                ldd "$(which bdsim)" | grep libG4 >&2
                exit 1
            fi
"""
        return f"""
            echo "Checking CVMFS"
            set +u
            if ! timeout 15 ls "{self.bdsim_path}" >/dev/null 2>&1; then
                echo "ERROR: Stack not found or not accessible: {self.bdsim_path}" >&2
                exit 1
            fi

            echo "Sourcing stack"
            if ! source "{self.bdsim_path}" >/dev/null 2>&1; then
                echo "ERROR: Failed to source stack: {self.bdsim_path}" >&2
                exit 1
            fi
            set -u
"""

    def convert_cmd(self):
        """
        convert.sh sources key4hep itself. In standalone mode key4hep is already
        active, and key4hep refuses a second setup in the same environment ("The
        Key4hep software stack is already set up, please start a new shell"), so
        the conversion runs in a clean environment instead. convert.sh itself is
        unchanged and still works as before in cvmfs mode.
        """
        # BDSIM_SOURCE is picked up by convert_edm4hep.py and stored as bdsimSource
        if self.bdsim_mode == "standalone":
            return f'env -i HOME="$HOME" TMPDIR="${{TMPDIR:-/tmp}}" PATH=/usr/bin:/bin BDSIM_SOURCE={BDSIM_SOURCE} bash convert.sh'
        return f"BDSIM_SOURCE={BDSIM_SOURCE} bash convert.sh"

    def make_script(self, runconfig):
        return f"""
            #!/bin/bash
            set -euo pipefail

            if [ $# -lt 1 ]; then
                echo "Usage: $0 SEED" >&2
                exit 1
            fi

            seed="$1"
            echo $seed
            # ngenerate comes from config/{self.lattice}/sample_config.py

            echo "Release:"
            cat /proc/version

            echo "Hostname:"
            hostname

            echo "Current working directory:"
            pwd

            echo "List current working dir"
            ls -lrt

{self.env_block()}
            echo "Unpack sandbox"
            if ! tar -xf sandbox.tar; then
                echo "ERROR: Failed to unpack sandbox.tar" >&2
                exit 1
            fi
            ls -lrt

            SECONDS=0

            echo "Running generation"
            python "run_{runconfig}.py" --seed="$seed" {self.ngenerate_arg}

            if [ ! -f "output_${{seed}}.root" ]; then
                echo "ERROR: Generation did not produce output_${{seed}}.root" >&2
                exit 1
            fi

            echo "Generation step done"
            gen_duration=$SECONDS

            echo "Running conversion"
            #python convert.py --input "output_${{seed}}.root"
            {self.convert_cmd()} "${{seed}}" {runconfig}

            if [ ! -f "output_${{seed}}_edm4hep.root" ]; then
                echo "ERROR: Conversion did not produce output_${{seed}}_edm4hep.root" >&2
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
            fOut.write(f'arguments      = $(SEED)\n')

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
            

        
            fOut.write(f'+JobBatchName = "BDSIM_{self.name}_v{subv}"\n')
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
        rundir = f"{self.outdir}/tmp/"      # inside the sample directory, wiped on every dryrun

        script_init = f"""
        set -e

        rm -rf {rundir}
        mkdir -p {rundir}
        cd {rundir}
        pwd

        cp {self.sandbox} .
        {"cp " + self.standalone_tarball + " ." if self.bdsim_mode == "standalone" else ""}
        ls -lrt

        """
        subprocess.run(["/bin/bash", "-c", script_init])

        script_sandbox = self.make_script(self.runconfig)
        with open(f"{rundir}/run.sh", "w") as tf:
            tf.write(script_sandbox)
        # clean environment, like a condor job: an already-sourced key4hep in the
        # calling shell makes key4hep's setup.sh refuse to run again
        clean_env = {k: os.environ[k] for k in ("HOME", "USER", "LOGNAME", "TMPDIR", "X509_USER_PROXY") if k in os.environ}
        clean_env["PATH"] = "/usr/bin:/bin"
        subprocess.run(["bash", "run.sh", "12345"], cwd=rundir, env=clean_env)


def main():
    if args.analysis:   # query only: no sandbox, no output directory
        analyse_sample(args)
        return

    producer = BDSIMProducer(args)
    if args.submit:
        producer.generate_submit()
    if args.dryrun:
        producer.dryrun()


if __name__ == "__main__":
    main()
