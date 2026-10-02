
import argparse
import pybdsim
import LCC_v3short as ltlt
import sample_config

if __name__ == "__main__":
    # Parse command-line arguments for xweight and yweight
    parser = argparse.ArgumentParser(description="Run MDIStudy with varying xweight and yweight.")
    parser.add_argument("--ngenerate", type=int, default=sample_config.NGENERATE["halo"], help="Number of particles")
    parser.add_argument("--userfile", type=int, default=3, help="1--> beam core simulation; 2--> beam halo simulation; 3--> beam halo with exponential density")
    parser.add_argument("--roundmask", type=float, default=3, help="1-->round mask 7mm; 2-->elliptical mask 7x8.5; else-->jaw mask")
    parser.add_argument("--xtail", type=int, default=10, help="Horizontal halo width in sigma")
    parser.add_argument("--ytail", type=int, default=513, help="Vertical halo width in sigma")
    parser.add_argument("--xweight", type=float, default=1.45, help="X-weight value")
    parser.add_argument("--yweight", type=float, default=0.043, help="Y-weight value")
    parser.add_argument("--withDip", type=int, required=False, default=1, help="Add post IP dipoles")
    parser.add_argument("--seed", type=int, required=False, default=123456, help="Seed number")
    parser.add_argument("--withSol", type=int, default=0, help="1-->Including solenoid 0-->Without solenoid")
    parser.add_argument("--withCorr", type=int, default=0, help="1-->Including correction 0-->No correction")
    args = parser.parse_args()

    # Create an MDIStudy instance including the specified beam halo xweight and yweight
    mdi_study = ltlt.MDIStudy(seed=args.seed, withSol=args.withSol, withCorr=args.withCorr, ngenerate=args.ngenerate, withBC1L=3, withDip=args.withDip, maskA=[15e-3, 15e-3, 7e-3, 8.5e-3], userfile=args.userfile, roundmask=args.roundmask, xtail=args.xtail, ytail=args.ytail, xweight=args.xweight, yweight=args.yweight)

    print("Run genGMAD")
    mdi_study.genGMAD() # Generate the BDSIM model

    print("Run runBDSIM")
    mdi_study.runBDSIM() # Run BDSim
