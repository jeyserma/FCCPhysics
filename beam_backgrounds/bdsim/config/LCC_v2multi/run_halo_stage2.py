"""
Stage 2 of the two-stage halo run: start at the split element, read the
replicated stage-1 dump as the beam, and track to the QD0AL sampler.
Run make_stage2_input.py first -- it writes GMAD/stage2_input.dat.

    source env.sh
    python run_halo_stage2.py --seed 12345
"""

import argparse
import os
import sys

import pybdsim
import LCC_v2multi as ltlt
import sample_config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LCC_v2multi halo, stage 2.")
    parser.add_argument("--ngenerate", type=int,
                        default=sample_config.NGENERATE["halo"],
                        help="Unused in stage 2; the count comes from the input file")
    parser.add_argument("--userfile", type=int, default=3,
                        help="1: gaussian core, 2: uniform halo, 3: exponential halo")
    parser.add_argument("--roundmask", type=float, default=3)
    parser.add_argument("--xtail", type=int, default=10)
    parser.add_argument("--ytail", type=int, default=513)
    parser.add_argument("--xweight", type=float, default=1.2)
    parser.add_argument("--yweight", type=float, default=0.06)
    parser.add_argument("--withDip", type=int, default=1)
    parser.add_argument("--seed", type=int, default=123456)
    parser.add_argument("--withSol", type=int, default=0)
    parser.add_argument("--withCorr", type=int, default=0)
    parser.add_argument("--slim", action="store_true",
                        help="Drop the tunnel from the stage-2 model. Much faster "
                             "to navigate, and every injected particle pays the "
                             "navigation cost on entry. Validate before using.")
    parser.add_argument("--species", choices=("both", "e+", "gamma"), default="both",
                        help="Which species to inject; 'both' (default) runs the two "
                             "passes back to back and writes output_<seed>_ep.root and "
                             "output_<seed>_gamma.root. They are always two BDSIM runs: "
                             "BDSIM applies the REFERENCE particle's mass to every line "
                             "regardless of pdgid, so one mixed run would lose every "
                             "photon below 511 keV. See README.md section 4.")
    args = parser.parse_args()

    species_list = ["e+", "gamma"] if args.species == "both" else [args.species]
    suffix = {"e+": "_ep", "gamma": "_gamma"}

    outputs = []
    for species in species_list:
        infile = os.path.join("GMAD", ltlt.STAGE2_INPUT[species])
        if not os.path.exists(infile):
            sys.exit(f"{infile} not found -- run make_stage2_input.py on the "
                     "stage-1 output first")
        # BDSIM must generate exactly as many primaries as there are lines, or it
        # will loop back to the start of the file and silently double-count.
        with open(infile) as f:
            nlines = sum(1 for line in f if line.strip() and not line.startswith("#"))
        print(f"stage 2 [{species}]: {nlines} particles in {infile}")

        mdi_study = ltlt.MDIStudy(
            seed=args.seed, withSol=args.withSol, withCorr=args.withCorr,
            ngenerate=nlines, withBC1L=3, withDip=args.withDip,
            maskA=[15e-3, 15e-3, 7e-3, 8.5e-3], userfile=args.userfile,
            roundmask=args.roundmask, xtail=args.xtail, ytail=args.ytail,
            xweight=args.xweight, yweight=args.yweight,
            stage=2, stage2_species=species, slimGeometry=args.slim,
            runKey=suffix[species] if len(species_list) > 1 else "",
        )
        print("Run genGMAD")
        mdi_study.genGMAD()
        print("Run runBDSIM")
        mdi_study.runBDSIM()
        outputs.append(f"output_{args.seed}"
                       f"{suffix[species] if len(species_list) > 1 else ''}.root")

    if len(outputs) > 1:
        print()
        print("Convert BOTH passes together -- they are two halves of one sample:")
        print("  python convert_edm4hep.py "
              + " ".join(f"--input {o}" for o in outputs) + " --runconfig halo")
