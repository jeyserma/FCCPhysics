"""
Stage 1 of the two-stage halo run: track the exponential halo from the start
of the tracked region up to the split element and sample EVERY species
crossing that plane. Writes output_<seed>.root, which make_stage2_input.py
turns into the stage-2 beam file.

    source env.sh
    python run_halo_stage1.py --seed 12345
"""

import argparse
import os
import sys

import pybdsim
import LCC_v2multi as ltlt
import sample_config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LCC_v2multi halo, stage 1.")
    parser.add_argument("--ngenerate", type=int,
                        default=sample_config.NGENERATE["halo"],
                        help="Number of primary positrons (stage 1 only; stage 2 "
                             "takes its count from the input file)")
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
    args = parser.parse_args()

    ngenerate = args.ngenerate

    mdi_study = ltlt.MDIStudy(
        seed=args.seed, withSol=args.withSol, withCorr=args.withCorr,
        ngenerate=ngenerate, withBC1L=3, withDip=args.withDip,
        maskA=[15e-3, 15e-3, 7e-3, 8.5e-3], userfile=args.userfile,
        roundmask=args.roundmask, xtail=args.xtail, ytail=args.ytail,
        xweight=args.xweight, yweight=args.yweight,
        stage=1,
    )

    print("Run genGMAD")
    mdi_study.genGMAD()

    print("Run runBDSIM")
    mdi_study.runBDSIM()
