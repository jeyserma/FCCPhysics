#!/bin/bash

ACC="FCCee_Z_GHC_V25p1"
PAR="CFG_DEF_128"

CURR_DIR=$(pwd)

TAG="${ACC}_${PAR}"
DIR="/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/visualization/"

GP_DIR="/work/submit/jaeyserm/fccee/FCCAnalyses/FCCPhysics/beam_backgrounds/guineapig/guinea-pig-15122025_dev_time/"



cd $GP_DIR
source env.sh

RUNDIR="$DIR/$TAG"

mkdir -p $RUNDIR
cd $RUNDIR

$GP_DIR/build/src/guinea --acc_file $CURR_DIR/acc.dat $ACC $PAR output

cd $CURR_DIR
