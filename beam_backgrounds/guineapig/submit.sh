#!/bin/bash

PARAMETER_SETS=(
Z128_2Tm_grids8_nbf
Z128_2Tm_grids7_nbf
Z128_2Tm_grids6_nbf
Z128_2Tm_grids5_nbf
Z128_2Tm_grids4_nbf
Z128_2Tm_grids3_nbf
Z128_2Tm_grids2_nbf
Z128_2Tm_grids1_nbf
)

for PARAM in "${PARAMETER_SETS[@]}"; do
    echo "Submitting ${PARAM}"

    python submit.py \
        --accelerator FCCee_Z_GHC_V25p1 \
        --parameter_set "${PARAM}" \
        --input_file cards/studies.dat \
        --storagedir /ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_studies \
        --submit \
        --njobs 5100 \
        --cms_pool
        
    echo ""
done

