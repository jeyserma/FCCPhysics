#!/bin/bash

PARAMETER_SETS=(
Z64_2T_grids8_nbf_noBoundary
Z256_2T_grids8_nbf_noBoundary
)

for PARAM in "${PARAMETER_SETS[@]}"; do
    echo "Merging ${PARAM}"

    python submit.py \
        --accelerator FCCee_Z_GHC_V25p1 \
        --parameter_set "${PARAM}" \
        --input_file cards/studies.dat \
        --storagedir /ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_studies \
        --merge 

    echo ""
done

