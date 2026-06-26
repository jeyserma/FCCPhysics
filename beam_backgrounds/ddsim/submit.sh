#!/bin/bash

PARAMETER_SETS=(
Z128_0T_grids8_vtx000
Z128_2T_grids8_vtx000
)


for PARAM in "${PARAMETER_SETS[@]}"; do
    echo "Submitting ${PARAM}"

    python submit.py \
        --input_type root \
        --input_source "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_studies/FCCee_Z_GHC_V25p1/${PARAM}/" \
        --input_name "${PARAM}" \
        --storagedir "/ceph/submit/data/group/fcc/ee/beam_backgrounds/ddsim/guineapig/ipc_studies/FCCee_Z_GHC_V25p1/" \
        --geometry CLD_o2_v07_stack \
        --submit \
        --cms_pool
        
    echo ""
done

