#!/bin/bash
# Conversion step of an LCC_v2 job: BDSIM ROOT -> EDM4hep.
#
#     bash convert.sh <seed> <runconfig> [collection]
#
# where <runconfig> is a key of CHARGE_FRACTION in sample_config.py (halo, core, ...).
#
# This runs in a SEPARATE stack from the generation. run_halo.py needs the BDSIM
# stack; convert_edm4hep.py needs key4hep for podio/edm4hep. Both live on CVMFS,
# so a single job can source one and then the other -- but not in the same shell,
# hence this wrapper. Nothing here depends on pybdsim: the BDSIM sampler is read
# through the StreamerInfo stored in the file itself.
#
# The event split and the normalisation constants come from sample_config.py.
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 SEED RUNCONFIG [COLLECTION]" >&2
    echo "  RUNCONFIG is a key of CHARGE_FRACTION in sample_config.py (halo, core, ...);" >&2
    echo "  it selects the charge fraction and is recorded in the output metadata." >&2
    exit 1
fi
seed="$1"
runconfig="$2"
collection="${3:-MCParticle}"

KEY4HEP_SETUP="/cvmfs/sw.hsf.org/key4hep/setup.sh"
KEY4HEP_RELEASE="2026-04-08"

# The BDSIM stack has already been sourced by the parent job script, and it sets
# PYTHONHOME to its own LCG Python 3.9. That leaks into key4hep's Python 3.13,
# which then cannot find its own stdlib and dies with
#     Fatal Python error: Failed to import encodings module
# Clear the Python and library search paths before layering key4hep on top. This
# is a dedicated subshell for the conversion, so nothing else is affected.
unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH

echo "Sourcing key4hep ${KEY4HEP_RELEASE}"
set +u
if ! source "${KEY4HEP_SETUP}" -r "${KEY4HEP_RELEASE}" >/dev/null 2>&1; then
    echo "ERROR: failed to source ${KEY4HEP_SETUP} -r ${KEY4HEP_RELEASE}" >&2
    exit 1
fi
set -u

python convert_edm4hep.py --input "output_${seed}.root" \
                          --output "output_${seed}_edm4hep.root" \
                          --runconfig "${runconfig}" \
                          --collection "${collection}"

if [ ! -f "output_${seed}_edm4hep.root" ]; then
    echo "ERROR: conversion did not produce output_${seed}_edm4hep.root" >&2
    exit 1
fi
