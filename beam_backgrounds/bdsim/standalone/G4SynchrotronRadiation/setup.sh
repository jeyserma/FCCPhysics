#!/bin/bash
# Environment for the standalone BDSIM + patched Geant4.
#
#   source setup.sh
#
# Relocatable: paths are taken from this file's location, so the same script
# works on the submit node and inside an unpacked job sandbox.
#
# ORDER MATTERS. key4hep ships an UNPATCHED Geant4 11.4.0 and puts its libraries
# on LD_LIBRARY_PATH. The patched libraries must come first, or the dynamic
# linker can load key4hep's libG4processes and the SR patch is silently not used.
# check_patch.sh verifies which one is actually loaded.

_SA_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# key4hep release the install was BUILT against (gcc, ROOT, CLHEP, Xerces ABI).
# Pinned on purpose: never source the untagged setup.sh, which follows the latest
# release and would silently break the shipped libraries.
KEY4HEP_SETUP="/cvmfs/sw.hsf.org/key4hep/setup.sh"
KEY4HEP_RELEASE="2026-04-08"

if ! timeout 15 ls "$KEY4HEP_SETUP" >/dev/null 2>&1; then
    echo "ERROR: key4hep not accessible: $KEY4HEP_SETUP" >&2
    unset _SA_HERE
    return 1
fi

unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH
set +u
if ! source "$KEY4HEP_SETUP" -r "$KEY4HEP_RELEASE" >/dev/null 2>&1; then
    echo "ERROR: failed to source $KEY4HEP_SETUP -r $KEY4HEP_RELEASE" >&2
    unset _SA_HERE
    return 1
fi
echo "key4hep release: $KEY4HEP_RELEASE"

# geant4.sh sets the G4*DATA variables (pointing at the key4hep datasets on
# CVMFS, which are the versions 11.4.2 needs) and prepends its lib dir.
pushd "$_SA_HERE/install/geant4/bin" >/dev/null
source ./geant4.sh
popd >/dev/null

export PATH="$_SA_HERE/install/bdsim/bin:$PATH"
export LD_LIBRARY_PATH="$_SA_HERE/install/geant4/lib64:$_SA_HERE/install/bdsim/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$_SA_HERE/install/bdsim/include/bdsim:$_SA_HERE/install/bdsim/include/bdsim/analysis:${ROOT_INCLUDE_PATH:-}"
# (nounset is left off: a sourced script should not change the caller's shell
#  options, and key4hep's own setup does not survive `set -u`.)

unset _SA_HERE
