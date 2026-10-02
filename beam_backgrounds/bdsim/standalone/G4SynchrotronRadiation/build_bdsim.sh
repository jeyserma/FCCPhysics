#!/bin/bash
# Build BDSIM (v1.7.8) against the patched Geant4 in install/geant4 and the
# key4hep 2026-04-08 ROOT / CLHEP / Xerces-C.
#
#   bash build_bdsim.sh [njobs]
#
# Run build_geant4.sh first.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NJ="${1:-32}"
SRC="$HERE/src/bdsim"
BLD="$HERE/build/bdsim"
INS="$HERE/install/bdsim"
G4INS="$HERE/install/geant4"

unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH
set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08 >/dev/null 2>&1
set -u

G4DIR=$(dirname "$(find "$G4INS" -name Geant4Config.cmake | head -1)")
if [ -z "$G4DIR" ]; then
    echo "no Geant4Config.cmake under $G4INS -- build Geant4 first" >&2
    exit 1
fi

mkdir -p "$BLD" "$INS"
cd "$BLD"

# Geant4_DIR pins find_package to the patched build. Without it CMake would find
# key4hep's unpatched 11.4.0 through CMAKE_PREFIX_PATH and build happily against
# the wrong Geant4.
cmake "$SRC" \
  -DCMAKE_INSTALL_PREFIX="$INS" \
  -DCMAKE_BUILD_TYPE=Release \
  -DGeant4_DIR="$G4DIR" \
  -DCMAKE_PREFIX_PATH="$G4INS;${CMAKE_PREFIX_PATH}" \
  -DUSE_GDML=ON \
  -DUSE_BOOST=OFF \
  -DUSE_HEPMC3=OFF \
  -DUSE_PYTHON_BINDINGS=OFF \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_INSTALL_RPATH="$G4INS/lib64;$INS/lib"

make -j "$NJ"
make install
echo "BDSIM installed in $INS"
