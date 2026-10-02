#!/bin/bash
# Build the patched Geant4 (gbroggi/geant4, branch G4SynchrotronRadiation_patch,
# on top of Geant4 11.4.2) against the key4hep 2026-04-08 toolchain.
#
#   bash build_geant4.sh [njobs]
#
# Choices, and why:
#   - key4hep gcc 14.2 / cmake / Xerces-C, so the result links cleanly with
#     key4hep's ROOT 6.38 when BDSIM is built on top.
#   - C++20: the standard key4hep's ROOT and Geant4 are built with. BDSIM takes
#     its standard from Geant4's flags, so this propagates.
#   - external CLHEP from key4hep (2.4.7.2): BDSIM does find_package(CLHEP REQUIRED),
#     so Geant4 must use the same one, or two CLHEPs end up in the same process.
#   - sequential (no MT): BDSIM tracks in serial.
#   - no Qt / OpenGL: batch jobs only; keeps the install small to ship.
#   - datasets NOT downloaded: 11.4.2 needs exactly the dataset versions key4hep
#     already has on CVMFS, so GEANT4_INSTALL_DATADIR points there and geant4.sh
#     sets the G4*DATA variables to CVMFS paths.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NJ="${1:-32}"
SRC="$HERE/src/geant4"
BLD="$HERE/build/geant4"
INS="$HERE/install/geant4"
G4DATA="/cvmfs/sw.hsf.org/key4hep/releases/2026-02-01/x86_64-almalinux9-gcc14.2.0-opt/geant4-data/11.4.0-vfilzd/share/geant4-data-11.4.0"

unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH
set +u
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08 >/dev/null 2>&1
set -u

mkdir -p "$BLD" "$INS"
cd "$BLD"

cmake "$SRC" \
  -DCMAKE_INSTALL_PREFIX="$INS" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=20 \
  -DGEANT4_BUILD_MULTITHREADED=OFF \
  -DGEANT4_USE_GDML=ON \
  -DGEANT4_USE_QT=OFF \
  -DGEANT4_USE_OPENGL_X11=OFF \
  -DGEANT4_USE_SYSTEM_CLHEP=ON \
  -DGEANT4_INSTALL_DATA=OFF \
  -DGEANT4_INSTALL_DATADIR="$G4DATA" \
  -DGEANT4_BUILD_TLS_MODEL=global-dynamic

make -j "$NJ"
make install
echo "Geant4 installed in $INS"
