#!/bin/bash
# Pack the run-time part of the install for shipping to condor jobs.
#
#   bash make_tarball.sh [output.tgz]      (default: bdsim_g4sr.tgz next to this script)
#
# In the job:   tar -xzf bdsim_g4sr.tgz && source setup.sh
#
# ~31 MB compressed, against 230 MB for the full install. Dropped:
#   install/geant4/include                 31 MB  compile-time only
#   install/geant4/share/Geant4/examples  ~114 MB  never used
#   */cmake, */pkgconfig                           compile-time only
# KEPT on purpose: install/bdsim/include (3.3 MB). Without it ROOT's cling fails
# to autoload GMAD::BeamBase and prints "Missing FileEntry for
# ../parser/beamBase.h" at startup. The output was still correct in testing, but
# an error on every job is exactly the noise that hides a real failure.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$HERE/bdsim_g4sr.tgz}"

tar -czf "$OUT" -C "$HERE" \
    --exclude='install/geant4/include' \
    --exclude='install/geant4/share/Geant4/examples' \
    --exclude='*/cmake' \
    --exclude='*/pkgconfig' \
    install setup.sh

echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
