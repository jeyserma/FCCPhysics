# Standalone BDSIM with the patched Geant4 SR process

BDSIM **v1.7.8** built on **Geant4 11.4.2 + the SR patch**
(`gitlab.cern.ch/gbroggi/geant4`, branch `G4SynchrotronRadiation_patch`,
commit `86e142e4c`), linked against the key4hep `2026-04-08` stack.

```bash
source standalone/G4SynchrotronRadiation/setup.sh
bdsim --file=GMAD/input.gmad --batch --ngenerate=1000 --outfile=out --seed=1
```

## Layout

| path | what |
|---|---|
| `src/geant4`, `src/bdsim` | source checkouts |
| `build/` | build trees (not needed at run time) |
| `install/geant4` | 212 MB, patched Geant4 |
| `install/bdsim` | 19 MB |
| `build_geant4.sh`, `build_bdsim.sh` | reproducible build, in that order |
| `setup.sh` | run-time environment, relocatable |

To ship to jobs: `bash make_tarball.sh` writes `bdsim_g4sr.tgz`, **31 MB**
compressed (the full install is 230 MB, but most of that is Geant4 headers and
examples). In the job: `tar -xzf bdsim_g4sr.tgz && source setup.sh`. Tested
unpacked in a fresh directory: all 35 libraries resolve there, the LCC_v2short
model runs, no errors. Everything else — compiler runtime, ROOT, CLHEP,
Xerces-C, **and the Geant4 datasets** — comes from CVMFS, so jobs need
`/cvmfs/sw.hsf.org`.

The `build/` trees contain absolute paths from before the move into this
directory; to recompile, delete `build/` and rerun the two build scripts.

## The patch

18 commits on top of the 11.4.2 import, touching only
`source/processes/electromagnetic/xrays/{include,src}/G4SynchrotronRadiation.*`.
There is also a `G4SynchrotronRadiation_fastPatch` branch in the same fork,
not used here.

## Build choices

- **key4hep toolchain** (gcc 14.2, cmake 3.31) and libraries: ROOT 6.38.04,
  CLHEP 2.4.7.2, Xerces-C 3.3.0.
- **C++20**, matching key4hep's ROOT; BDSIM inherits it from Geant4's flags.
- **External CLHEP** (`GEANT4_USE_SYSTEM_CLHEP=ON`). BDSIM does
  `find_package(CLHEP REQUIRED)`; with Geant4's internal CLHEP two copies would
  end up in one process.
- **Sequential** Geant4, GDML on, no Qt/OpenGL.
- **Datasets not downloaded.** 11.4.2 needs exactly the versions key4hep has
  (G4EMLOW 8.8, G4NDL 4.7.1, PhotonEvaporation 6.1.2, ...), so
  `GEANT4_DATA_DIR` points at the CVMFS copy.

## Verified

1. **The patched Geant4 is what gets loaded, not key4hep's 11.4.0.** key4hep puts
   its unpatched Geant4 on `LD_LIBRARY_PATH`; `setup.sh` prepends ours.
   `ldd bdsim`: every `libG4*` from `install/geant4`, zero from key4hep.
2. **The loaded library contains the patch.** `libG4processes_core.so` exports
   `G4SynchrotronRadiation::GetRandomEnergySR(double,double,double,double&)` and
   `fLowestPhotonEnergy`; key4hep's has only the original 3-argument version.
3. **It runs the FCC model.** `LCC_v2short` halo, 5000 primaries, seed 777, same
   input for both:

   | | photons at sampler / primary | median E | E > 100 keV / primary |
   |---|---|---|---|
   | CVMFS BDSIM 1.7.7 / Geant4 10.7.2.3 | 2.018 | 7.99 keV | 0.169 |
   | this build, 1.7.8 / Geant4 11.4.2 + patch | 2.180 | 6.97 keV | 0.177 |

   31 s, 528 MB peak RSS, no warnings.

   The ~8 % difference mixes three changes (Geant4 10.7 → 11.4, the SR patch,
   BDSIM 1.7.7 → 1.7.8). Isolating the patch needs an unpatched 11.4.2 build
   with the same flags.

## pybdsim

key4hep already ships `pybdsim` 3.6.1 and `pymadx`, so after `source setup.sh`
the whole job — `genGMAD()` (pybdsim), `bdsim`, and `convert_edm4hep.py` (podio)
— can run in **one** environment. `pybdsim.Run.Bdsim` calls `bdsim` from `PATH`,
which is this build. This removes the old BDSIM-stack / key4hep split and the
`PYTHONHOME` workaround in `convert.sh`, but that has not been wired into
`submit.py` yet.
