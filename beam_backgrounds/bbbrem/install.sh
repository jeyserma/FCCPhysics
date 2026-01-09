
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29


git clone https://gitlab.cern.ch/hbu/bbbrem.git
cd bbbrem
cmake .
make
cd ../

