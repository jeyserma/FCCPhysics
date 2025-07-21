
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10

git clone https://gitlab.cern.ch/clic-software/guinea-pig.git guinea-pig_devTrackingWindow
cd guinea-pig_devTrackingWindow
mkdir gp
export INSTALL_DIR=$('pwd')/gp
./configure -prefix=$INSTALL_DIR
make
make install
