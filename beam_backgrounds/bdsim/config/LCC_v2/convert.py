
import numpy as np
import pandas as pd
import sys
import pybdsim
import pymadx
import ROOT
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from pathlib import Path
import re
import argparse

if __name__ == "__main__":

    # Parse command-line arguments for xweight and yweight
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True, help="ROOT input file")
    args = parser.parse_args()


    infile = args.input
    outfile = infile.replace(".root", ".hepevt")

    # Read the dist file
    d = pybdsim.Data.Load(infile)
    part2 = pybdsim.Data.SamplerData(d, -1) # Only one sampler, hence take the last one, 0 is the initial positron disitribution
    part = pd.DataFrame([part2.data['x'], part2.data['xp'], part2.data['y'], part2.data['yp'], part2.data['energy'], part2.data['p'], part2.data['zp'],
    part2.data['T'], part2.data['partID'], part2.data['mass']]).T
    part.columns = ['x', 'xp', 'y', 'yp', 'E', 'p', 'zp', 't', 'partID', 'm']
    part['r'] = np.sqrt(part['x']**2 + part['y']**2)
    # Straight-line extrapolation 8.4 m downstream of the sampler, i.e. to s=+6.0 m.
    # 'xp'/'yp'/'zp' are momentum direction cosines, so the slopes are xp/zp and yp/zp.
    # LCC_v1 used xp and yp directly; identical for zp~1, but wrong-signed for the
    # ~0.1% of photons travelling backwards, which are now removed anyway.
    part['rp'] = np.sqrt((part['x']+part['xp']/part['zp']*(2.4+6))**2 + (part['y']+part['yp']/part['zp']*(2.4+6))**2)
    # Filtering
    # in energy > 2keV
    # in radius removing photons beyond R=18mm, should not have an affect with the smaller sampling surface
    # removing all photons R<9mm at 8.4m downstream the sampler equivalent to s=+6.0m
    # keeping only forward-going photons (zp>0); the ~0.1% travelling back up the
    # beamline cannot reach the detector and the extrapolation above is not valid for them
    part = part[(part['partID']==22)&(part['zp']>0)&(part['E']>2e-6)&(part['r']<18e-3)&(part['rp']>9e-3)]
    part.reset_index(drop=True, inplace=True)
    result = part

    # Conversion from sampler units to hepevt units
    # 'time' is in ns AT THE SAMPLER, converted to mm/c w.r.t. the IP
    # 'p' is in GeV/c
    # 'E' is in GeV
    # 'xp', 'yp', 'zp' are in fractional part of the momentum p, convert to GeV/c
    # 'x', 'y' are in m, convert to mm
    # z is -2.4m, location of the sampler
    # Vectorised with respect to LCC_v1, which called result.iloc[i] ten times per
    # photon; the output format is unchanged.
    partID = result['partID'].to_numpy().astype(int)
    p = result['p'].to_numpy()
    px, py, pz = result['xp'].to_numpy()*p, result['yp'].to_numpy()*p, result['zp'].to_numpy()*p
    E = result['E'].to_numpy()
    x, y = 1e3*result['x'].to_numpy(), 1e3*result['y'].to_numpy()
    z = -2.4e3
    t = 1e3*(result['t'].to_numpy()*1e-9*2.997924588e8 - (d.model['QD0AL']['SEnd']+2.4))

    # Creating file in hepevt format following https://hugonweb.com/hepevt/
    # Each photon is a separate vertex of the single event in the file
    # <Status> <PDG ID> <1st Mother> <2nd Mother> <1st Daughter> <2nd Daughter> <Px> <Py> <Pz> <E> <Mass> <x> <y> <z> <t>
    # where Px/Py/Pz are in GeV/c, E is in GeV, and M is in GeV/c^2. x/y/z are in mm and t is in mm/c
    with open(outfile, "w") as f:
        # Write total number of particles in the file
        f.write(f'{len(result)}\n')
        f.writelines(
            f'1 {partID[i]} 0 0 0 0 {px[i]} {py[i]} {pz[i]} {E[i]} {0} {x[i]} {y[i]} {z} {t[i]}\n'
            for i in range(len(result))
        )
