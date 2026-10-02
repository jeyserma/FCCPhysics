
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
    part['rp'] = np.sqrt((part['x']+part['xp']*(2.4+6))**2 + (part['y']+part['yp']*(2.4+6))**2)
    # Filtering
    # in energy > 2keV
    # in radius removing photons beyond R=18mm, should not have an affect with the smaller sampling surface
    # removing all photons R<9mm at 8.4m downstream the sampler equivalent to s=+6.0m
    part = part[(part['partID']==22)&(part['E']>2e-6)&(part['r']<18e-3)&(part['rp']>9e-3)]
    part.reset_index(drop=True, inplace=True)
    result = part
    # Creating file in hepevt format following https://hugonweb.com/hepevt/
    # Loop through the photons and add each to a new vertex
    f = open(outfile, "w")
    lines = [
            # Write total number of events in the file
            f'{len(result)}\n'
        ]
    f.writelines(lines)
    for i in range(len(result)):
        # Conversion from sampler units to hepevt units
        # 'time' is in ns AT THE SAMPLER, converted to mm/c w.r.t. the IP
        # 'p' is in GeV/c
        # 'E' is in GeV
        # 'xp', 'yp', 'zp' are in fractional part of the momentum p, convert to GeV/c
        # 'x', 'y' are in m, convert to mm
        # z is -2.4m, location of the sampler
        partID = int(result.iloc[i]['partID'])
        px, py, pz, E = result.iloc[i]['xp']*result.iloc[i]['p'], result.iloc[i]['yp']*result.iloc[i]['p'], result.iloc[i]['zp']*result.iloc[i]['p'], result.iloc[i]['E']
        x, y, z, t = 1e3*result.iloc[i]['x'], 1e3*result.iloc[i]['y'], -2.4e3, 1e3*(result.iloc[i]['t']*1e-9*2.997924588e8-(d.model['QD0AL']['SEnd']+2.4))
        lines = [
            # Create an active photon particle (PDG ID for photon is 22)
            # <Status> <PDG ID> <1st Mother> <2nd Mother> <1st Daughter> <2nd Daughter> <Px> <Py> <Px> <E> <Mass> <x> <y> <z> <t>
            # where Px/Py/Pz are in GeV/c, E is in GeV, and M is in GeV/c^2. x/y/z are in mm and t is in mm/c
            f'1 {partID} 0 0 0 0 {px} {py} {pz} {E} {0} {x} {y} {z} {t}\n'
        ]
        f.writelines(lines)
    
    f.close()
