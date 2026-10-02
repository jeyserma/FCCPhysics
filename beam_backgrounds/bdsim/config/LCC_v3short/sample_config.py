"""
The few constants that define the LCC_v2 samples and their normalisation.

Imported from both stacks -- by run_halo.py / run_core.py under the BDSIM stack
and by convert_edm4hep.py under key4hep -- so it must stay dependency-free:
standard library only. That is also why it is not part of LCC_v2.py, which
cannot be imported at all under key4hep.

Everything else stays in the file it belongs to: the halo tail shape in
run_halo.py, the masks / solenoid / lattice options in LCC_v2.py.

Note that BUNCH_INTENSITY and CHARGE_FRACTION are not simulation inputs -- BDSIM
never sees them. They are applied afterwards to turn simulated counts into rates
per bunch crossing.
"""

# Positrons per bunch, FCC-ee Z.
BUNCH_INTENSITY = 2.02e11

# Fraction of the bunch charge each sample represents.
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}

# Primary positrons per job. The core is run in bigger jobs because it yields
# ~23x fewer photons per primary, so a 1M-primary core file is still only ~0.5 MB
# and this halves the job count. Both are ~2 ms/primary, i.e. ~8 ms on the grid.
NGENERATE = {"halo": 1000000, "core": 1000000}

# Number of EDM4hep events one output file is split into.
EVENTS_PER_FILE = 1
