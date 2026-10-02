"""
The few constants that define the LCC_v1p1 samples and their normalisation.

Imported from both stacks -- by run_halo.py / run_core.py under the BDSIM stack
and by convert_edm4hep.py under key4hep -- so it must stay dependency-free:
standard library only. That is also why it is not part of LCC_v1p1.py, which
cannot be imported at all under key4hep.

Everything else stays in the file it belongs to: the halo tail shape in
run_halo.py, the masks / solenoid / lattice options in LCC_v1p1.py.

Note that BUNCH_INTENSITY and CHARGE_FRACTION are not simulation inputs -- BDSIM
never sees them. They are applied afterwards to turn simulated counts into rates
per bunch crossing.
"""

# Positrons per bunch, FCC-ee Z.
BUNCH_INTENSITY = 2.02e11

# Fraction of the bunch charge each sample represents.
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}

# Primary positrons per job.
#
# The halo value is much smaller than LCC_v2's 500000 because this config keeps
# the ORIGINAL v1 sampler, whose init_size = 250*ngenerate makes memory scale
# steeply: measured peak RSS of genGMAD is
#     479 MB + 16.6 kB per primary   ->  2.05 GB at 100k, 3.63 GB at 200k
# so a 4000 MB request caps this at ~215k. 200000 is also exactly what the
# original v1 campaign used, which makes this a like-for-like reference.
#
# This does NOT affect the comparison with LCC_v2: the weight is
# chargeFraction*bunchIntensity / sum(nPrimaries), so a campaign of more, smaller
# jobs normalises identically -- it just needs more jobs for the same statistics.
#
# The core is unaffected (gausstwiss is generated inside BDSIM and never calls
# generate_4d_distribution), so it keeps the v2 value.
NGENERATE = {"halo": 200000, "core": 1000000}

# Number of EDM4hep events one output file is split into.
EVENTS_PER_FILE = 1
