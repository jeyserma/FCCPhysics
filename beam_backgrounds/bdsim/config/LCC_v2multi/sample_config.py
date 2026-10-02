"""
The few constants that define the LCC_v2multi samples and their normalisation.

Same role as LCC_v2/sample_config.py, plus the two-stage ("dual sampling")
parameters. Imported from both stacks, so it must stay dependency-free:
standard library only.

Normalisation with two stages
-----------------------------
Stage 1 generates NGENERATE["halo"] primary positrons and dumps everything
crossing the split plane. Each stage-2 run replays that whole dump once, so its
output represents NGENERATE primaries -- exactly like a single-stage file.

Run stage 2 K times with different seeds and the K outputs sum to K * NGENERATE
by themselves. That is the whole trick: nothing carries a per-particle weight,
nothing has to be told what K is, and convert_edm4hep.py, the metadata,
sr_weight() and ddsim all keep working unchanged.
"""

# Positrons per bunch, FCC-ee Z.
BUNCH_INTENSITY = 2.02e11

# Fraction of the bunch charge each sample represents.
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}

# Primary positrons per stage-1 job.
NGENERATE = {"halo": 500000, "core": 1000000}

# Number of EDM4hep events one output file is split into.
EVENTS_PER_FILE = 1

# ---------------------------------------------------------------- two-stage --

# NOTE: there is deliberately no N_REPLAY knob.
#
# The replay factor K is simply how many times you run stage 2 on the same
# stage-1 dump with different --seed. Replicating inside the bridge instead
# would be computationally identical but multiply the stage-2 input files by K
# (251 GB -> 1.25 TB for a w=1 stage 1) and prevent the passes from going to
# separate batch jobs. Stage-2 setup is ~9 s, so K separate runs cost nothing
# extra.
#
# Each stage-2 output therefore carries nPrimariesTotal = NGENERATE, and the K
# files sum to NGENERATE * K on their own -- sr_weight() and read_sample()
# already do that summation, so nothing needs to be told what K is.
#
# How large should K be? Replaying only averages down the *emission* noise; the
# trajectory spread is fixed by the stage-1 primary count. For the halo that
# split is 69 % / 31 % of the per-primary variance, so with stage 1 at w=1:
#
#     K =  1  ->  1.00 bunch crossings of statistical power
#     K =  5  ->  2.24     (~10 % more CPU -- the sweet spot)
#     K = 10  ->  2.65
#     K -> inf ->  3.25    <-- hard ceiling, for any K
#
# S of the QD0AL sampler in the SINGLE-STAGE (full) beamline, in metres.
# Stage 2's own line is only 6.8 m long, so its Model tree reports a much
# smaller S. The converter needs the full-line value to reconstruct the time of
# the reference particle at the IP, because stage-2 particles are injected with
# their absolute stage-1 arrival time.
SAMPLER_S_FULL_M = 286.865599
