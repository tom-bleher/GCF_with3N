Systematic variant of config/genie_override: identical except the Intranuke
cascade-target NuclearModel is genie::LocalFGM/Correlated (local Fermi gas
with a 20% SRC high-momentum tail, SRC-Fraction=0.2) instead of the
baseline genie::LocalFGM/Default (no tail, matching N. Wright's published
setup, used for all paper figures).

The tail is arguably more physical for the cascade medium -- the residual
nucleus still contains SRC pairs -- but it only affects rescattering
kinematics and in-cascade Pauli blocking, a second-order effect on the
paper's shape observables. Kept as a documented systematic:

    GENIE_OVERRIDE_DIR=$PWD/config/genie_override_srctail ./run_2N_parallel.sh ...
