Variant of config/genie_override: identical except that the Intranuke
cascade-target NuclearModel is genie::LocalFGM/Correlated (local Fermi gas
with a 20% SRC high-momentum tail, SRC-Fraction=0.2) instead of the
baseline genie::LocalFGM/Default (no tail).

    GENIE_OVERRIDE_DIR=$PWD/config/genie_override_srctail ./run_2N_parallel.sh ...

WHICH ONE MATCHES N. WRIGHT?  Correlated -- i.e. THIS directory.

An earlier version of this file claimed the opposite. That was wrong.
In her fork (github.com/nwright38/GCF_GENIE @ fc45782a),
config/G18_10a/ModelConfiguration.xml:37 reads

    <param type="alg" name="NuclearModel"> genie::LocalFGM/Correlated </param>

with the LocalFGM/Default line on the next line commented out, and the
per-nucleus SpectralFunc1d overrides sitting inside a comment block. So
C-12 does use LocalFGM/Correlated, and config/LocalFGM.xml:349 gives that
param_set SRC-Fraction = 0.2.

This is the model Intranuke draws its cascade *target* nucleons from, so it
sets the momentum an FSI knock-out nucleon can carry. That matters most for
the mean-field sample, where the second detected proton can ONLY be an FSI
product: without the tail those knock-outs are capped near the local Fermi
momentum, with it they reach ~1 GeV/c. Fig 8 of the paper (Wright et al.,
arXiv:2104.05090, doi:10.48550/arXiv.2104.05090) turns on the MF
pair c.m. distribution being much BROADER than the GCF one, so this is a
first-order setting for that figure, not the second-order systematic the
old note assumed.

run_natalie_samples.sh therefore points at this directory. The baseline
config/genie_override (no tail) is kept as the systematic variant, and is
still what the 3N paper's samples use.
