# 3N wavefunction comparison — new Faddeev tables vs pD-only vs old AV8

Run 2026-10-05, after merging `origin/master` (`aeb725f` + `008f78a`) into local master
(merge commit `ca2a240`).

## What changed upstream

* `aeb725f` replaced `array/{AV8,N2LO,G3,AV4}_1D.dnst` with five Faddeev 3He densities
  (`u=1` AV18, `u=2` AV18+UIX, `u=3` AV18+TM', `u=4` CD-Bonn, `u=5` CD-Bonn+TM'),
  stored as the RAW unpolarized density rho(theta_ab, k_cm, k_rel).
* `008f78a` added the pD-only (deuteron-channel) counterparts as `u=6..10`.
* Only non-comment code change in `QEGenerator_3N.cc`: the table filenames in
  `fill_array()` (plus an unrelated `fabs()` on the delta-function Jacobian).
  The phase-space Jacobian `J = p_a^2 p_b^2 sin(theta_ab)` was ALREADY in place at
  `8df93bc`, so old-vs-new differences here are a TABLE effect, not a Jacobian effect.

## Conventions (verified in `get_rho()`)

Sampling builds the triplet as `v1 = v_a`, `v2 = -v_a/2 + v_b`, `v3 = -v_a/2 - v_b`,
so P_cm = 0 and `p_a = p1`, `p_b = (p2-p3)/2`, theta_ab = angle(p_a, p_b).
Before the table lookup the code PERMUTES so the **neutron** occupies the a-slot and the
**two protons** form the pair (selected by isospin, not by index); that map has |det| = 1
and needs no Jacobian. So the table axes are:
`k_cm` = neutron momentum, `k_rel` = half the pp relative momentum.

Note for anyone reading the trees: `p1_gen` (`pLead_pre` in the old code) is the OUTGOING
lead, already carrying q. The initial struck momentum is `p1_gen - q` — identical to the
p_miss vector. All wavefunction variables below use the initial momenta.

## Runs

6 GeV, 12C, kF = 0.22 GeV/c, sigma_CM = 0.15, Q^2 > 1 GeV^2; FSI = hN, independent transport.
PWIA (`-n`): 10M events each. FSI: 2M events each. AV8 baseline generated with the
pre-replacement build (commit `8df93bc`) pointed at `AV8_1D.dnst`.

## 1. Table level (`tables_summary.txt`, `tables_marginals.png`, `tables_tails.png`)

| table | pD/full | f(p_b>0.22 GeV) | f(p_b>0.25 GeV) | <p_b> fm^-1 |
|---|---|---|---|---|
| AV18 (u=1)        | —     | 4.67 % | 3.78 % | 0.523 |
| AV18+UIX (u=2)    | —     | 5.47 % | 4.41 % | 0.552 |
| CD-Bonn (u=4)     | —     | 3.05 % | 2.19 % | 0.504 |
| pD AV18 (u=6)     | 0.531 | 6.40 % | 5.25 % | 0.529 |
| pD AV18+UIX (u=7) | 0.534 | 7.33 % | 5.99 % | 0.561 |
| pD CD-Bonn (u=9)  | 0.519 | 3.94 % | 2.94 % | 0.503 |
| AV8 (old)         | —     | 3.24 % | 2.60 % | 0.503 |

* The three full tables agree in total normalisation to 1e-4 — they are consistently normalised.
* pD carries ~52-53 % of the full normalisation (matches the commit message) but a LARGER
  fraction of its strength at high p_b: pD is harder than the full wavefunction.
* AV8's integral is 4e-4 of AV18's — meaningless, its convention was never pinned down.
  Treat AV8 as shape-only.

## 2. Yields — sum(weight) per generation attempt (`samples_summary.txt`)

| sample | PWIA ratio to AV18 | FSI ratio to AV18 |
|---|---|---|
| AV18 (u=1)       | 1.000 | 1.000 |
| CD-Bonn (u=4)    | 0.863 | 0.843 |
| pD AV18 (u=6)    | 0.764 | 0.727 |
| pD CD-Bonn (u=9) | 0.682 | 0.670 |
| AV8 (old)        | 3e-4 (scale meaningless) | 3e-4 |

**Key number: pD AV18 gives 76 % of the AV18 yield, not the 53 % of its table
normalisation.** The SRC selection (k > kF, Q^2 > 1) sits in the high-momentum region where
the deuteron channel carries a disproportionate share — exactly what the table tail
fractions predict. FSI does not change this much (0.76 -> 0.73).

## 3. Shapes, pre-FSI (`pwia_observables.png`, `pwia_tails.png`, `pwia_shape_ratios.png`)

| sample | <p1_init> | <p_b> | <p_a> | f(p1>0.4) |
|---|---|---|---|---|
| AV18 (u=1)       | 0.4396 | 0.3507 | 0.4584 | 0.459 |
| CD-Bonn (u=4)    | 0.3754 | 0.2999 | 0.3910 | 0.313 |
| pD AV18 (u=6)    | 0.4166 | 0.3142 | 0.4693 | 0.420 |
| pD CD-Bonn (u=9) | 0.3703 | 0.2844 | 0.4065 | 0.322 |
| AV8 (old)        | 0.4210 | 0.3486 | 0.4208 | 0.403 |

* **Potential choice matters more than full-vs-pD.** AV18 -> CD-Bonn moves <p1> by
  -0.064 GeV/c and f(p1>0.4) from 0.46 to 0.31; AV18 -> pD AV18 moves them by only
  -0.023 and 0.46 -> 0.42.
* **pD redistributes between the two Jacobi axes:** pair relative momentum DOWN
  (<p_b> 0.351 -> 0.314) while the spectator-neutron momentum goes UP
  (<p_a> 0.458 -> 0.469). The full wavefunction puts more strength in the pp pair, the
  pD channel more in the neutron.
* theta_ab is ~90 deg in every table (the angular distribution is symmetric about 90 deg).
* **Old AV8 sits close to AV18 in shape** (<p1> 0.421 vs 0.440, <p_b> 0.349 vs 0.351) but
  with a visibly lower <p_a> (0.421 vs 0.458). So switching tables changes the SRC yield
  scale a lot and the shapes only moderately.

## 4. Post-FSI (`fsi_observables.png`, `fsi_tails.png`, `fsi_shape_ratios.png`)

| sample | <n protons > 0.35> | <n nucleons > 0.35> | <n particles> |
|---|---|---|---|
| AV18 (u=1)       | 1.629 | 2.786 | 6.72 |
| CD-Bonn (u=4)    | 1.441 | 2.485 | 6.80 |
| pD AV18 (u=6)    | 1.541 | 2.752 | 6.75 |
| pD CD-Bonn (u=9) | 1.471 | 2.550 | 6.88 |

FSI shifts the pre-FSI means by <= 0.01 GeV/c in every case, i.e. the wavefunction
differences survive the cascade essentially unchanged. The hard-nucleon multiplicity
tracks the hardness of the initial state (AV18 highest, CD-Bonn lowest), as expected.

## Bottom line

1. The pD tables are NOT a 53 % rescale of the full ones in the SRC regime — they are
   harder, and deliver ~73-76 % of the yield.
2. The dominant systematic between these wavefunctions is the POTENTIAL (AV18 vs CD-Bonn,
   ~14-16 % in yield and a clear shape shift), not the pD channel selection.
3. Old AV8 is shape-wise close to the new AV18 table, so previous 3N shape results should
   not move much; absolute normalisation is not comparable and was the reason for the
   replacement.
4. FSI is essentially transparent to the choice: it rescales yields by a few percent and
   leaves shapes alone.

## Files

* `compare_tables.py` / `tables_*.png` / `tables_summary.txt` — density tables, no generator.
* `compare_samples.py` / `{pwia,fsi}_*.png` / `samples_summary.txt` — generated samples.
* Samples themselves: `../../events/wf_test/{pwia,fsi}_u{1,4,6,9}.root`, `*_av8.root`.
