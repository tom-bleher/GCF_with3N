# Reproducing Wright et al. (arXiv:2104.05090)

Every known difference between our 2N GCF generator + GENIE chain and hers, as
of 2026-08-30. Her code is `github.com/nwright38/GCF_GENIE` @ `fc45782a`
(2021-03-29) — the commit immediately preceding the 2021-04-11 arXiv posting.

## Where we now match

| Item | Value (both) |
|---|---|
| ¹²C AV18 contacts | Cpn1 = 15.876, Cpn0 = 1.244, Cpp0 = 1.140 |
| Pair-type mix | pn 88.2 % / same-isospin 11.8 % |
| σ_cm | 0.150 GeV/c per Cartesian component |
| Off-shell struck-nucleon energy | E₁ = M_A − √(P_cm² + M_{A−2}²) − E₂, ground-state A−2 mass |
| E* | 0 (both) |
| Recoil nucleon | on-shell, p₂ = P_cm − p₁ |
| Q² floor | ≥ 1.0 GeV² |
| SRC pair vertex | ρ¹(r)·r², GENIE modified-harmonic-oscillator density (ap = 1.69, alf = 1.08) |
| Both nucleons | one shared GENIE record, same vertex |
| FSI model | hA |
| Intranuke cascade medium | `genie::LocalFGM/Correlated`, SRC-Fraction = 0.2 |

Note on the contacts: her later commit `a03f6fc1` (2021-08-21, **four months
after the paper**) changes them to Cpn1 = 46.367, Cpn0 = 3.633, Cpp0 = 50
("increase number of pp pairs"), and comments out the `AV18_pp0` table in the
same commit. That set does **not** apply to the published figures.

## Remaining differences

| Item | Hers | Ours | Assessment |
|---|---|---|---|
| **SRC onset** | k > k_F(**r**), the local Fermi momentum, i.e. AV18 tail above a flat Fermi gas | selectable: constant k_F (`fixed`) or k_F(r) (`local`) | Real modelling difference. `local` matches her and makes our SRC and mean-field samples complementary — our MF ceiling is already k_F(r), so a constant floor leaves k_F(r) < k < k_F in neither sample (78 % of vertices have k_F(r) < 0.25; median 0.212). But it extends the GCF's back-to-back pair picture down to k ≈ 0.08 GeV/c, well below the ~0.3 GeV/c where the paper itself says the GCF is validated. Kept switchable; `fixed` is the systematic. |
| **SRC trigger** | re-tests \|p₁\| = \|k + ½P_cm\| > k_F(r) *after* adding the c.m. motion | tests k_rel only (`src_require_p1` available) | Her code decides the branch twice and the two disagree near the boundary, which is what leaks pairs into her mean-field class. `fc45782a` is literally a revert of that line and `3258056a` changes it again, so she was unsure too. The GCF split belongs on the relative momentum. |
| **Pair-type sampling** | type drawn first, then \|k\| from that type's table renormalised per type | flat isospin ×4, weighted by C·φ²(k)·k² | Differs only by 1/N_type; the three integrals are 1.08 / 1.11 / 1.10, so ≤ 3 %. |
| **Elastic form factors** | BBA07 | Kelly | 1–3 % over Q² = 1–5 GeV². Every panel is area-normalised and the transparency is a ratio of two samples sharing the cross section, so it cancels. Not changed. |
| **P_cm truncation** | Gaussian truncated at ±0.5 GeV/c by the `TF1` range (3.33σ) | untruncated | Removes 0.087 % per component. Far below our statistical precision. Not changed. |
| **Vertex r_max** | 3R = 9.62 fm | 2.5R = 8.01 fm | The 2.5R–3R shell holds 9.5×10⁻⁹ of the vertex weight; even with the largest possible I_SRC enhancement (~40×) that is 4×10⁻⁷. Not changed. |
| **Transport order** | SRC recoil is added at Module-2 and the struck nucleon's daughter at Module-5, so the **recoil is transported first** | order randomised per event | Both nucleons leave the same vertex simultaneously; no order is physical. Ours avoids systematically giving one leg the fullest medium. |
| **Pauli blocking** | stock GENIE `PauliBlocker` module in the QEL chain | none at generator level | Cannot fire: it tests \|p_out\| < k_F(r) on the pre-FSI knocked-out nucleon, which carries 1–2 GeV/c under the Q² ≥ 1 GeV² floor. Intranuke's internal blocking applies in both. |
| **Weighting** | unweighted (weight ≡ 1), normalisation from the σ splines | weighted MC, normalisation from Σw / attempts | Different technique, both valid. Ours can separate the SRC and mean-field cross sections; hers cannot. |

## Two things to ask the authors

1. **Which build produced the published samples?** At `fc45782a` a
   GeV-vs-fm⁻¹ unit mismatch in `4πp²φ²` (`LocalFGM::ProbDistro`) makes the
   AV18 tail integrate to ≈0.008 instead of ≈1.1, so only ~1 % of nucleons
   land above k_F rather than the nominal 20 % — the generator is ~99 %
   mean-field and cannot be run as-is to make a GCF sample. The surrounding
   commits oscillate around exactly this (`204c1bbc`/`e5eb1c0f` zero the
   Fermi-gas part, `fc45782a` reverts that, `3258056a` adds a rejection loop),
   so the GCF sample was probably made with a GCF-only build — in which case
   the constant divides out and the published figures are unaffected. Worth
   confirming rather than assuming.

2. **Were cascade knock-out protons counted as the recoil?** Her FULL
   (e,e'pp) effective transparency is 0.73, i.e. below 1: FSI *removes* pp
   events. Ours is above 1 — FSI *adds* them — and the paper's own T+SCX
   already contains a large SCX gain (T_NN = 0.44 → 0.83). The analysis repo
   would settle whether the recoil was taken as any detected second proton or
   as the SRC partner from the GENIE record.

## Status of the reproduction

Passing:
- **Q² sanity check.** Pre-cut PWIA and FULL agree to 1.0 % mean / 3.4 % max —
  statistical only, as they must be: Q² is fixed by the scattered electron and
  our FSI never modifies the primary weight. They separate only after the
  acceptance and the SRC cuts, which is the transparency effect.
- **Fig 2 left** now shows the low-p_lead FSI shoulder the paper has. This came
  from the event-definition fix.
- **PWIA and T+SCX** ratio curves sit on the data in Fig 5.

Open:
- Effective transparencies (`fixed` / `local` onset vs paper):
  (e,e'p) 0.696 / 0.866 vs **0.58**; (e,e'pp) 1.284 / 1.166 vs **0.73**.
- Fig 5 FULL is too high with the `fixed` onset; with `local` it tracks the
  data to p_miss ≈ 0.75 and then runs away, though that region is
  statistics-limited at 2M.
- Fig 8: the mean-field c.m. distribution is not yet broader than the GCF one.
  The `LocalFGM/Correlated` cascade medium helps a great deal — yield above
  |p_cm| = 0.5 GeV/c goes from 2.5 % to 16.2 % for MF while GCF barely moves —
  but on only ~100 selected events at 2M. Needs the 60M production sample.

Composition of the selected (e,e'pp) sample after SRC cuts (`fixed` onset),
from `part_parentRole` and `part_rescatterCode`:

```
recoil is the untouched SRC partner      37.7 %
recoil is a product of the partner       45.7 %   (includes n -> p charge exchange)
recoil is a product of the struck chain  15.6 %   (medium knock-out; pure background)
```

## Validation of the |p_rel| importance sampling

`GCF_PREL_IS=0` (GCF p_rel Importance Sampling) reverts `SRC_analysis_2N` to
the old uniform |p_rel| draw with everything else identical, so the change can
be audited. 3M PWIA events each way, `fixed` onset:

| | IS on | IS off |
|---|---|---|
| N_eff / N | 5.69 % | 2.29 % |
| Σw / attempts | 8.5298 ± 0.021 | 8.4386 ± 0.032 |

**Shapes agree.** Comparing un-normalised per-attempt yields (bins independent,
so nothing can hide in a normalisation constraint) gives chi²/ndf of 1.81
(k_rel), 1.11 (p_miss), 1.90 (Q²), 0.60 (x_B), with no trend across bins — the
pulls are a common positive offset, not a slope. Splitting each sample in half
and comparing like with like shows this statistic naturally scatters over
0.64–1.83, so those values are unremarkable.

**Integrals differ by 1.08 % (2.38σ).** Bootstrap errors match the analytic
√Σw² for both samples (ratio 1.00–1.02), and the running estimate does not
drift upward with sample size, so this is not an unsampled-tail artefact. The
estimator is unbiased algebraically: drawing k with density f_hat/dF and
weighting by dF·(f/f_hat) gives ∫f dk exactly, whatever the accuracy of the
tabulated cumulative. 2.4σ across a handful of observables is not compelling,
but it is not fully settled either.

**It cannot affect any published number.** Every figure panel is area-normalised
to the data, and the effective transparency is a ratio of two samples generated
with the same setting, so a common normalisation offset cancels identically.
The setting must simply be kept the same within a PWIA/FSI pair — which
`run_natalie_samples.sh` guarantees.

