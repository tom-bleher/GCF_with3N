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

---

# Her analysis code (`GENIE-GCF-Analysis`)

All items verified directly against the source; file:line given.

## Which file is which

1. **`genie_gcf_loop.C` produced the paper's FULL curve.** Its output path
   (`:628`, `..._501_1bil_reco_f.root`) is verbatim the file
   `gst_plotter.C:1309` opens as the FULL input. `_501_` = 5.014 GeV beam,
   tune G18_10a (hA).
2. `gcf_generator_analysis_EM.cpp` = PWIA baseline, reads the standalone GCF
   generator tree `genT`. No GENIE.
3. `gcf_genie_analysis.C` = an **older** GENIE analysis with a completely
   different sample definition (lead/recoil strictly by daughter-link
   provenance, `:1248`, so FSI knock-outs are excluded by construction).
   Superseded -- do not confuse with item 1.
4. `gst_plotter.C` = plotting and normalisation.

## Event building (`genie_gcf_loop.C`)

5. The final-state list contains **protons only** (`:1506`); each is smeared
   1% on the magnitude, then must pass the proton fiducial **and** have a
   non-zero acceptance-map value -- the map is a **binary filter** here, not a
   weight (`:1514-1516`).
6. **The lead is simply the highest-momentum accepted proton, with no cuts.**
   Her first loop intends to apply the SRC cuts to lead candidates, but
   `:1618` compares against `pLead_max`, a histogram axis limit of 2.4 GeV
   (`:701`, active since `data_bins = true` at `:641`), instead of the running
   `lead_max`. The loop almost never fires; the fallback (`:1628-1640`)
   decides.
7. **The recoil is any other detected proton above 0.35 GeV/c**, falling back
   to any other proton if none qualify (`:1652-1687`). Provenance is never
   enforced -- a cascade knock-out is as valid a recoil as the SRC partner.
8. **With 2 or 3 candidates she runs a lottery** (`assign_recoil`,
   `:540-623`): each candidate's map value is an independent detection
   probability, one random number is thrown, and **if two or more come out
   detected the entire event is discarded, its (e,e'p) entry included**
   (`:1728-1731`). With 4+ the dispatch is commented out and the event
   silently becomes (e,e'p)-only.
9. The recoil map weight is applied **only** when there was exactly one
   candidate (`:1691`).
10. **"Goodbye Absorption"** (`:1554-1557`): any event where either SRC
    nucleon has rescatter fate 4 is discarded (27% of our weighted sample).
11. **Analysis-level truth cut, k_rel > 0.25 GeV/c** (`:1824-1828`). `pn`
    there is the struck nucleon's **truth** momentum (status-11 particle,
    `:1476-1479`), reassigned to the smeared p_miss only on the line *after*
    the cut. For the mean-field sample `true_recoil` is the outgoing nucleon,
    so her prel is ~0.5|q| and the cut never fires there.
12. Fate legend (`:1541-1552`): 0 Undefined, 1 No Interaction, 2 SCX,
    3 Inelastic, 4 Absorption, 5 Knock out, 6 Compound nucleus,
    7 Pion Production. Confirmed independently: fate 1 was the only value
    whose final-state particle matched its input nucleon's momentum.

## Cuts and constants

13. `pass_selection_cuts` (`:512-524`): xB > 1.2, theta_pq < 25 deg,
    **0.62 < |p_N|/|q| < 0.96**, 0.4 < p_miss < 1.0, m_miss < 1.1,
    p_recoil > 0.35 (epp only). **The 0.96 disagrees with the paper's text,
    which says 0.92.**
14. **No Q^2 cut anywhere** -- Q^2 appears only as a histogram axis.
15. Smearing: protons 1%, electrons 1.25%, magnitude only (`:637-638`).
16. SCX probabilities (`:422-436`) identical to ours; applied as a per-event
    random flip rather than a weight.

## Normalisation

17. **The published (e,e'pp) T+SCX transparency 0.83 is hard-coded**, not
    simulated: `gst_plotter.C:1356`, `transparency = (.44/.53)` = 0.8302
    = T_NN/T_N.
18. Every curve is area-normalised to the data on the (e,e'p) p_miss
    histogram (`gst_plotter.C:1360`), so the (e,e'p) yield carries no
    absolute information.
19. The GENIE event weight is never used; the only per-event weight is the
    acceptance-map value.

## Bugs in her code that affect the published numbers

20. `pLead_max` vs `lead_max` (`:1618`) -- see item 6.
21. **`events_pass_ep` over-counted ~7x**: the increment sits in an `else`
    attached to the m_miss test *inside* the cut ladder (`:1971-1976`), so it
    fires once per surviving stage. `events_pass_epp` fires once. Any
    transparency from those printouts is wrong, by different factors for ep
    and epp.
22. **Out-of-bounds read**: `norm_factors` has 3 elements
    (`gst_plotter.C:1360`) but `:1172-1181` indexes `norm_factors[j+6]` for
    j = 0,1,2 -- undefined behaviour in the Fig 6 / Fig 8 scaling, alongside
    hand-tuned rebin factors.
23. `p0_id = 1;` immediately followed by `if(p0_id == -3)` (`:1486-1490`)
    makes her p0_id/p1_id provenance indices dead.
24. `gRandom->SetSeed(i)` re-seeds per event when SCX is on (`:1860`), after
    the smearing was drawn.
25. phi-frame inconsistency: protons rotated by +pi before the fiducial and
    map filter (`:1514`), but the lead's map *weight* uses the unsmeared,
    unrotated vector (`:1853`).
26. `Sumw2()` is never called in either GENIE macro, so her weighted GENIE
    histograms carry sqrt(N) errors.

## Reproduction with her conventions

`--convention wright`, `local` SRC onset. **Corrected 2026-08-31** after the
60M production: the earlier 2M numbers were not converged.

| statistics | (e,e'p) Full | (e,e'pp) Full |
|---|---|---|
| 2M validation sample | 0.594 +- 0.053 | 0.740 +- 0.152 |
| first 2M of the 60M sample | 0.607 +- 0.055 | 0.920 +- 0.160 |
| 30M of the 60M sample | 0.612 +- 0.015 | **1.100 +- 0.059** |
| paper | 0.58 | 0.73 |

The quoted errors are 1/sqrt(N_eff) added in quadrature over the PWIA and FSI
weight sums, N_eff = (Sum w)^2 / Sum w^2 computed on the *selected* bucket
(inside the acceptance, after the SRC cuts). They are the honest errors on this
estimator and had not been quoted before.

**The (e,e'pp) transparency is weight-tail limited.** In the selected (e,e'pp)
FSI bucket at 2M there are ~800 surviving events but only N_eff ~ 35-60: the
single largest-weight event carries 4-10 % of the total and the top ten carry
34-43 %. So the 2M value has a ~20 % error, and the apparent agreement with her
0.73 was a ~2 sigma downward fluctuation, not a reproduction. At 30M the bucket
has 11941 events with N_eff = 533 (top event 1.4 %, top ten 8.2 %) and the value
is stable: 1.06 (2M subset) -> 1.14 (8M) -> 1.14 (30M).

Conclusion: **(e,e'p) is reproduced** (0.61 vs 0.58); **(e,e'pp) is not**
(1.10 +- 0.06 vs 0.73, a 6 sigma gap). FSI still adds pp pairs for us and
removes them for her. The cascade medium is not the cause: srctail 1.186 vs
no-tail 1.138 at 30M.

The `local` onset still matters for a physical reason unrelated to the analysis
cut k_rel > 0.25: with it the pair-vertex distribution is rho^n(r) r^2 I_SRC(r),
which is surface-biased because I_SRC(r) grows as kF(r) -> 0, and surface pairs
suffer less FSI. That is a genuine consequence of her "AV18 tail above a local
Fermi gas" ansatz, not a normalisation choice.

Note for future runs: N_eff, not the raw selected-event count, is what sets the
precision of any weighted transparency here. Quote it.

## Where the (e,e'pp) excess comes from (2026-08-31)

Fate-by-fate decomposition of the selected (e,e'pp) FSI yield, `--convention
wright`, no-tail medium, 8M events (3245 selected, so the fractions carry a few
% relative from the weight tail). `rec_role` = which chain the detected recoil
descends from, `rec_fate` = the hA fate stamped on its ancestor:

| recoil origin | fate | yield share |
|---|---|---|
| partner chain | NoInteraction (the true SRC partner) | 47.5 % |
| partner chain | **ChargeExchange (n -> p)** | **39.3 %** |
| partner chain | Inelastic | 5.9 % |
| struck chain | Inelastic | 4.5 % |
| struck chain | PiProd | 1.7 % |
| struck chain | NoInteraction | 1.2 % |

Applying those shares to the 30M transparency 1.098:

| component | contribution | Glauber reference |
|---|---|---|
| genuine pp pair survives | 0.52 | T_NN = 0.44 |
| n -> p charge exchange gain | **0.43** | none -- T+SCX cannot add pairs |
| other cascade products | 0.15 | none |
| total | 1.10 | paper 0.73 |

Two conclusions:

1. **Our cascade is only ~18 % more transparent than the Glauber numbers the
   paper uses** (pp survival 0.52 vs T_NN 0.44), consistent with the (e,e'p)
   channel (0.613 vs T_N 0.53, +16 %). That part is a modest, defensible
   difference.
2. **The excess is n -> p charge exchange of the SRC partner.** GCF is 88.2 %
   pn / 5.9 % pp, so the partner is almost always a neutron; hA converts it to
   a detected proton often enough to promote pn pairs into the (e,e'pp)
   channel. A Glauber T+SCX treatment cannot do this -- it only removes flux
   and flips isospin as a substitution, never as a gain.

The implied per-pair conversion rate is not obviously wrong: CX yield /
untouched-pp yield = 39.3/47.5 = 0.83, and scaling by the pair-type mix gives
0.83 x (5.9/88.2) ~ 0.055 per pn pair, against the paper's Glauber SCX
probabilities of 0.035-0.048. Same order.

But it does not fully close the gap: even rescaling our pp survival down to the
Glauber T_NN would give 0.44/0.52 x 1.10 = 0.93, still well above 0.73. So her
hA must produce far less charge-exchange gain in this channel than ours does,
despite both being GENIE hA. **That is now the sharpest open question.**

Also measured (30M, her conventions): accepted protons per event

| protons | PWIA | FSI |
|---|---|---|
| 0 | 19964034 | 16888948 |
| 1 |  9110708 | 10358219 |
| 2 |   925258 |  2376437 |
| >=3 |      0 |   376396 |

PWIA cannot exceed 2 (the pair has two nucleons) and indeed is exactly 0. FSI
gives 2.57x more two-proton events. That is proton multiplication.

## src_require_p1 (Wright's second branch test) -- measured, no effect

30M each, no-tail medium, everything else pinned to the baseline:

| | (e,e'p) | (e,e'pp) |
|---|---|---|
| baseline, wright | 0.613 | 1.138 |
| src_require_p1=1, wright | 0.633 | 1.098 |
| baseline, detector | 0.834 | 1.678 |
| src_require_p1=1, detector | 0.838 | 1.546 |

The switch removes 7.6 % of the SRC strength (Sum w / attempts 39.5 vs 42.7)
and moves the transparencies within their errors. Lever closed.

## 2026-09-13 — revised detector convention and the {absorption} x {recoil rule} grid

Defaults of `natalie_paper_plots.py --convention detector` (all switchable):

- CLAS fiducials and acceptance map evaluated at **phi + 180 deg**
  (`--acc-phi-rotation`), as `genie_gcf_loop.C` does for the electron
  (1425-1426) and every proton (1514). The fiducial parameters and proton gaps
  are per sector, so this is not a symmetry: it changes the fiducial verdict
  for 13 % of protons and 5 % of electrons.
- **Absorbed events thrown away** (weight 0) when either generated nucleon has
  hA fate 4 (`--absorption-veto`, default; `--no-absorption-veto` keeps them).
  Her code vetoes absorption everywhere.
- (e,e'p): no recoil requirement.
- (e,e'pp) recoil rule `--epp-recoil`:
  `unique` (exactly one other accepted proton >= 0.35 GeV/c, else not (e,e'pp)
  but still (e,e'p)); `all` (each such proton gives its own (e,e'pp) entry);
  `exactly2` (exactly two accepted protons in total -- the topology of the CLAS
  data file, whose 605 events all have two protons); `highest` (old rule).
  The CLAS papers (Cohen PRL 121 092501; Duer PRL 122 172502) define the recoil
  as a second proton >= 350 MeV/c and publish no multi-recoil rule.
- `--convention wright` now applies her rotation as well: gate rotated, while
  the lead weight (1853) and the lottery (1697) read the map at the unsmeared,
  unrotated momentum, as her code does.

Full effective transparencies, 30M SRC / 60M MF, no-tail medium:

| recoil rule | absorbed events | (e,e'p) | (e,e'pp) |
|---|---|---|---|
| unique | thrown away | 0.684 | 1.286 |
| all | thrown away | 0.684 | 1.304 |
| unique | kept | 0.829 | 1.653 |
| all | kept | 0.829 | 1.759 |
| paper | | 0.58 | 0.73 |

**Multi-recoil events come almost entirely from absorption.** Weighted fraction
of SRC-cut (e,e'p) events with >= 2 possible recoils: GCF 0.07 % (absorbed
thrown away) vs 0.32 % (kept); MF 0 % vs 0.92 %. With absorbed events thrown
away the recoil rule is irrelevant (unique vs all: 1.4 %).

**Fig 8.** hA's nucleon "inelastic" fate is a 2->2 collision with a partner
drawn from the Fermi distribution (`HAIntranuke2018::InelasticHA`,
`INUKE-DoFermi = true`), so for a single rescatter MF p_cm = p_i + p_t, the sum
of two initial Fermi momenta: narrow by momentum conservation. Per recoil-fate
class (unique rule, 8M GCF / 30M MF):

| sample | recoil from | share | <\|p_cm\|> | frac > 0.5 |
|---|---|---|---|---|
| MF | struck: inelastic | 42 % | 0.239 | 0.0 % |
| MF | struck: pion production | 33 % | 0.257 | 2.7 % |
| MF | struck: absorption | 25 % | 0.486 | 44.2 % |
| GCF | partner untouched | 34 % | 0.280 | 7.5 % |
| GCF | partner charge exchange | 26 % | 0.371 | 18.4 % |
| GCF | absorption (partner / struck) | 8 / 8 % | 0.50 / 0.44 | 41 / 39 % |

Absorption is the only channel with a broad p_cm. With absorbed events thrown
away (her convention) MF has 1.2 % above |p_cm| = 0.5 against 11.5 % for GCF:
**MF is narrower than GCF, the opposite of the paper's statement.** Keeping
absorbed events makes MF about as broad as GCF, never broader. Her Fig 8/9/10
plotting code is not in her repository (`gst_plotter.C` has no mean-field
input), so her MF curve cannot be checked against code.

Correction to an earlier statement in this file's history: the back-to-back
MF correlation is not produced by the SRC cuts; it follows from momentum
conservation in the 2->2 knock-out.
