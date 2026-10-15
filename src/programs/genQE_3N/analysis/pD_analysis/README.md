# pD analysis

**Target is ¹²C.** `QEGenerator_3N.cc:49-50` calls `SetTargetNucleus(12, 6)`;
ejecting a ppn triplet leaves A−3 = ⁹Be. Confirmed from the samples themselves:
`E_lead + E₂ + E₃ − ν = 2.78172` GeV against `M(¹²C) − M(⁹Be) = 2.78214` (the
0.4 MeV gap is E*). The density tables are **³He** 3N densities, used as the
universal 3N SRC density exactly the way GCF uses the deuteron for 2N pairs —
so "³He" below refers only to the table, never to the nucleus being scattered off.

## Two different "proton vs pn pair" angles

With **a** = the odd proton's momentum, there are two candidates for **b**:

| **b** | what it is | status |
|---|---|---|
| `(p_p - p_n)/2` | the pair's **relative** momentum | carries wavefunction information — observable 1 |
| `p_p + p_n` | the pair's **total** momentum | **identically 180°** when the three momenta sum to zero — observable 2 |

The second is fixed by momentum conservation: in the **3N rest frame**
`p1+p2+p3 = 0`, so the pair's total momentum is exactly `-p_odd`. Verified with
`use_CM` off (`pwia_u1_noCM.root`): min = max = mean = 180.0000°, every event,
every potential. It becomes non-trivial only when the triplet is *not* at rest —
which in ¹²C it genuinely is not (the 3N cluster moves and recoils against ⁹Be)
— and when the cascade scrambles the final state. So it is a CM/FSI diagnostic,
not a wavefunction observable.

Observable 1: the opening angle between the **odd proton** and the **relative
momentum of the remaining (p,n) pair**, computed three ways — from the density
table alone, from PWIA events, and from the full generator with FSI — for every
potential, full wavefunction and pD component only.

## Can you choose "the" proton?

The ejected triplet is ppn. Picking one proton as the odd nucleon always leaves a (p,n) pair,
so the decomposition exists; the question is which proton.

**From the wavefunction: no, and provably so.** The two protons are identical.
In the table's own Jacobi set the two protons sit in slots 2 and 3, and swapping
them is exactly `theta -> 180 - theta`. Every table is symmetric under that to
~1e-6 (storage rounding), so the two choices give the *same* distribution. The
A/B overlay in column 3 of `opening_angle_density.png` confirms this
numerically: `max|A-B|/max(A) = 3e-9`.

**From events: yes.** The electron picks the struck nucleon through
`sigma_eN` and the kinematics, which is not symmetric between the two protons.
So "odd = struck proton" and "odd = recoil proton" are genuinely different
measurements, and the size of the difference is physics, not labelling.

Note the pair here is (p,n) — *not* identical particles — so its relative
momentum has a well defined sign. We use `(p_proton - p_neutron)/2` throughout.

## Conventions

The tables are stored in the Jacobi set `QEGenerator_3N::get_rho()` permutes
into (`QEGenerator_3N.cc:332-400`): slot a = the **neutron**, slots 2,3 = the
two **protons**, so the stored axes are `k_a = |p_n|`,
`k_b = |p_p1 - p_p2|/2`, `theta = angle(k_a, k_b)`. For the proton-odd
decomposition the three momenta are rebuilt exactly as the generator does
(`QEGenerator_3N.cc:205-215`, P_cm = 0):

    v1 = (0,0,k_a)              the neutron
    vb = k_b (sin th, 0, cos th)
    v2 = -v1/2 + vb             proton 1
    v3 = -v1/2 - vb             proton 2
    choice A: odd = v2, b_A = (v3 - v1)/2      choice B: odd = v3, b_B = (v2 - v1)/2

The angle is rotation invariant, so it is a function of `(k_a, k_b, theta)`
alone — no Euler-angle sampling needed. Each grid point carries the Cartesian
weight `rho * k_a^2 k_b^2 sin(theta)`.

**Generator support.** `fkMinInternal = 0.25` (`QEGenerator_3N.cc:27`) makes the
generator discard any configuration with *any* of the three internal momenta
below 0.25 GeV/c (`QEGenerator_3N.cc:224`, before the CM boost). That is a
triple cut, not a cut on one nucleon, and it changes the distribution
completely — the broad plateau of the full wavefunction becomes a peak at 90°.
Row 2 of the density figure applies it, and that is the row to compare with
events.

**Events.** PWIA uses the generated initial momenta exactly as they would be
measured: `p_miss = p1_gen - q` (`p1_gen` is the *outgoing* lead, so q must be
removed), `N1_type = proton` required. No CM correction is applied. These
samples were generated with `use_CM` on, `sigCM = 0.15` GeV/c per component
(`QEGenerator_3N.cc:26`, the literature ¹²C value) giving `<|v_cm|> = 0.24`
GeV/c. **This CM is physical**: the 3N cluster moves inside ¹²C and recoils
against the ⁹Be remnant. Removing it via `p_i -> p_i - P_tot/3` shifts `<theta>`
by only 0.3°, and it cannot be done on the FSI side at all (`P_tot` is not
measurable), so both rows use the raw momenta. `pwia_u1_noCM.root` (CM switched
off, an artificial control) gives 82.4° against 83.0° for `pwia_u1`. FSI uses the measured final state from the
post-cascade list: detected = p or n with |p| > 0.25, lead = detected nucleon
closest in angle to q (required to be a proton), recoil p/n = highest-|p|
detected proton/neutron among the rest, `p_miss = p_lead - q`. No CM removal is
possible there. All samples: 8 < theta_e < 45 deg, Q² > 1 GeV².

## What the plots show

- **PWIA reproduces the wavefunction.** On the generator support the means agree
  to ≲1.5° (AV18 82.7° events vs 84.0° table; pD AV18 76.9° vs 77.6°).
- **The two choices agree in PWIA** (≈1.5° apart) but **split under FSI**
  (AV18: 87.8° struck vs 97.7° recoil). The cascade and the detection threshold
  treat the struck and recoil protons differently.
- **CM motion.** The samples carry a 0.15 GeV/c-per-component CM, the physical
  motion of the 3N cluster in ¹²C. It barely touches this observable (0.3°), but
  it dominates any observable sensitive to the three momenta summing to zero —
  see observable 2. `pwia_u1_noCM.root` is the CM-off control.
- **FSI washes the distribution out**: the 90° peak is suppressed and the
  back-angle region filled in.
- **pD-only sits at smaller angles** than the full wavefunction throughout
  (74.5° vs 82.6° on full support; 77.6° vs 84.0° on the generator support).
- **Three-body forces matter more than the NN potential here** on the generator
  support: AV18 84.0°, AV18+UIX 83.6°, AV18+TM 87.0°, CD-Bonn 84.6°,
  CD-Bonn+TM 89.3°.

## Normalisation

Every curve is normalised **individually to unit area** over its plotted range
(`h / sum(h) / bin_width`), with the event weights as histogram weights — so the
y axis is a weighted probability density and `sum(h)` is `sum(w)`, not a count.
The two curves in a panel, and the grey wavefunction reference, are each
normalised separately: **the plots compare shapes only** and carry no
information about relative rates between pairings or potentials.

Entries outside the plotted range are dropped from both numerator and `sum(h)`,
so strictly it is unit area *within the window*. Measured fractions outside:

| observable | PWIA | FSI |
|---|---|---|
| `angle(p_p,p_n)`, 0–180° | 0.00% | 0.00% |
| `prel`, 0–1.2 GeV/c, odd = struck p | ≤0.01% | 0.6–1.0% |
| `prel`, odd = recoil p | ≤0.05% | 4.2–5.5% |
| `M_miss − M(⁹Be)`, 0–1.5 GeV | 0.00% | 3.0–3.4% |

The quoted means are computed over the **full** sample, not the window, so they
are not biased by that truncation (and may sit slightly off the drawn curve's
own mean for the FSI cases). `pair_observables_summary.txt` carries the
out-of-range fraction per curve.

## Caveat: effective statistics

The event samples are weight-dominated. Raw entries are millions but
`N_eff = (sum w)^2 / sum w^2` is **1.0k–2.5k for the FSI samples** and
**5.7k–14.5k for PWIA** (see `opening_angle_events_summary.txt`). The density
curves are exact (no Monte Carlo), the PWIA means are sound, but the FSI shapes
are visibly noisy and potential-to-potential differences in the FSI row should
not be read as physics without more events.

## Observable 2: angle to the pair's TOTAL momentum

`pair_total_angle.py` → `pair_total_angle.png`, `pair_total_angle_summary.txt`.

- **CM-off control**: exactly 180.0000° in every event (the artificial limit).
- **PWIA**: 145.7°–148.9°, all four potentials nearly identical and the two
  choices of odd proton agreeing to <0.5°. This width is set entirely by
  `|v_cm| / |p_odd|`, i.e. it measures the 3N cluster's motion in ¹²C; the small
  spread between potentials just reflects how hard the odd proton is, not the
  shape of the wavefunction.
- **FSI**: pushed much further from 180° (struck 134.9°–138.6°, recoil
  127.3°–130.5°), and the struck/recoil split reappears at ~8–10°, as in
  observable 1.

## Observables 3–5: the (p,n) pair itself

`pair_observables.py` → `pair_relmom.png`, `pair_opening_angle.png`,
`missing_mass_pd.png`, `pair_observables_summary.txt`.

**3. Relative momentum `|(p_p − p_n)/2|`.** Wavefunction (generator support):
AV18 0.452, CD-Bonn 0.382, pD AV18 0.435, pD CD-Bonn 0.367 GeV/c. PWIA events:
0.37, 0.32, 0.37, 0.32. The ordering is reproduced but the events sit ~20%
*below* the table — unlike the opening angles, which agreed to a few degrees.
That is expected in direction: the grey curve is the bare density, while the
events carry `sigma_eN`, the δ-function Jacobian and the θ_e/Q² cuts, all of
which fall with momentum and so pull the mean down. It does mean the
"PWIA reproduces the wavefunction" statement holds for *angular* observables and
only qualitatively for *magnitudes*. FSI pushes it up (0.39–0.54).

**4. Opening angle `angle(p_p, p_n)`.** Wavefunction 122.7°–125.1°; PWIA
119.4°–122.9° (agrees to ~3°); FSI splits by pairing — 102.6°–107.2° for
odd = struck p against 115.6°–118.8° for odd = recoil p.

**5. Missing mass of ¹²C(e,e'pd)X**, with
`p_miss = p_lead + (p_p+p_n) − q`, `E_miss = ν + M(¹²C) − E_lead − E_d`,
`E_d = sqrt(M_d² + |p_pair|²)`. The undetected system is ⁹Be, so the
distribution has a hard lower edge at `M(⁹Be) = 8.3928` GeV and we plot the
excess above it — the energy cost of calling the free (p,n) pair a deuteron,
zero for a perfect deuteron. Only one pairing is meaningful (the pair is the two
recoils), so there is no struck/recoil choice here.

| sample | ⟨M_miss − M(⁹Be)⟩ | fraction < 50 MeV |
|---|---|---|
| PWIA AV18 | 0.161 GeV | 16.7% |
| PWIA CD-Bonn | 0.121 GeV | 19.8% |
| PWIA pD AV18 | 0.158 GeV | 15.9% |
| PWIA pD CD-Bonn | 0.123 GeV | 17.3% |
| FSI AV18 | 0.367 GeV | 9.0% |
| FSI CD-Bonn | 0.330 GeV | 10.3% |
| FSI pD AV18 | 0.363 GeV | 8.2% |
| FSI pD CD-Bonn | 0.309 GeV | 9.5% |

Two things stand out. **The pD tables are no more deuteron-like than the full
ones** (AV18 0.161 vs pD AV18 0.158; CD-Bonn 0.121 vs pD CD-Bonn 0.123) — the
pD selection does not concentrate the pair near the deuteron on this observable.
**The NN potential does matter**: CD-Bonn is consistently ~25% closer to the
deuteron than AV18, full or pD. This repeats the pattern from the earlier
wavefunction comparison, where the potential mattered more than the pD choice.
FSI roughly doubles the excess and halves the deuteron-like fraction.

## Files

| file | what |
|---|---|
| `pd_jacobi.py` | shared helper: table loading, Jacobi transformation, cell-integrated histograms |
| `opening_angle_from_density.py` | → `opening_angle_density.png`, `opening_angle_density_summary.txt` (all 10 tables) |
| `opening_angle_from_events.py` | → `opening_angle_events.png`, `opening_angle_events_summary.txt` (u=1,4,6,9) |
| `pair_total_angle.py` | observable 2 → `pair_total_angle.png`, `pair_total_angle_summary.txt` |
| `pair_observables.py` | observables 3–5 → `pair_relmom.png`, `pair_opening_angle.png`, `missing_mass_pd.png`, `pair_observables_summary.txt` |

Events only exist for u=1 (AV18), u=4 (CD-Bonn), u=6 (pD AV18), u=9 (pD CD-Bonn)
in `events/wf_test/`; the density figure covers all five potentials × {full, pD}.
