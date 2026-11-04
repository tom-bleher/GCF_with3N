# GCF_with3N: what the generator computes

A short account of the physics in the generators and the FSI stage, with the
expressions as they are coded. File and line references are to the code at
the time of writing; the code wins where they disagree. For how to run
things, see [PIPELINE.md](PIPELINE.md).

## Contents
- [The picture](#the-picture)
- [Electron–nucleon cross section](#electronnucleon-cross-section)
- [2N: struck nucleon from an SRC pair](#2n-struck-nucleon-from-an-src-pair)
- [3N: struck nucleon from an SRC triplet](#3n-struck-nucleon-from-an-src-triplet)
- [MF: single mean-field nucleon](#mf-single-mean-field-nucleon)
- [Event weights and normalization](#event-weights-and-normalization)
- [Final-state interactions (GENIE)](#final-state-interactions-genie)
- [Assumptions and known limitations](#assumptions-and-known-limitations)

## The picture

Quasi-elastic electron scattering on a nucleus, A(e,e′N…), in the **plane-wave
impulse approximation (PWIA)** with one-photon exchange: the virtual photon is
absorbed by one nucleon, which leaves the nucleus; the rest of the nucleus is a
spectator. The cross section factorizes into

    dσ = σ_eN(off-shell) × [probability of finding the struck nucleon with
                            momentum p and energy E in the initial state]

The second factor is where the **Generalized Contact Formalism (GCF)** enters.
Its assumption: at short distance (high relative momentum, above a cutoff of
order the Fermi momentum) the nucleus looks like a collection of correlated
pairs (or triplets) whose relative wave function is universal — the same
two-body function `φ_NN(k)` for every nucleus — multiplied by a
nucleus-dependent *contact* `C_NN` that counts how many such pairs there are.
The pair as a whole carries a centre-of-mass momentum drawn from a Gaussian.

The generators sample the initial-state configuration (pair or triplet
momenta, CM motion), fix the electron kinematics by energy–momentum
conservation, and weight the event by the cross section and the initial-state
density. **Final-state interactions** are then applied, as a separate step or
inline, by transporting the outgoing nucleons through GENIE's intranuclear
cascade in the residual nucleus.

## Electron–nucleon cross section

`src/libraries/cross_sections/eNCrossSection.cc`, following
T. de Forest, Nucl. Phys. A 392 (1983) 232
([doi:10.1016/0375-9474(83)90124-0](https://doi.org/10.1016/0375-9474(83)90124-0)).

Off-shell struck nucleon with momentum **p** and energy Ē = √(p² + m²) ("free"
energy), on-shell outgoing nucleon **p′**, physical momentum transfer
q = (ω, **q**). The cross section is de Forest's eq. (13):

    σ_eN = σ_Mott [ (Q²/|q|²)² W_C + (Q²/(2|q|²) + tan²(θ/2)) W_T
                    + (Q²/|q|²) (Q²/|q|² + tan²(θ/2))^½ W_I cos φ
                    + (Q²/|q|² cos²φ + tan²(θ/2)) W_S ]

with σ_Mott = 4α² E′² cos²(θ/2) / Q⁴ and φ the angle between the scattering
plane and the (**q**, **p′**) plane. The responses W_C, W_T, W_S, W_I come
from the nucleon current evaluated with free spinors and free form factors
(Kelly parametrization, PRC 70, 068202,
[doi:10.1103/PhysRevC.70.068202](https://doi.org/10.1103/PhysRevC.70.068202);
G_E and G_M for proton and neutron), in one of two
off-shell prescriptions:

- **cc1** (`-c cc1`, the default): de Forest's eq. (17), the simpler current
  with q̄ = p′ − p̄.
- **cc2** (`-c cc2`): eq. (18), the current with σ_μν q_ν.
- **onshell**: the free elastic eN cross section at the same (E, θ), no
  off-shell prescription.

Form factors are evaluated at the physical Q² (de Forest's and SIMC's choice;
GENIE/Wright
[doi:10.48550/arXiv.2104.05090](https://doi.org/10.48550/arXiv.2104.05090)
evaluate them at Q̄²). For a free nucleon at rest all three
options coincide; this is the check used in #32, which found and fixed three
errors in the coded cc1/cc2 expressions (metric sign, a parenthesis in W_C,
the longitudinal factor). Until that fix is merged, the coded cc2 is wrong by
orders of magnitude and cc1 by 5–30 %.

## 2N: struck nucleon from an SRC pair

`src/libraries/generator/gcfGenerator.cc`, `src/programs/genQE_FSI/QEGeneratorFSI.cc`,
nuclear input in `src/libraries/nucleus/gcfSRC.cc`, `gcfNucleus.cc`.

**Initial state.** A pair of nucleons (pp, pn or nn) with

- centre-of-mass momentum **P**_cm, each component Gaussian with width σ_CM
  (nucleus-dependent, e.g. 0.15 GeV/c for ¹²C; `-s` overrides);
- relative momentum **k** = (**p**₁ − **p**₂)/2 with magnitude drawn from
  k² S_NN(k) between k_cut (default 0.25 GeV/c, `-k`; or the local Fermi
  momentum k_F(r) at the pair's position with `onset = local`) and 1 GeV/c,
  and isotropic direction;
- S_NN(k) = Σ C_NN^(l) |φ_NN^(l)(k)|², the contact-weighted universal
  functions tabulated per NN interaction (AV18 by default; N2LO, AV4 and NV
  available; N3LO has tables but no contacts, so it aborts). The contacts are the "Scale and Scheme" k-space values
  (Cruz-Torres et al., Nature Phys. 17, 306 (2021),
  [doi:10.1038/s41567-020-01053-7](https://doi.org/10.1038/s41567-020-01053-7)).

The struck nucleon has **p**₁ = **P**_cm/2 + **k**, its partner **p**₂ =
**P**_cm/2 − **k**; the A−2 residual recoils with −**P**_cm and mass
m_{A−2} + E*, E* an optional excitation energy (`-E`, `-M`). The struck
nucleon is off-shell: E₁ = m_A − E_{A−2} − E₂.

**Lepton kinematics.** Q² (within a window) and the electron azimuth φ_k are
sampled; the scattered-electron energy is the root of energy conservation for
an on-shell outgoing nucleon **p**′ = **p**₁ + **q**. The delta function is
removed with the Jacobian J (coded in `solve_lepton_kinematics`; checked
analytically in the audit).

**Weight** (per generation attempt):

    w = [pair-type factor] × σ_CM-Gaussian (implicit in the sampling)
        × Δφ_rel Δcosθ_rel Δk  ·  k² S_NN(k) / (2π)³          (relative motion)
        × ΔQ² Δφ_k  /  J                                     (lepton phase space)
        × σ_eN(E, k, p₁)                                    (cross section)
        × [radiation factor]                                (with -O)

With `--prel-is` the relative momentum is drawn from k² S(k) itself and the
weight carries the integral over the sampled range instead (same expectation
values, far lower variance: N_eff/N goes from ~2.5 % to order 1).

**Options.** Peaking-approximation radiative correction (`-O`); Coulomb
correction of the electron and proton energies by the nucleus' Coulomb
potential (`-C`); randomized nuclear parameters for systematics (`-r`);
light-cone variables (`-l`); a deuteron mode (A = 2).

## 3N: struck nucleon from an SRC triplet

`src/programs/genQE_3N/QEGenerator_3N.cc`, density tables in `genQE_3N/array/`.

**Initial state.** Three nucleons (ppn or pnn) whose internal momenta follow
the three-body density ρ(p_a, p_b, θ_ab) of **³He**, computed from Faddeev
wave functions for several interactions (AV18, AV18+UIX, AV18+TM, CD-Bonn,
CD-Bonn+TM; `-pD` variants keep only the deuteron-channel part). Here
**p**_a is the momentum of the odd nucleon relative to the pair and **p**_b
the pair's relative momentum, θ_ab the angle between them. The same ³He
density is used for the triplet in every target nucleus — the 3N analogue of
the universal pair function — and the triplet as a whole gets a Gaussian CM
momentum (σ_CM = 0.15 GeV/c by default, or √(3/5) k_F with `-K`), shared
equally among the three nucleons. Internal momenta below k_F (default 0.25
GeV/c) are rejected: that is the SRC condition. The A−3 residual recoils with
−**P**_cm.

The tables are normalized as a Cartesian density per d³p_a d³p_b:
8π² ∫ ρ p_a² p_b² sinθ_ab dp_a dp_b dθ_ab = 1/3 (the isospin assignment; the
three choices of which nucleon is the neutron sum to 1).

**Sampling and weight.** p_a, p_b (0–5 GeV/c), θ_ab and the orientation
angles are sampled uniformly; the electron direction within an angular window
(`-t`), its energy fixed by energy conservation with Jacobian
J_δ = E′₁ / |E′₁ − **p**′₁·**k̂**|. The weight is

    w = Δ(angles) Δp_a Δp_b Δθ_ab  ·  p_a² p_b² sinθ_ab ρ(p_a, p_b, θ_ab)
        × J_δ × σ_eN × ½ (2π)⁻⁷ C t

where the last factor, with C = 0.12 and t = 2, is a placeholder: the
absolute normalization of the 3N channel is **not** defined relative to the 2N
one, so 2N/3N comparisons must be normalized by other means. The uniform
sampling over a density spanning many decades makes the weights very uneven
(N_eff/N ≈ 0.4 %); an importance sampler like the 2N one would help.

## MF: single mean-field nucleon

`QEGeneratorFSI::generate_event_MF`. One nucleon is drawn from a Fermi gas:
either a global one (|p| uniform in a sphere of radius k_F = 0.25 GeV/c) or a
local one (`--fermi-gas local`: position from the one-body density, then |p|
from the local Fermi sphere with k_F(r) ∝ ρ(r)^{1/3}, as in GENIE's
LocalFGM). The A−1 residual recoils with −**p**; there is no correlated
partner. This sample is meant to complement the SRC sample: with the local
onset the two tile the momentum distribution without gap or overlap.

Known: the MF weight is normalized per nucleon and divided by (2π)³ a second
time, so it is not on the same absolute scale as the SRC weight (#34); the
samples can be compared only after normalizing each.

## Event weights and normalization

Every generator keeps producing until the requested number of events with
weight > 0 has been written, and counts every attempt in `hAttempts`. The
quantity to compare between samples is the weighted yield per attempt,
Σw / attempts (printed by `gcf_generate`, shown by `gcf quicklook`). Events
rejected at generation time (zero weight: kinematic cuts, k < k_F, Q² window)
are part of the attempts, so cuts applied at generation and cuts applied later
give the same normalization. Weights vary over orders of magnitude, so the
effective sample size N_eff = (Σw)²/Σw² matters more than the event count.

## Final-state interactions (GENIE)

`src/libraries/fsi/GenieFSIHelpers.cc`, `FSIStage.cc`; inline in
`QEGeneratorFSI::ApplyFSI` (2N, MF) and `QEGeneratorFSI_3N::ApplyFSI` (3N);
as a stage in `gcf_fsi`.

The outgoing nucleons are transported through **GENIE's INTRANUKE 2018**
cascade (vendored GENIE R-3_06_02):

- **Model:** hN (full hadron–nucleon cascade: each step is an NN or πN
  collision with free cross sections and medium corrections) or hA (an
  effective hadron–nucleus model with parametrized fates: charge exchange,
  inelastic, absorption, pion production). In hN a nucleon's fates are
  elastic, inelastic (pion production) and compound-nucleus formation; what hA
  calls nucleon absorption or charge exchange arises in hN from sequences of
  collisions, not as a separate fate.
- **Medium:** the residual nucleus, A minus the knocked-out nucleons and Z
  minus the knocked-out protons (A−2, A−3 or A−1). The cascade is skipped if
  the residual has A < 3. GENIE's cascade nucleons are drawn from a local
  Fermi gas (project override `config/genie_override/*.xml`:
  `NuclearModel = LocalFGM`), consistent with the generator's local Fermi gas.
- **Starting point:** all nucleons of the cluster start from one point,
  drawn from ρ(r)ⁿ r² with ρ the one-body density of the *full* target and
  n = 1 (single nucleon), 2 (pair; configurable, Wright
  [doi:10.48550/arXiv.2104.05090](https://doi.org/10.48550/arXiv.2104.05090)
  uses 1) or 3
  (triplet): the cluster formed before the knockout, and its probability
  scales with the density to the power of its size. In the 2N generator the
  vertex is drawn during generation because, with the local onset, k_F at the
  vertex enters the sampling; a PWIA file passes it on to the FSI stage.
- **Record:** a GENIE event record in lepton–nucleus mode with a fake
  electron, the target, the remnant (PDG mass) and the cluster nucleons as
  "hadrons in the nucleus"; `TransportHadrons` steps each one through the
  medium. In **shared** mode all nucleons go into one record in random
  order, so a nucleon absorbed during one cascade depletes the remnant seen by
  the next; in **independent** mode each gets its own record and the full
  remnant.
- **Output:** GENIE's stable final state verbatim as a flat particle list
  (`part_*`), with the fate of the input nucleon each particle descends from
  (`rescatterCode`; 1 = left without interacting), the mother-chain depth
  (`generation`, ≥ 1 for every cascade output) and the role of its ancestor
  (`parentRole`). No particle is labelled "lead" or "recoil": that is a
  detection-level choice made in the analysis. No Pauli veto is applied
  after the cascade; in this GENIE version the NN-elastic Pauli block inside
  the cascade is off (#10).
- **Channels:** the rate of each fate can be scaled or switched off
  (`gcf_fsi --channel`, `[fsi] off`), through GENIE's own per-channel
  factors, for studies of which interaction does what.
- **Known issues fixed on branches:** an unbounded retry inside GENIE that
  overflowed the stack for nucleons left with no allowed fate (#1), and the
  Coulomb correction not being applied to the post-FSI list (#35).

## Assumptions and known limitations

- PWIA and factorization: no interference between nucleons, no
  meson-exchange currents, no two-body currents in the primary vertex.
- The universal ³He three-body density stands in for the triplet in every
  nucleus; the absolute 3N normalization is a placeholder (see above).
- The pair CM motion is Gaussian and isotropic; the SRC region starts at a
  sharp cutoff (fixed k_cut or local k_F).
- The off-shell nucleon current is a prescription (cc1/cc2); form factors at
  the physical Q².
- GENIE INTRANUKE is validated for A > 6; light residuals (A_res < 7) get a
  warning. The cascade medium is GENIE's generic density for the residual's
  A, while the vertex is drawn from the full target's density.
- Open physics questions are tracked as issues on the fork: #9 (a Jacobian in
  the 3N plotting tools), #10 (Pauli blocking), #34 (MF normalization), #35
  (Coulomb with FSI), and the conventions noted in #32.
