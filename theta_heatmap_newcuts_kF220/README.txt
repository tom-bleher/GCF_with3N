theta12-theta23 heatmaps — JACOBIAN-FIX RERUN (2026-07-14)
============================================================
12C, 6 GeV, GENIE hN FSI, kF = 220 MeV/c, sigmaCM = sqrt(3/5)kF (3N) /
sqrt(2/5)kF (2N) / 0 (MF). Selection: exactly 2 measured nucleons with
k > kF (eq2); p3 reconstructed as -(pmiss + p_recoil). 150 x 150 bins,
color floor 1e-6 of each panel max. Independent color scale per panel.

BASE CUTS (everywhere): 8 < theta_e < 45 deg, Q2 > 1 GeV^2, eq2 tag
(+ theta_pq < 12 in the nocuts plots). NO xB CUT IN BASE. xB windows
appear ONLY as explicitly chosen cuts: 0.6<xB<1.2 in the newcuts sets,
0.8<xB<1.2 in frontiercuts, and xB<1.2 inside the OLD SRC cut set.
MF heatmaps and no-FSI region plots use the xB-free base (rebuilt).
Numbers in the *_summary.txt files.

CODE CONFIGURATION (all samples generated 2026-07-14 with it):
- 3N generator jacobian fix (QEGenerator_3N, user edit of Jul 12/13).
- No post-cascade Pauli/kF veto anywhere (GENIE-internal blocking only).
- Fate-code absorption (validated) with ejecta promotion: an absorbed leg's
  saved momentum is its leading nucleon ejecta (experiment measures the
  final state); the absorbed flag is bookkeeping only.
- FSI transport: shared depleting remnant (single GENIE record), with
  RANDOMIZED transport order (second fix of 2026-07-14): GENIE transports
  the record in insertion order while sharing the depleting remnant, so a
  fixed lead-first order gave the lead too much FSI and the last recoil too
  little (validated: pre-fix rescatter fractions lead 0.5125, N2 0.4327,
  N3 0.4234; post-fix 0.5006 / 0.4304 / 0.4308 — N2/N3 symmetric). The
  order is now shuffled per event and results demuxed by parentRole.
  MF (single nucleon) and PWIA (no FSI) are unaffected by this fix.

SAMPLES (all regenerated 2026-07-14 with the fixed code):
- 3N shared FSI: 40M, generated with -q 1:1e9 (Q2 > 1 GeV^2 enforced at
  generation, before FSI). Valid for all plots here (they all cut Q2 > 1);
  do NOT reuse for Q2 < 1 studies. Pass fractions in summaries are
  relative to this restricted sample.
- 2N shared FSI: 60M (tree includes generation-failure placeholders with
  weight = 0; analyses mask weight > 0).
- MF global / local FG: 40M each.
- 3N PWIA (no FSI): 20M; 3N PWIA no-CM: 20M.
NOTE: the Jul-13 files with the same names were a different (RNTuple)
production; they are preserved as *_rntuple_jul13.root and NOT used here.

RESULTS (final, order-randomized samples; *_highres.pdf variants use
300 x 300 bins with a 1e-8 color floor)
  sharedFSI_pmiss_0.25-0.90 (theta_pq<12, 0.6<p/q<0.95, 0.6<xB<1.2,
    0.25<pmiss<0.90):        (L/R)_3N = 1.156   (L/R)_2N = 0.244
  standard_plots/:
    nocuts (theta_e 8-45, Q2>1, eq2, theta_pq<12; NO xB cut):
                             (L/R)_3N = 1.151   (L/R)_2N = 0.218
    OLD-SRC-cuts (thpq<8, p/q>0.75, 0.25<pmiss<0.9, xB<1.2):
                             (L/R)_3N = 1.154   (L/R)_2N = 0.346
    MF fake-3N (pp): global L/R base 0.055 / srccuts 0.041
                     local  L/R base 0.050 / srccuts 0.056
    3N no-FSI TRUE-angle region plots:
      with CM:  L 4.70% R 4.80% BR 4.23% star 11.03%  (L/R = 0.98)
      no CM:    L 8.13% R 8.33% BR 7.04% star 22.02%  (L/R = 0.98)

PHYSICS NOTE: the jacobian fix removed the strong 3N left-skew seen in all
earlier rounds. The no-FSI opening-angle distribution WITH CM smearing is
now symmetric (L/R = 0.98; it was 0.46 pre-fix), and the 3N+FSI L/R is
~1.06 at these cuts (vs 1.7-7.5 in the pre-fix rounds). The huge-weight
region-L events of earlier samples came from the 1/(p1*p2*p3)-type
divergence the fix addresses. 2N is essentially unchanged (~0.2-0.3).
Consequence: at the current cuts the 3N/2N region-L contrast is now
f_L 0.087 vs 0.028 (ratio ~3) and L/R 1.06 vs 0.22 -- much weaker than
pre-fix. Cut re-optimization on these samples is recommended if region-L
discrimination is the goal.

INDEPENDENT-FSI SET (2026-07-15): indepFSI_* folders mirror the sharedFSI
plots with identical cuts, but each nucleon is transported in its own GENIE
record through an undepleted A-2/A-3 remnant (-I / fsi_indep=1). The
transport-order randomization is a no-op in this mode (no shared remnant).
Samples: 3N_FSIindep_hN_40M_12C_kF220_6GeV_Q2gt1.root (40M),
events_2N_hNindep_6GeV_kF220_60M.root (60M), same code base otherwise.
  shared vs indep, (L/R)_3N / (L/R)_2N:
    pmiss 0.25-0.90:  shared 1.156/0.244   indep 1.101/0.221
    frontiercuts:     shared 2.524/0.070   indep 2.057/0.097
    pmiss>0.55:       shared 2.649/0.942   indep 3.354/0.643
    nocuts:           shared 1.151/0.218   indep 1.087/0.212
    srccuts_OLD:      shared 1.154/0.346   indep 1.185/0.327

AB_oldcode/: 5M reference with the OLD committed code (per-nucleon FSI,
kill-on-absorption, 220-MeV veto, kF=0.25, sigCM=0.15, pre-fix jacobian):
nocuts (L/R)_3N=9.71 / (L/R)_2N=0.018; old SRC cuts 5.25 / 0.032 —
reproduces the historical ~10 / ~1:100 (which therefore also carried the
jacobian artifact).

Removed as stale (pre-fix physics): indepFSI folders and
sharedFSI_pmiss_gt0.55 — rebuildable on request (needs a 3N indep and/or
2N indep regeneration with the fixed code).
