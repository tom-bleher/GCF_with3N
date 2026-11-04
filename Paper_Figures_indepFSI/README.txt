PAPER FIGURES — INDEPENDENT FSI, proton-based (e,e'pp) selection
==================================================================
12C, 6 GeV, GENIE hN FSI (independent transport), kF = 220 MeV/c.
Fig5/6/7 regenerated 2026-10-22 with make_paper_figures.sh (Fig5 = the
_triangles variant of plot_fig5_theta_heatmap_pair.py) on the 2026-08-29
production of run_paper_production.sh (read from /home/ilayws/Ilay-Generators):
  3N+FSI : genQE_3N/events/events_3N_100M.root  (100M, -A 12 -Z 6 -f hN -K 0.22
           -s 0.15 -I -q 1:1000000, CM on; events tree built by
           make_events_from_genT.py --ebeam 6.0 --kF 0.22)
  2N+FSI : genQE_FSI/events/misc/2N_hN_100M_6GeV_kF220_indep.root  (100M, AV18,
           sigCM 0.15, fsi_indep=1)
  MF     : genQE_FSI/events/misc/MF_localFG_hN_100M_6GeV_kF220.root  (100M,
           local Fermi gas, wf_mode 4)
The stored lead/recoil/nProtonAboveKF of all three trees were checked to match
the current make_events_from_genT.build_events rule (first 300k entries each).
Fig1/2/4 are unchanged (no-FSI configuration plots, not affected by the
selection). Fig3 not regenerated (needs a C++ histogram refresh).

SELECTION (identical for 3N, 2N and MF; the paper's "two high-momentum
protons are detected and a recoil high-momentum nucleon is reconstructed"):
  detected = final-state PROTON with |p| > kF; exactly 2 (nProtonAboveKF == 2);
             neutrons are neither detected nor vetoed
  lead     = detected proton closest in angle to q; recoil = the other one
  p3       = -(p_lead - q + p_recoil), |p3| > kF
  theta_pq < 8 deg, p_lead/|q| > 0.75, 0.25 < pmiss < 0.9 GeV/c, xB < 1.2,
  Q2 > 1 GeV^2, theta_e < 45 deg (at 6 GeV, Q2 > 1 already implies
  theta_e > 9.6 deg, so the old 8 deg lower cut changes nothing)
  theta12 = angle(p_lead - q, p_recoil), theta23 = angle(p_recoil, p3);
  L/R = triangles with A = 135 deg, K = 4.

REGION FRACTIONS (weighted, of selected events). Uncertainties are MC
statistical: the larger of the sum-of-w^2 estimate (printed by the scripts)
and a 20-batch jackknife over the 150-worker productions.
  3N+FSI : f_L = 14.48 +- 0.31 %  f_R =  8.21 +- 0.11 %  L/R = 1.762 +- 0.044
  2N+FSI : f_L =  7.11 +- 0.17 %  f_R =  3.21 +- 0.08 %  L/R = 2.21  +- 0.09
  MF     : f_L =  5.23 +- 0.30 %  f_R = 12.35 +- 0.47 %  L/R = 0.424 +- 0.030
  Double ratio (L/R)_3N / (L/R)_2N = 0.80 +- 0.04
  (L/R)_2N / (L/R)_MF              = 5.2  +- 0.4
  Fig6: R(x_L) = (1 + x_L)/(1/R_3N + x_L/R_2N), x_L = N_L,2N / N_L,3N:
        R(0) = 1.762, R(0.1) = 1.795 +- 0.042, R(0.5) = 1.890
  Pass fraction of the full selection (weighted, vs all generated):
        3N 4.39e-2, 2N 1.32e-2, MF 9.91e-4  (2N/MF = 13)

With this selection the 2N+FSI background is skewed to the LEFT as much as the
3N signal (L/R 2.2 vs 1.8). 83% of the 2N weight in L is generated
(lead p, partner n) pairs, so the second "detected" proton is an FSI ejectum,
not the SRC partner; such ejecta are soft. LR_recoil_scan.txt shows the L/R of all three
samples vs an extra |p_recoil| threshold (e.g. > 0.40 GeV/c: 3N 0.72,
2N 0.27, MF 0.11).

SUPERSEDED (shipped 2026-10-15, from the OLD nucleon-counting load_3N of
plot_theta_heatmap_3N_vs_2N.py at 87f3d0c: exactly 2 nucleons of either
isospin above kF incl. secondaries, lead = generator N1 of either isospin,
no |p3| cut; 40M 3N / 60M 2N / 40M MF samples of 2026-07-15, since deleted):
  3N+FSI : f_L=2.96%  f_R=6.57%   L/R=0.451
  2N+FSI : f_L=0.80%  f_R=7.50%   L/R=0.107
  MF     : f_L=0.59%  f_R=22.2%   L/R=0.0267  (proton-based, no |p3| cut)
  Double ratio 4.2, R(x_L=0.1) = 0.349.
