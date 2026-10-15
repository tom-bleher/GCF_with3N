PAPER FIGURES — jacobian-fixed generator, INDEPENDENT FSI, DETECTION-BASED 2N
==============================================================================
12C, 6 GeV, GENIE hN FSI, kF=220 MeV/c. SRC cuts (Fig5/Fig7):
theta_pq<8 deg, p_lead/|q|>0.75, 0.25<pmiss<0.9 GeV/c, xB<1.2.
Base: 8<theta_e<45 deg, Q2>1 GeV^2, eq2 tag (exactly 2 measured k>kF).
Fig1/2/4 use the NO-CM sample; Fig5/Fig7 use CM.

DETECTION-BASED 2N (fix of 2026-07-20): the two nucleons whose opening angles are
formed are the lead (fixed by detector cuts, consistent with q) and the UNIQUE
above-kF nucleon = recoil_post if it survived above kF, else the leading above-kF
secondary. (Previously the below-kF recoil_post was used even when a secondary was
the actually-detected nucleon — 18.8% of events — inflating region L.)

Fig1  no-FSI 3N configuration map (theta12-theta23), symmetry axes. [add cartoons]
Fig2  ground-state momenta by config.
Fig3  x_B vs P_miss — NOT regenerated (needs C++ histogram refresh).
Fig4  x_B by config.
Fig5  3N+FSI vs 2N+FSI (independent, detection-based). L/R: 3N 1.185, 2N 0.107.
Fig6  R=N_L/N_R vs x_L, with R_3N=1.185, R_2N=0.107.
Fig7  MF fake-3N (local FG). L/R 0.056.

REGION FRACTIONS (of cut-passing events), independent FSI, paper cuts:
  3N+FSI :  f_L=8.10%  f_R=6.83%   L/R=1.185
  2N+FSI :  f_L=0.80%  f_R=7.50%   L/R=0.107   (detection-based)
  MF     :  f_L=1.26%  f_R=22.3%   L/R=0.056
Double ratio (L/R)_3N / (L/R)_2N = 1.185/0.107 ~ 11.

2N L/R by treatment (paper cuts): with old veto 0.045 | below-kF-recoil (buggy)
0.327 | detection-based (correct) 0.107.
