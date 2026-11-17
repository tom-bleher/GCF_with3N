#ifndef __QE_GENERATOR_FSI_H__
#define __QE_GENERATOR_FSI_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/eNCrossSection.hh"
#include "fsi/GenieFSIHelpers.hh"
#include <vector>

// QEGeneratorFSI — QE generator with Final State Interactions (FSI)
//
// Extends gcfGenerator to transport the outgoing SRC nucleon pair (or a
// single mean-field nucleon) through GENIE's intranuclear cascade (hN or hA
// mode).
//
// hN (HNIntranuke2018): Full intranuclear cascade with hadron-nucleon
//   cross sections and medium corrections. More accurate.
// hA (HAIntranuke2018): Effective cascade using hadron-nucleus cross sections.
//   Faster (~5-10x), suitable for systematics studies.
//
// The FSI output is purely a particle list (GetFinalStateParticles()):
// GENIE's stable final state verbatim, with no lead/recoil designation and
// no absorption special-casing. The generated (pre-FSI) nucleons are
// returned through the generate_event out-parameters unchanged. Which
// final-state particle is "lead" or "recoil" is a detection-level choice
// made downstream from types and momenta alone. Pauli blocking is handled
// inside GENIE's cascade; no additional momentum cut is applied here.

class QEGeneratorFSI: public gcfGenerator
{
 public:
  QEGeneratorFSI(double E, gcfNucleus * thisInfo, eNCrossSection * thisCS,
                 TRandom3 * thisRand);
  ~QEGeneratorFSI();

  // Event generators (same interface as QEGenerator). The out-parameters
  // hold the GENERATED (pre-FSI) nucleons; the post-FSI final state is in
  // GetFinalStateParticles().
  void generate_event(double &weight, int &lead_type, int &rec_type,
                      TLorentzVector &vk_target,
                      TLorentzVector &vLead_target,
                      TLorentzVector &vRec_target,
                      TLorentzVector &vAm2_target);
  void generate_event(double &weight, int &lead_type, int &rec_type,
                      TLorentzVector &vk_target,
                      TLorentzVector &vLead_target,
                      TLorentzVector &vRec_target,
                      TLorentzVector &vAm2_target,
                      double &Estar);
  void generate_event(double &weight, int &lead_type, int &rec_type,
                      TLorentzVector &vk_target,
                      TLorentzVector &vLead_target,
                      TLorentzVector &vRec_target,
                      TLorentzVector &vAm2_target,
                      TVector3 &vRel_target,
                      TLorentzVector &q_target,
                      double &Estar);
  void generate_event_lightcone(double &weight, int &lead_type, int &rec_type,
                                 TLorentzVector &vk_target,
                                 TLorentzVector &vLead_target,
                                 TLorentzVector &vRec_target,
                                 TLorentzVector &vAm2_target);
  void generate_event_lightcone(double &weight, int &lead_type, int &rec_type,
                                 TLorentzVector &vk_target,
                                 TLorentzVector &vLead_target,
                                 TLorentzVector &vRec_target,
                                 TLorentzVector &vAm2_target,
                                 double &Estar);

  // ---- Single-nucleon mean-field (MF) mode ----
  // Models a single struck nucleon drawn from a Fermi-gas momentum
  // distribution (no correlated pair), transported through GENIE FSI in the
  // A-1 residual. No gen-level recoil is emitted: rec_type is set to a
  // non-nucleon sentinel and vRec is zero, so any second proton in an
  // (e,e'pp) analysis must originate from the FSI particle list — the
  // correct MF picture. FSI is still toggled by EnableFSI().
  enum FermiGasMode { kGlobalFG, kLocalFG };
  void generate_event_MF(double &weight, int &lead_type, int &rec_type,
                         TLorentzVector &vk_target,
                         TLorentzVector &vLead_target,
                         TLorentzVector &vRec_target,
                         TLorentzVector &vAm2_target);
  // kGlobalFG: flat n(p)=3/(4 pi kF^3) up to kF.
  // kLocalFG : density-folded n(p) from a local Fermi gas, kF(r) set by the
  //            one-body density (built lazily, cached per nucleus).
  void SetFermiGasMode(FermiGasMode m) { fFGMode = m; }
  // Fermi momentum ceiling for the GLOBAL Fermi gas (GeV/c). For the local
  // Fermi gas the ceiling is kF at central density (computed internally).
  void SetFermiMomentumFG(double kF_GeV) { fkF = kF_GeV; fLocalFG_built = false; }
  // Density value n_FG(p) used as the MF "wavefunction" (for bookkeeping in
  // the analysis tree). Normalised so integral n_FG d^3p = 1.
  double get_n_FG(double p) { return eval_n_FG(p); }

  // ---- SRC onset: where the correlated-pair branch begins (2N pair mode) ----
  //
  // kFixedKF (default): k_rel > a constant k_F, namely the sampling floor the
  //   caller configured through set_pRel_range(). Simple and defensible -- the
  //   GCF is a short-distance factorization -- but it leaves a hole in the
  //   momentum distribution: the mean-field sample is capped at the LOCAL
  //   k_F(r), so nucleons with k_F(r) < k < k_F belong to neither sample.
  //
  // kLocalKF: k_rel > k_F(r), the local Fermi momentum at the sampled pair
  //   vertex. This is N. Wright's convention (GCF_GENIE, LocalFGM::ProbDistro:
  //   flat Fermi gas below k_F(r), AV18 universal function above), and it makes
  //   the SRC and mean-field samples exactly complementary -- together they
  //   tile the whole momentum distribution with no gap and no overlap.
  //   Caveat: k_F(r) runs down to ~0.08 GeV/c at large radius, well below the
  //   ~0.3 GeV/c where the GCF is validated against ab-initio calculations.
  enum SrcOnset { kFixedKF, kLocalKF };
  void SetSrcOnset(SrcOnset m) { fSrcOnset = m; }
  SrcOnset GetSrcOnset() const { return fSrcOnset; }

  // Additionally require |p1| = |k_rel + P_cm/2| > k_F(r) AFTER the c.m.
  // motion is added, i.e. Wright's second test (FermiMover.cxx:182,
  // "FermiMomentum() < Momentum()"). Off by default: her code decides the
  // SRC/mean-field branch twice -- once on k_rel when choosing whether to add
  // P_cm/2, then again on |p1| when choosing whether to emit a recoil -- and
  // the two disagree near the boundary, which is where her one-way leak into
  // the mean-field class comes from. The GCF split belongs on the RELATIVE
  // momentum. Provided so the size of that leak can be measured.
  // Events failing the test get weight = 0 (they belong in the MF sample);
  // they still count as generation attempts, so attempts-based luminosity
  // normalization stays correct.
  void SetSrcRequireP1(bool b = true) { fSrcRequireP1 = b; }
  bool SrcRequireP1() const { return fSrcRequireP1; }

  // k_F(r) at the pair vertex of the last generated 2N event (GeV/c).
  double GetLastLocalKF() const { return fLastLocalKF; }

  // FSI on/off toggle (enabled by default)
  void EnableFSI(bool enable = true) { doFSI = enable; }

  // Select FSI model (default: kHN2018)
  void SetFSIModel(FSIModel m) { fFSIModel = m; }

  // Independent-FSI mode: transport each nucleon through its own GENIE
  // record with an undepleted A-2 remnant, instead of the default shared
  // depleting remnant in a single record.
  void SetIndependentFSI(bool i = true) { fIndependentFSI = i; }
  bool IndependentFSI() const { return fIndependentFSI; }

  // Post-FSI final state of the last generated event: a flat particle list
  // (parentRole = 0/1 for descendants of the generated lead/recoil,
  // diagnostic only). With FSI disabled (or skipped for an unphysical
  // residual) this is simply the generated nucleon(s).
  const std::vector<FSIParticle> & GetFinalStateParticles() const {
    return fParticles;
  }

  // Negative particles dropped by the -C correction of the final-state list
  // (CoulombCorrectFinalState), summed over all events so far. Only events with
  // weight > 0 carry a particle list, so this counts stored events only.
  long GetNCoulombDropped() const { return fNCoulombDropped; }

  // Fermi-gas density n_FG(|p1|) of the last generate_event_MF call,
  // evaluated at the TRUE sampled struck-nucleon momentum (free of the
  // radiative/Coulomb reconstruction ambiguity). Valid for MF events only.
  double GetLastMFRho() const { return fLastMFRho; }

  // Interaction vertex of the last generated event -- the one the cascade
  // uses, or would use with FSI on. Drawn during generation (it fixes k_F(r)
  // for the local onset / local Fermi gas), so a PWIA-only run must hand it
  // to a separate FSI stage rather than let that stage re-sample it.
  //   2N pair events: the SRC pair vertex (false if none was drawn, e.g. the
  //                   light-cone path);
  //   MF events:      the struck-nucleon position.
  bool GetLastSrcVertex(TLorentzVector &x4) const {
    if (fSrcVertexValid) x4 = fSrcVertex;
    return fSrcVertexValid;
  }
  const TLorentzVector & GetLastMFPosition() const { return fLastMFPosition; }

 private:
  eNCrossSection * myCS;

  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;

  // FSI state
  bool   doFSI;
  FSIModel fFSIModel;
  bool   fIndependentFSI = false;
  int    fA, fZ;           // nucleus A and Z (set from gcfNucleus in ctor)
  std::vector<FSIParticle> fParticles;
  long fNCoulombDropped = 0;

  // SRC onset state (2N pair mode)
  SrcOnset fSrcOnset     = kFixedKF;
  bool     fSrcRequireP1 = false;
  double   fLastLocalKF  = 0.;   // k_F(r) at the last pair vertex (GeV/c)
  // Pair vertex of the current event. Sampled at the TOP of generate_event so
  // k_F(r) is available before the relative momentum is drawn; the same
  // vertex is then handed to the cascade, so the FSI entry point is unchanged
  // -- only the order in which it is drawn.
  TLorentzVector fSrcVertex;
  bool fSrcVertexValid = false;

  // ---- kLocalKF vertex sampling ----
  // With a fixed onset the SRC strength a vertex carries, I_SRC =
  // int_{kF}^{pRelmax} k^2 S(k) dk, is the same everywhere, so drawing the
  // vertex from rho^n(r) r^2 and leaving I_SRC in the weight is fine.
  // With the LOCAL onset I_SRC becomes I_SRC(r) and grows by orders of
  // magnitude as kF(r) -> 0 at large radius, so that recipe puts almost all
  // the weight on a handful of rare far-out vertices (measured: 3% of events
  // carried 48% of the weight, N_eff/N = 0.25%).
  //
  // Instead draw the vertex from rho^n(r) r^2 I_SRC(r) directly. That is both
  // the variance fix and the physically correct pair-vertex distribution --
  // with a radius-dependent onset, SRC pairs are no longer distributed like
  // the one-body density. The leftover weight factor is then the constant
  // <I_SRC>_{rho^n r^2}, independent of r.
  void build_local_onset_vertex_table(int lead_type, int rec_type);
  TLorentzVector sample_local_onset_vertex(int lead_type, int rec_type,
                                           double &kF_r, double &I_src_r);
  bool   fLocalVtx_built[3] = {false, false, false};
  double fLocalVtx_rmax[3]  = {0., 0., 0.};   // fm
  double fLocalVtx_fmax[3]  = {0., 0., 0.};   // max of rho^n r^2 I_SRC(r)
  double fLocalVtx_meanI[3] = {0., 0., 0.};   // <I_SRC> over rho^n r^2

  // Single-nucleon MF state
  static constexpr int kNoRecoil = 0;  // rec_type sentinel: no recoil nucleon
  FermiGasMode fFGMode = kGlobalFG;
  double fkF = 0.25;                    // global-FG Fermi momentum (GeV/c)
  double fLastMFRho = 0.;               // density at last MF struck-nucleon
                                        // (r,|p|): n_local(p|r) for LFG,
                                        // 3/(4 pi kF^3) for global FG.
  TLorentzVector fLastMFPosition;       // sampled (r) used as FSI entry point
  bool   fLocalFG_built = false;        // local-FG marginal table cached for fA?
  double fLocalFG_kmax = 0.;            // local-FG momentum ceiling = kF(r=0)
  std::vector<double> fLocalFG_p;       // local-FG marginal table: p grid (GeV/c)
  std::vector<double> fLocalFG_n;       // local-FG marginal table: n_MF(p)

  // Transport the generated nucleon(s) through GENIE and fill fParticles
  // with the resulting flat final-state particle list. Does not modify the
  // generated types/momenta or the event weight.
  void ApplyFSI(int lead_type, int rec_type,
                const TLorentzVector &vLead_target,
                const TLorentzVector &vRec_target);

  // Single-nucleon variant (mean-field knockout, A-1 residual).
  void ApplyFSI_single(int lead_type, const TLorentzVector &vLead_target);

  // Seed fParticles with the generated nucleon(s) unchanged (the final state
  // when FSI is off or skipped).
  void SeedFinalState(int lead_type, int rec_type,
                      const TLorentzVector &vLead_target,
                      const TLorentzVector &vRec_target);

  // Coulomb correction (-C) of the final-state particle list, the same shift
  // the generated nucleons get: each charged particle leaving the nucleus
  // gains q*deltaE (q in units of e). A negative particle with less kinetic
  // energy than that cannot escape and is dropped.
  void CoulombCorrectFinalState(double deltaE);

  // Shared electron/Jacobian/cross-section kinematics (extracted so both the
  // 2N pair path and the single-nucleon MF path use identical physics).
  // Given the struck-nucleon 3-momentum v1, its off-shell 4-vector v1_target
  // (with energy E1) and type, solves the scattered-electron kinematics and
  // fills vk_target/vLead_target/q_target and the weight. Returns false (and
  // sets weight=0) on any kinematic rejection.
  bool solve_lepton_kinematics(double &weight, int lead_type,
                               const TVector3 &v1,
                               const TLorentzVector &v1_target, double E1,
                               TLorentzVector &vk_target,
                               TLorentzVector &vLead_target,
                               TLorentzVector &q_target);

  // Sample the struck-nucleon position r AND momentum v1 from the Fermi-gas.
  // In local-FG mode the position is sampled from the one-body density and
  // |p| is then drawn from the local Fermi sphere at that r (per-isospin
  // k_F^q(r) — N_q = Z for protons, A-Z for neutrons); r and v1 carry the
  // physically correct (r,p) correlation. In global-FG mode r is sampled
  // independently of p and |p| is uniform in [0, fkF].
  bool sample_FG_momentum(int nucleon_type, double &weight, TVector3 &v1);

  // Fermi-gas density n_FG(p) (normalised so integral n d^3p = 1).
  double eval_n_FG(double p) const;
  // Build (and cache) the local-Fermi-gas momentum table for nucleus fA.
  void build_local_FG_table();
};

#endif
