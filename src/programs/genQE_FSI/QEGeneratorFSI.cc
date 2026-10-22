#include "QEGeneratorFSI.hh"
#include "constants.hh"
#include "helpers.hh"
#include "fsi/GenieFSIHelpers.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <TDatabasePDG.h>
#include <TMath.h>

using namespace std;

QEGeneratorFSI::QEGeneratorFSI(double E, gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{
  myCS = thisCS;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  QSqmin = 1.0;
  QSqmax = 10.0;
  phikmin = 0.;
  phikmax = 2.*M_PI;

  // FSI defaults
  doFSI        = true;
  fFSIModel    = kHN2018;
  fA           = thisInfo->get_A();
  fZ           = thisInfo->get_Z();

  fsi::ResolveGenieXMLPath();
}

QEGeneratorFSI::~QEGeneratorFSI()
{
}

// ---------------------------------------------------------------------------
// SeedFinalState — the final state when FSI is off or skipped: the generated
// nucleon(s) unchanged. ApplyFSI replaces this with GENIE's cascade output
// when the cascade actually runs.
void QEGeneratorFSI::SeedFinalState(int lead_type, int rec_type,
                                    const TLorentzVector &vLead_target,
                                    const TLorentzVector &vRec_target)
{
  fParticles.clear();
  auto seed = [&](int pdg, const TLorentzVector &p4, int role) {
    if (pdg != pCode && pdg != nCode) return;
    FSIParticle part;
    part.pdg = pdg; part.p4 = p4; part.rescatterCode = -1; part.parentRole = role;
    fParticles.push_back(part);
  };
  seed(lead_type, vLead_target, 0);
  seed(rec_type, vRec_target, 1);
}

void QEGeneratorFSI::CoulombCorrectFinalState(double deltaE)
{
  std::vector<FSIParticle> kept;
  kept.reserve(fParticles.size());
  for (FSIParticle &part : fParticles) {
    const TParticlePDG *pd = TDatabasePDG::Instance()->GetParticle(part.pdg);
    const double q = pd ? pd->Charge() / 3. : 0.;   // TParticlePDG: units of |e|/3
    if (q < 0. && part.p4.E() - part.p4.M() <= -q * deltaE) { ++fNCoulombDropped; continue; }
    if (q != 0.) coulombCorrection(part.p4, q * deltaE);
    kept.push_back(part);
  }
  fParticles.swap(kept);
}

// ---------------------------------------------------------------------------
// ApplyFSI — transport the generated pair through GENIE's intranuclear
// cascade (hN or hA mode) and store the resulting flat final-state particle
// list in fParticles. The generated types/momenta and the event weight are
// NOT modified: the cascade output is purely a particle list, and any
// lead/recoil identification happens downstream at detection level.
void QEGeneratorFSI::ApplyFSI(int lead_type, int rec_type,
                               const TLorentzVector &vLead_target,
                               const TLorentzVector &vRec_target)
{
  // Residual nucleus after removing the struck/knocked-out SRC pair.
  const int removedP = (lead_type == pCode ? 1 : 0) + (rec_type == pCode ? 1 : 0);
  const int removedN = (lead_type == nCode ? 1 : 0) + (rec_type == nCode ? 1 : 0);
  const int Z_res = fZ - removedP;
  const int N_res = (fA - fZ) - removedN;
  const int A_res = Z_res + N_res;

  // No meaningful residual medium left to traverse.  Skip FSI for A_res<3:
  // GENIE's intranuke cascade is not validated for A<3 (produces NaN in
  // cross-section splines), and a ≤2-nucleon residual provides negligible
  // rescattering in any case.
  if (A_res < 3 || Z_res < 0 || N_res < 0) {
    return;
  }

  // Light residual nucleus warning (one-time)
  if (A_res < 7) {
    static bool warned = false;
    if (!warned) {
      std::cerr << "Warning: GENIE intranuke applied to light residual nucleus "
                << "(A_res=" << A_res << "). GENIE is validated for A>6; "
                << "results for A_res<7 should be treated with caution." << std::endl;
      warned = true;
    }
  }

#ifdef USE_FSI
  // Position of the SRC pair (both nucleons originate from the same
  // location), weighted by rho^n of the FULL nucleus since the pair was
  // formed before knockout. generate_event() has normally already drawn it
  // (it needs k_F(r) before the relative momentum); sample here only for
  // callers that have not, e.g. the light-cone path.
  const TLorentzVector x4_src =
      fSrcVertexValid ? fSrcVertex : fsi::SampleSRCPosition(fA, myRand);

  std::vector<FSIInputNucleon> fsi_inputs;
  fsi_inputs.reserve(2);
  auto push = [&](int pdg, const TLorentzVector &p4, int role) {
    if (pdg != pCode && pdg != nCode) return;
    FSIInputNucleon in;
    in.pdg = pdg; in.p4 = p4; in.x4 = x4_src; in.parentRole = role;
    fsi_inputs.push_back(in);
  };
  push(lead_type, vLead_target, 0);
  push(rec_type, vRec_target, 1);
  if (fsi_inputs.empty()) return;

  std::vector<FSIParticle> particles;
  bool ran = false;
  if (fIndependentFSI) {
    // Independent mode: each nucleon gets its own GENIE record and sees the
    // full undepleted A-2 remnant, regardless of what happens to the other.
    for (const auto &in : fsi_inputs) {
      std::vector<FSIInputNucleon> one(1, in);
      ran |= fsi::ApplyGenieFSI(A_res, Z_res, one, particles, myRand, fFSIModel);
    }
  } else {
    // Shared mode: both nucleons in ONE GENIE record so absorption of one
    // leg depletes the remnant seen by the other. GENIE's TransportHadrons
    // iterates the record in insertion order, so a fixed lead-first order
    // would systematically give the lead the fullest medium; both nucleons
    // leave the same vertex simultaneously and no order is physical, so the
    // transport order is always randomized.
    if (fsi_inputs.size() == 2 && myRand->Integer(2) == 1) {
      std::swap(fsi_inputs[0], fsi_inputs[1]);
    }
    ran = fsi::ApplyGenieFSI(A_res, Z_res, fsi_inputs, particles,
                             myRand, fFSIModel);
  }

  // Replace the seeded pre-FSI nucleons with GENIE's final state. If the
  // cascade could not run at all, keep the generated nucleons.
  if (ran) fParticles = particles;

  // No post-cascade Pauli veto: GENIE's cascade already Pauli-blocks each
  // collision internally (blocked outcomes are rerolled), so the final
  // state needs no further momentum cut. Matches the 3N generator.
#endif
}


// ---------------------------------------------------------------------------
// kLocalKF vertex sampling
//
// I_SRC(r) = int_{kF(r)}^{pRelmax} k^2 S(k) dk is the SRC strength a vertex at
// radius r carries once the onset follows the local Fermi momentum. It grows
// by orders of magnitude as kF(r) -> 0 at large radius, so sampling the vertex
// from rho^n r^2 and leaving I_SRC(r) in the weight concentrates the sample's
// entire weight on a few far-out vertices. Draw the vertex from
// rho^n(r) r^2 I_SRC(r) instead: the leftover weight factor is then the
// constant <I_SRC> and the vertex distribution is the physically correct one
// for a radius-dependent onset.
// ---------------------------------------------------------------------------
void QEGeneratorFSI::build_local_onset_vertex_table(int lead_type, int rec_type)
{
  const int cls = is_class(lead_type, rec_type);
  if (fLocalVtx_built[cls]) return;

  const int    n_rho = fsi::GetSRCVertexRhoPower();
  const double rmax  = 2.5 * 1.4 * TMath::Power(fA, 1./3.);   // fm, as fsi::SampleSRCPosition
  const int    nbin  = 2000;
  const double dr    = rmax / nbin;

  double fmax = 0., num = 0., den = 0.;
  for (int i = 1; i <= nbin; i++) {
    const double r   = i * dr;
    const double rho = fsi::NuclearDensity(fA, r);
    if (rho <= 0.) continue;
    const double base = TMath::Power(rho, n_rho) * r * r;    // rho^n r^2
    const double kF_r = fsi::LocalFermiMomentum(
        fA, (lead_type == pCode) ? fZ : (fA - fZ), r);
    const double I    = pRel_integral(lead_type, rec_type, kF_r, pRelmax);
    if (base * I > fmax) fmax = base * I;
    num += base * I * dr;
    den += base * dr;
  }

  fLocalVtx_rmax[cls]  = rmax;
  fLocalVtx_fmax[cls]  = fmax * 1.05;               // safety margin
  fLocalVtx_meanI[cls] = (den > 0.) ? num / den : 0.;
  fLocalVtx_built[cls] = true;
}

TLorentzVector QEGeneratorFSI::sample_local_onset_vertex(int lead_type,
                                                        int rec_type,
                                                        double &kF_r,
                                                        double &I_src_r)
{
  const int cls = is_class(lead_type, rec_type);
  build_local_onset_vertex_table(lead_type, rec_type);

  const int    n_rho = fsi::GetSRCVertexRhoPower();
  const double rmax  = fLocalVtx_rmax[cls];
  const double fmax  = fLocalVtx_fmax[cls];
  const int    N_q   = (lead_type == pCode) ? fZ : (fA - fZ);

  double r = 0.;
  kF_r = 0.; I_src_r = 0.;
  for (int guard = 0; guard < 100000; guard++) {
    r = rmax * myRand->Rndm();
    const double rho = fsi::NuclearDensity(fA, r);
    if (rho <= 0.) continue;
    const double kF = fsi::LocalFermiMomentum(fA, N_q, r);
    const double I  = pRel_integral(lead_type, rec_type, kF, pRelmax);
    const double f  = TMath::Power(rho, n_rho) * r * r * I;
    if (myRand->Rndm() * fmax <= f) { kF_r = kF; I_src_r = I; break; }
  }

  const double cosTheta = 2. * myRand->Rndm() - 1.;
  const double sinTheta = TMath::Sqrt(1. - cosTheta * cosTheta);
  const double phi = 2. * TMath::Pi() * myRand->Rndm();
  return TLorentzVector(r * sinTheta * TMath::Cos(phi),
                        r * sinTheta * TMath::Sin(phi),
                        r * cosTheta, 0.);
}

void QEGeneratorFSI::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  double Estar;
  generate_event(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,Estar);
}

void QEGeneratorFSI::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &Estar)
{
  TVector3 vRel_target;
  TLorentzVector q_target;
  generate_event(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,vRel_target,q_target,Estar);
}

void QEGeneratorFSI::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, TVector3 &vRel_target, TLorentzVector &q_target, double &Estar)
{
  fParticles.clear();
  fSrcVertexValid = false;
  fLastLocalKF = 0.;
  weight = 1.;

  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // ---- Pair vertex and the local Fermi momentum there --------------------
  // Drawn BEFORE the relative momentum so that k_F(r) can define the SRC
  // onset. ApplyFSI reuses this same vertex, so the FSI entry point is
  // unchanged -- only the order in which it is sampled.
  // N_q follows the struck nucleon's isospin (Z for a proton, A-Z for a
  // neutron), matching LocalFGM::GenerateNucleon in GCF_GENIE.
  // ---- SRC onset -------------------------------------------------------
  // kFixedKF: vertex from rho^n(r) r^2 as before, and the caller-configured
  //   [pRelmin, pRelmax] is left untouched (wf_mode 0 configures its own).
  // kLocalKF: the floor becomes kF(r), so the vertex must be drawn first --
  //   and drawn from rho^n(r) r^2 I_SRC(r), see sample_local_onset_vertex.
  double I_src_r = 0.;
  if (fSrcOnset == kLocalKF) {
    fSrcVertex = sample_local_onset_vertex(lead_type, rec_type,
                                           fLastLocalKF, I_src_r);
    fSrcVertexValid = true;
    if (fLastLocalKF <= 0. || I_src_r <= 0.) {
      static bool complained = false;
      if (!complained) {
        std::cerr << "ERROR: SRC onset 'local' could not sample a pair vertex "
                  << "(kF(r)=" << fLastLocalKF << ", I_SRC=" << I_src_r
                  << "; built without USE_FSI?). Events will be rejected."
                  << std::endl;
        complained = true;
      }
      weight = 0.; return;
    }
    if (fLastLocalKF >= pRelmax) { weight = 0.; return; }
    set_pRel_range(fLastLocalKF, pRelmax);
    set_pRel_cut(fLastLocalKF);
  } else {
    fSrcVertex = fsi::SampleSRCPosition(fA, myRand);
    fSrcVertexValid = true;
    fLastLocalKF = fsi::LocalFermiMomentum(
        fA, (lead_type == pCode) ? fZ : (fA - fZ), fSrcVertex.Vect().Mag());
  }

  double mAm2 = get_mAm2(lead_type, rec_type, Estar);

  TVector3 v1, vRec;
  decay_function(weight, lead_type, rec_type, v1, vRec);

  if (weight <= 0.) return;

  // Undo the per-vertex I_SRC(r) that decay_function put into the weight (the
  // vertex was already drawn proportional to it) and replace it by the
  // constant <I_SRC>. Net effect: expectation values unchanged, the r-driven
  // weight spread gone.
  if (fSrcOnset == kLocalKF) {
    weight *= fLocalVtx_meanI[is_class(lead_type, rec_type)] / I_src_r;
  }

  // Optional second test on the struck-nucleon momentum AFTER the c.m. motion
  // is added (v1 = P_cm/2 + k_rel). Off by default -- see the header.
  if (fSrcRequireP1 && v1.Mag() <= fLastLocalKF) { weight = 0.; return; }

  TVector3 vAm2 = - v1 - vRec;
  double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);

  vRel_target = 0.5*(v1 - vRec);

  double Erec = sqrt(sq(mN) + vRec.Mag2());
  vRec_target.SetVect(vRec);
  vRec_target.SetT(Erec);

  double E1 = mA - EAm2 - Erec;
  TLorentzVector v1_target(v1,E1);

  if (!solve_lepton_kinematics(weight, lead_type, v1, v1_target, E1,
                               vk_target, vLead_target, q_target))
    return;

  if (weight > 0.) {
    SeedFinalState(lead_type, rec_type, vLead_target, vRec_target);
    if (doFSI) ApplyFSI(lead_type, rec_type, vLead_target, vRec_target);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    if (rec_type == pCode) coulombCorrection(vRec_target, deltaECoul);
    CoulombCorrectFinalState(deltaECoul);
  }
}

// ---------------------------------------------------------------------------
// solve_lepton_kinematics — scattered-electron kinematics + Jacobian + e-N
// cross section, given the struck-nucleon momentum v1, its off-shell 4-vector
// v1_target (energy E1) and type. Extracted verbatim from the 2N pair path so
// the single-nucleon mean-field path uses identical physics. Returns false
// (weight=0) on any kinematic rejection.
bool QEGeneratorFSI::solve_lepton_kinematics(double &weight, int lead_type,
                                             const TVector3 &v1,
                                             const TLorentzVector &v1_target, double E1,
                                             TLorentzVector &vk_target,
                                             TLorentzVector &vLead_target,
                                             TLorentzVector &q_target)
{
  double p1_minus = E1 - v1.Z();

  TVector3 vbeam_int = vbeam;

  if(doCoul) {
    double deltaECoul = calcCoulombEnergy();
    vbeam_int.SetMag(vbeam_int.Mag() + deltaECoul);
  }

  vbeam_int = (doRad ? radiateElectron(vbeam_int) : vbeam_int);
  double Ebeam_int = vbeam_int.Mag();
  TLorentzVector vbeam_int_target(vbeam_int,Ebeam_int);

  double QSqmax_kine = 2*Ebeam_int*p1_minus;
  if (QSqmax_kine < QSqmin) { weight=0.; return false; }

  double QSq = QSqmin + (min(QSqmax,QSqmax_kine) - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (min(QSqmax,QSqmax_kine) - QSqmin) * (phikmax - phikmin);

  double k_minus = QSq/(2*Ebeam_int);
  double plead_minus = p1_minus - k_minus;
  if (plead_minus < 0.) { weight=0.; return false; }

  double p1_plus = E1 + v1.Z();
  double virt = v1_target.Mag2() - sq(mN);

  double A = k_minus/p1_minus;
  double c = A*(p1_plus*k_minus - 2*Ebeam_int*plead_minus - virt);

  double delta_phi = phik - v1.Phi();
  double p1_perp = v1.Perp();
  double y = p1_perp*cos(delta_phi);
  double b = -2*A*y;
  double D = sq(b) - 4*c;
  if (D < 0.) { weight=0.; return false; }

  double k_perp = (- b + sqrt(D))/2.;

  if (sqrt(D) < -b) {
    weight *= 2.;
    if (myRand->Rndm() > 0.5) k_perp = (- b - sqrt(D))/2.;
  }

  double k_plus = sq(k_perp)/k_minus;
  double kz = (k_plus - k_minus)/2.;
  TVector3 vk_int(k_perp*cos(phik),k_perp*sin(phik),kz);
  double Ek_int = vk_int.Mag();
  TLorentzVector vk_int_target;
  vk_int_target.SetVect(vk_int);
  vk_int_target.SetT(Ek_int);

  q_target = vbeam_int_target - vk_int_target;

  vLead_target = v1_target + vbeam_int_target - vk_int_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();

  TVector3 vk = (doRad ? radiateElectron(vk_int) : vk_int);
  double Ek = vk.Mag();
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);

  double J = 2.*Ebeam_int*Ek_int*fabs(1. - (v1.Z() + y * tan(vk.Theta()/2.) + Ebeam_int - Ek_int)/Elead);
  weight *= 1./J;

  weight *= myCS->sigma_eN(Ebeam_int, vk_int, vLead, (lead_type==pCode));

  if (doRad) weight *= radiationFactor(Ebeam, Ek_int, QSq);

  return true;
}

// ---------------------------------------------------------------------------
// generate_event_MF — single struck nucleon drawn from a Fermi gas, no
// correlated recoil, transported through the A-1 residual.
void QEGeneratorFSI::generate_event_MF(double &weight, int &lead_type, int &rec_type,
                                       TLorentzVector& vk_target,
                                       TLorentzVector &vLead_target,
                                       TLorentzVector &vRec_target,
                                       TLorentzVector &vAm2_target)
{
  fParticles.clear();
  weight = 1.;

  // One struck nucleon; pick its isospin (factor 2 — no partner to choose).
  lead_type = (myRand->Rndm() > 0.5) ? pCode : nCode;
  rec_type  = kNoRecoil;   // no correlated recoil nucleon in the mean field
  weight *= 2.;

  // Struck-nucleon position AND momentum from the Fermi gas (Algorithm A in
  // LFG mode: sample r first, then |p| from the local Fermi sphere at that r).
  TVector3 v1;
  if (!sample_FG_momentum(lead_type, weight, v1)) { weight = 0.; return; }

  // A-1 residual ground-state mass; guard nuclei without tabulated A-1 mass.
  double mAm1 = get_mAm1(lead_type);
  if (mAm1 <= 0.) { weight = 0.; return; }

  // Off-shell struck energy: the A-1 residual recoils with -p1.
  double EAm1 = sqrt(v1.Mag2() + sq(mAm1));
  double E1 = mA - EAm1;
  TLorentzVector v1_target(v1, E1);

  // A-1 residual 4-vector (spectator).
  vAm2_target.SetVect(-v1);
  vAm2_target.SetT(EAm1);

  // No recoil nucleon in the final state.
  vRec_target.SetVect(TVector3(0.,0.,0.));
  vRec_target.SetT(0.);

  TLorentzVector q_local;
  if (!solve_lepton_kinematics(weight, lead_type, v1, v1_target, E1,
                               vk_target, vLead_target, q_local))
    return;

  if (weight > 0.) {
    SeedFinalState(lead_type, kNoRecoil, vLead_target, vRec_target);
    if (doFSI) ApplyFSI_single(lead_type, vLead_target);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    // rec_type is the no-recoil sentinel: nothing to correct.
    CoulombCorrectFinalState(deltaECoul);
  }
}

// ---------------------------------------------------------------------------
// sample_FG_momentum — Algorithm A in local-FG mode: sample the struck
// nucleon position r from rho^1(r), compute the per-isospin local Fermi
// momentum k_F^q(r), then draw |p| uniformly in [0, k_F^q(r)] and weight by
// the local density n_local(p|r) = 3/(4 pi k_F^q(r)^3) so that r and p
// carry the physically correct joint distribution.  In global-FG mode
// r is sampled independently and |p| is uniform in [0, fkF] with weight
// n_global(p) = 3/(4 pi fkF^3).  The position is cached in fLastMFPosition
// so ApplyFSI_single can reuse it as the FSI entry point.
bool QEGeneratorFSI::sample_FG_momentum(int nucleon_type, double &weight,
                                        TVector3 &v1)
{
  // Position (rho^1-weighted, isotropic). Used as the FSI entry point.
  fLastMFPosition = fsi::SampleMFPosition(fA, myRand);

  // Per-isospin local Fermi momentum at the sampled r.
  double kmax = 0.;
  if (fFGMode == kLocalFG) {
    const double r_fm = fLastMFPosition.Vect().Mag();
    const int N_q = (nucleon_type == pCode) ? fZ : (fA - fZ);
    kmax = fsi::LocalFermiMomentum(fA, N_q, r_fm);
    // Lazily build the marginal table for eval_n_FG / get_n_FG bookkeeping.
    if (!fLocalFG_built) build_local_FG_table();
  } else {
    kmax = fkF;
  }
  if (kmax <= 0.) return false;

  // Isotropic direction, |p| uniform in [0, kmax].
  const double phi      = 2.*M_PI*myRand->Rndm();
  const double cosTheta = -1. + 2.*myRand->Rndm();
  const double theta    = acos(cosTheta);
  const double p        = kmax * myRand->Rndm();
  v1.SetMagThetaPhi(p, theta, phi);

  // Local (conditional) density at the sampled (r, p): a uniform Fermi
  // sphere of radius kmax has n_local(p) = 3 / (4 pi kmax^3) for p<=kmax.
  const double n_local = 3.0 / (4.0*M_PI*pow(kmax,3));

  // Phase space (dphi=2pi, dcos=2, dp=kmax) * p^2 * n_local, times the number
  // of nucleons of this isospin. n_local (and the local-FG joint (r, p)
  // density) is normalized to ONE nucleon; multiplying by N_q puts the
  // mean-field weight on the same footing as the SRC weight, whose
  // k^2 S(k)/(2 pi)^3 integrates to the number of pairs: both then read
  // "number of nucleons (pairs) x density x cross section" per attempt.
  const int N_q = (nucleon_type == pCode) ? fZ : (fA - fZ);
  weight *= (2.*M_PI) * 2. * kmax * sq(p) * n_local * N_q;

  // Record the actual per-event density used in the weight.
  fLastMFRho = n_local;
  return weight > 0.;
}

// ---------------------------------------------------------------------------
// eval_n_FG — Fermi-gas momentum density n_FG(p), normalised to int n d^3p = 1.
double QEGeneratorFSI::eval_n_FG(double p) const
{
  if (fFGMode == kGlobalFG) {
    if (p > fkF) return 0.;
    return 3.0 / (4.0*M_PI*fkF*fkF*fkF);
  }
  // Local FG: interpolate the cached table.
  if (fLocalFG_p.empty()) return 0.;
  if (p <= fLocalFG_p.front()) return fLocalFG_n.front();
  if (p >= fLocalFG_p.back())  return 0.;
  int lo = 0, hi = static_cast<int>(fLocalFG_p.size()) - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi)/2;
    if (fLocalFG_p[mid] <= p) lo = mid; else hi = mid;
  }
  double t = (p - fLocalFG_p[lo]) / (fLocalFG_p[hi] - fLocalFG_p[lo]);
  return fLocalFG_n[lo] + t*(fLocalFG_n[hi] - fLocalFG_n[lo]);
}

// ---------------------------------------------------------------------------
// build_local_FG_table — tabulate the local-Fermi-gas momentum distribution
//   n(p) = (1/N) int d^3r Theta(kF(r) - p),  kF(r) = (3 pi^2 (A/2) rho(r))^{1/3}
// (per-isospin kF for a symmetric nucleus: GENIE's rho is normalised to 1,
// so the species-q nucleon density is N_q rho ~ (A/2) rho -- the same
// convention as fsi::LocalFermiMomentum used by the event sampler).
// with rho from GENIE's density model, normalised so int n d^3p = 1. Cached
// per nucleus; the momentum ceiling is kF at central density.
void QEGeneratorFSI::build_local_FG_table()
{
  const double hbarc = 0.1973269804;   // GeV*fm
  const int    nr    = 400;
  const double Rmax  = 2.5 * 1.4 * pow((double)fA, 1./3.);  // fm (matches sampler)
  const double dr    = Rmax / nr;

  std::vector<double> r(nr), kF(nr);
  double kFmax = 0.;
  for (int i = 0; i < nr; i++) {
    r[i] = (i + 0.5) * dr;
    double rho = fsi::NuclearDensity(fA, r[i]);   // fm^-3
    if (rho < 0.) rho = 0.;
    double kF_fm = pow(3.0*M_PI*M_PI*(0.5*fA)*rho, 1./3.);  // (3 pi^2 (A/2) rho)^{1/3} [fm^-1]
    kF[i] = hbarc * kF_fm;                          // GeV/c
    if (kF[i] > kFmax) kFmax = kF[i];
  }
  fLocalFG_kmax = kFmax;

  // Normalisation N = int d^3r (4pi/3) kF(r)^3.
  double norm = 0.;
  for (int i = 0; i < nr; i++)
    norm += 4.*M_PI*r[i]*r[i]*dr * (4.*M_PI/3.)*pow(kF[i],3);

  const int np = 200;
  fLocalFG_p.assign(np, 0.);
  fLocalFG_n.assign(np, 0.);
  for (int j = 0; j < np; j++) {
    double p = (np > 1) ? kFmax * j/(np-1) : 0.;
    fLocalFG_p[j] = p;
    double g = 0.;                                  // g(p) = int d^3r Theta(kF-p)
    for (int i = 0; i < nr; i++)
      if (kF[i] > p) g += 4.*M_PI*r[i]*r[i]*dr;
    fLocalFG_n[j] = (norm > 0.) ? g/norm : 0.;
  }
  fLocalFG_built = true;
}

// ---------------------------------------------------------------------------
// ApplyFSI_single — transport ONE struck nucleon through the A-1 residual and
// store GENIE's flat final-state particle list in fParticles. The generated
// type/momentum and the event weight are not modified.
void QEGeneratorFSI::ApplyFSI_single(int lead_type,
                                     const TLorentzVector &vLead_target)
{
  if (lead_type != pCode && lead_type != nCode) return;

  // Residual nucleus after removing the single struck nucleon (A-1).
  const int removedP = (lead_type == pCode ? 1 : 0);
  const int removedN = (lead_type == nCode ? 1 : 0);
  const int Z_res = fZ - removedP;
  const int N_res = (fA - fZ) - removedN;
  const int A_res = Z_res + N_res;
  if (A_res < 3 || Z_res < 0 || N_res < 0) return;

#ifdef USE_FSI
  // Reuse the position sampled in sample_FG_momentum so the (r, p)
  // correlation set by the local Fermi gas is preserved in the cascade.
  // In global-FG mode the position is independent of |p| (still rho^1).
  std::vector<FSIInputNucleon> one(1);
  one[0].pdg = lead_type;
  one[0].p4 = vLead_target;
  one[0].x4 = fLastMFPosition;
  one[0].parentRole = 0;

  std::vector<FSIParticle> particles;
  if (fsi::ApplyGenieFSI(A_res, Z_res, one, particles, myRand, fFSIModel)) {
    fParticles = particles;
  }

  // No post-cascade Pauli veto here either: GENIE's cascade handles Pauli
  // blocking internally (see ApplyFSI above).
#endif
}

void QEGeneratorFSI::generate_event_lightcone(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  double Estar;
  generate_event_lightcone(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,Estar);
}

void QEGeneratorFSI::generate_event_lightcone(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &Estar)
{
  fParticles.clear();
  weight = 1.;

  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  double mAm2 = get_mAm2(lead_type, rec_type, Estar);

  double alpha1, alphaRec;
  TVector2 v1_perp, vRec_perp;
  decay_function_lc(weight, lead_type, rec_type, alpha1, v1_perp, alphaRec, vRec_perp);

  if (weight <= 0.) return;

  double alphaCM = alpha1 + alphaRec;
  double alphaAm2 = Anum - alphaCM;
  TVector2 vAm2_perp = -1.*v1_perp - vRec_perp;

  TVector3 vbeam_int = vbeam;

  if(doCoul) {
    double deltaE = calcCoulombEnergy();
    vbeam_int.SetMag(vbeam_int.Mag() + deltaE);
  }

  vbeam_int = (doRad ? radiateElectron(vbeam_int) : vbeam_int);
  double Ebeam_int = vbeam_int.Mag();
  TLorentzVector vbeam_int_target(vbeam_int,Ebeam_int);

  double QSq = QSqmin + (QSqmax - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (QSqmax - QSqmin) * (phikmax - phikmin);

  double p1_minus = mbar*alpha1;
  double pRec_minus = mbar*alphaRec;
  double pAm2_minus = mbar*alphaAm2;
  double pRec_plus = (sq(mN) + vRec_perp.Mod2())/pRec_minus;
  double pAm2_plus = (sq(mAm2) + vAm2_perp.Mod2())/pAm2_minus;
  if (mAm2 == 0) pAm2_plus = 0.;
  double p1_plus = mA - pRec_plus - pAm2_plus;
  double virt = p1_plus*p1_minus- v1_perp.Mod2() - sq(mN);

  double E1 = 0.5*(p1_plus + p1_minus);
  double p1_z = 0.5*(p1_plus - p1_minus);
  double A = 0.5*(QSq - virt);
  double a = p1_plus*p1_minus;
  double b = -2*E1*A;
  double c = sq(A) - QSq*sq(p1_z);
  double D = sq(b) - 4*a*c;

  if (D<0.) { weight = 0.; return; }

  double omega1 = (-b - sqrt(D))/(2.*a);
  double omega2 = (-b + sqrt(D))/(2.*a);
  double omega;

  bool omega1Valid = (omega1>=0.) && (omega1<=Ebeam_int) && ((omega1*E1-A)*p1_z > 0);
  bool omega2Valid = (omega2>=0.) && (omega2<=Ebeam_int) && ((omega2*E1-A)*p1_z > 0);

  if ((!omega1Valid) && (!omega2Valid)) { weight=0.; return; }
  if (!omega1Valid) omega = omega2;
  else if (!omega2Valid) omega = omega1;
  else { omega = (gRandom->Rndm()>0.5)? omega1 : omega2; weight *= 2.; }

  double Ek_int = Ebeam_int - omega;
  double cosThetak = 1. - QSq/(2.*Ebeam_int*Ek_int);

  if (fabs(cosThetak) > 1.) { weight=0.; return; }

  double sinThetak = sqrt(1. - sq(cosThetak));
  double kz = Ek_int*cosThetak;
  double kperp = Ek_int*sinThetak;
  TVector3 vk_int(kperp*cos(phik),kperp*sin(phik),kz);

  TLorentzVector vk_int_target;
  vk_int_target.SetVect(vk_int);
  vk_int_target.SetT(Ek_int);

  TVector3 vq = vbeam_int - vk_int;
  double rot_phi = vq.Phi();
  double rot_theta = vq.Theta();

  TVector3 v1(v1_perp.X(),v1_perp.Y(),p1_z);
  v1.RotateY(rot_theta);
  v1.RotateZ(rot_phi);
  TLorentzVector v1_target(v1,E1);

  double ERec = 0.5*(pRec_plus + pRec_minus);
  double pRec_z = 0.5*(pRec_plus - pRec_minus);
  TVector3 vRec(vRec_perp.X(),vRec_perp.Y(),pRec_z);
  vRec.RotateY(rot_theta);
  vRec.RotateZ(rot_phi);
  vRec_target.SetVect(vRec);
  vRec_target.SetT(ERec);

  double EAm2 = 0.5*(pAm2_plus + pAm2_minus);
  double pAm2_z = 0.5*(pAm2_plus - pAm2_minus);
  TVector3 vAm2(vAm2_perp.X(),vAm2_perp.Y(),pAm2_z);
  vAm2.RotateY(rot_theta);
  vAm2.RotateZ(rot_phi);
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);

  vLead_target = v1_target + vbeam_int_target - vk_int_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();

  TVector3 vk = (doRad ? radiateElectron(vk_int) : vk_int);
  double Ek = vk.Mag();
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);

  double J = 2.*Ebeam_int*Ek_int*fabs(1. - (vLead.Dot(vq)*omega)/(vq.Mag2()*Elead));
  weight *= 1./J;

  weight *= myCS->sigma_eN(Ebeam_int, vk_int, vLead, (lead_type==pCode));

  if (doRad) weight *= radiationFactor(Ebeam, Ek_int, QSq);

  if (weight > 0.) {
    SeedFinalState(lead_type, rec_type, vLead_target, vRec_target);
    if (doFSI) ApplyFSI(lead_type, rec_type, vLead_target, vRec_target);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    if (rec_type == pCode) coulombCorrection(vRec_target, deltaECoul);
    CoulombCorrectFinalState(deltaECoul);
  }
}
