#include "QEGeneratorFSI.hh"
#include "constants.hh"
#include "helpers.hh"
#include "fsi/GenieFSIHelpers.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
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
  fFermiMomentum = 0.220; // GeV/c — typical for medium/heavy nuclei (C-12)

  fsi::ResolveGenieXMLPath();
}

QEGeneratorFSI::~QEGeneratorFSI()
{
}

void QEGeneratorFSI::SetFSITuning(double fermiMom_MeV)
{
  fFermiMomentum = fermiMom_MeV * 1.e-3; // convert MeV/c -> GeV/c
}

bool QEGeneratorFSI::CheckPauliBlocking(const TLorentzVector &p4) const
{
  return (p4.P() < fFermiMomentum);
}

// ---------------------------------------------------------------------------
// ApplyElasticScatter
// Perform a proper 2-body elastic scatter of p4 off a target nucleon at rest.
// Conserves 4-momentum exactly. Scattering angle is sampled isotropically in
// the CM frame — the same approximation used in GENIE's ElasHN when the
// nucleon-nucleon angular distribution is not separately parameterised.
// On return p4 contains the updated 4-momentum of the scattered nucleon.
void QEGeneratorFSI::ApplyElasticScatter(TLorentzVector &p4)
{
  // Target nucleon at rest in the lab frame
  TLorentzVector p4_target(0., 0., 0., mN);

  // Total 4-momentum
  TLorentzVector p_tot = p4 + p4_target;

  // Boost vector to the CM frame
  TVector3 beta = p_tot.BoostVector();

  // Boost the incident nucleon to CM
  TLorentzVector p4_cm = p4;
  p4_cm.Boost(-beta);

  double p_cm = p4_cm.P();
  double E_cm = p4_cm.E();

  // Sample an isotropic direction in the CM frame
  double cosTheta = 2.0 * myRand->Rndm() - 1.0;
  double sinTheta = TMath::Sqrt(1.0 - cosTheta * cosTheta);
  double phi      = 2.0 * TMath::Pi() * myRand->Rndm();

  TVector3 p_cm_new(p_cm * sinTheta * TMath::Cos(phi),
                    p_cm * sinTheta * TMath::Sin(phi),
                    p_cm * cosTheta);

  TLorentzVector p4_cm_new(p_cm_new, E_cm);

  // Boost back to the lab frame
  p4_cm_new.Boost(beta);
  p4 = p4_cm_new;
}

// ---------------------------------------------------------------------------
// ApplyFSI — main FSI dispatcher
//
// Uses GENIE's intranuclear cascade (hN or hA mode) to transport each
// outgoing nucleon through the residual nucleus. Pauli-blocks the final
// state if either nucleon momentum < p_Fermi.
void QEGeneratorFSI::ApplyFSI(int &lead_type, int &rec_type,
                               TLorentzVector &vLead_target,
                               TLorentzVector &vRec_target,
                               double &weight)
{
  fLastFSIEventStats = FSIEventStats();
  fLastFSISecondaries.clear();
  fLeadPreFSI = vLead_target;
  fRecPreFSI  = vRec_target;
  fLeadAbsorbed = false;
  fRecAbsorbed  = false;
  if (!doFSI || weight <= 0.) return;

  const int lead_type_in = lead_type;
  const int rec_type_in = rec_type;

  // Residual nucleus after removing the struck/knocked-out SRC pair.
  const int removedP = (lead_type_in == pCode ? 1 : 0) + (rec_type_in == pCode ? 1 : 0);
  const int removedN = (lead_type_in == nCode ? 1 : 0) + (rec_type_in == nCode ? 1 : 0);
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

  // Always transport through the actual A-2 residual nucleus.
  const int A_transport = A_res;
  const int Z_transport = Z_res;

#ifdef USE_FSI
  const TLorentzVector vLead_before = vLead_target;
  const TLorentzVector vRec_before  = vRec_target;

  // Sample position once for the SRC pair (both nucleons originate from the
  // same location).  Weighted by rho^2 of the FULL nucleus, since the pair
  // was formed before knockout.
  const TLorentzVector x4_src = fsi::SampleSRCPosition(fA, myRand);

  // Transport BOTH outgoing nucleons through a SHARED, depleting remnant in a
  // single GENIE call. Previously each nucleon got its own GHepRecord, so
  // nucleon 2 saw an undepleted A-2 nucleus even when nucleon 1 had been
  // absorbed during its cascade -- biasing FSI rates high. With a single
  // record GENIE's Intranuke2018::TransportHadrons walks the pair through
  // one fRemnA/fRemnZ that decrements correctly across the two cascades.
  std::vector<fsi::FSIInputNucleon> fsi_inputs;
  fsi_inputs.reserve(2);
  if (lead_type == pCode || lead_type == nCode) {
    fsi::FSIInputNucleon in;
    in.pdg = lead_type;
    in.p4 = vLead_target;
    in.x4 = x4_src;
    in.parentRole = 0;
    fsi_inputs.push_back(in);
  }
  if (rec_type == pCode || rec_type == nCode) {
    fsi::FSIInputNucleon in;
    in.pdg = rec_type;
    in.p4 = vRec_target;
    in.x4 = x4_src;
    in.parentRole = 1;
    fsi_inputs.push_back(in);
  }

  std::vector<fsi::FSIOutputNucleon> fsi_outputs;
  if (fIndependentFSI) {
    // Independent mode: each nucleon gets its own GENIE record and sees the
    // full undepleted A-2 remnant, regardless of what happens to the other.
    fsi_outputs.resize(fsi_inputs.size());
    for (size_t i = 0; i < fsi_inputs.size(); ++i) {
      int pdg = fsi_inputs[i].pdg;
      TLorentzVector p4 = fsi_inputs[i].p4;
      fsi_outputs[i].survived =
        fsi::ApplyGenieFSIToNucleon(A_transport, Z_transport, pdg, p4,
                                    fsi_inputs[i].x4, myRand,
                                    fsi_inputs[i].parentRole,
                                    fLastFSISecondaries, fFSIModel);
      fsi_outputs[i].pdg = pdg;
      fsi_outputs[i].p4  = p4;
    }
  } else {
    fsi::ApplyGenieFSIToNucleons(A_transport, Z_transport,
                                  fsi_inputs, fsi_outputs,
                                  myRand, fLastFSISecondaries, fFSIModel);
  }

  // Demux outputs back into lead/recoil.
  //
  // Absorption (fate = abs/cmp) does NOT kill the event: the GCF
  // primary-vertex cross section is still valid, and any surviving partner /
  // FSI secondaries are still physically observable. The absorbed flag
  // records that the ORIGINAL nucleon ceased to exist in the cascade, but
  // the experiment only sees the final state: if the absorption ejected a
  // nucleon, the detector would reconstruct THAT as the leg. So the leg's
  // saved momentum becomes its leading nucleon ejecta (removed from the
  // secondaries so multiplicity counts don't double-count it); only if no
  // nucleon ejecta exists is the leg zeroed.
  auto promote = [&](int role, int &pdg_out, TLorentzVector &p4_out) -> bool {
    int best = -1; double bestP = -1.;
    for (size_t i = 0; i < fLastFSISecondaries.size(); ++i) {
      if (fLastFSISecondaries[i].parentRole != role) continue;
      if (fLastFSISecondaries[i].pdg != pCode &&
          fLastFSISecondaries[i].pdg != nCode) continue;
      if (fLastFSISecondaries[i].p4.P() > bestP) {
        bestP = fLastFSISecondaries[i].p4.P(); best = (int)i;
      }
    }
    if (best < 0) return false;
    pdg_out = fLastFSISecondaries[best].pdg;
    p4_out  = fLastFSISecondaries[best].p4;
    fLastFSISecondaries.erase(fLastFSISecondaries.begin() + best);
    return true;
  };
  size_t out_idx = 0;
  if (lead_type == pCode || lead_type == nCode) {
    if (out_idx >= fsi_outputs.size() || !fsi_outputs[out_idx].survived) {
      fLeadAbsorbed = true;
      if (!promote(0, lead_type, vLead_target)) {
        lead_type = 0;
        vLead_target = TLorentzVector(0., 0., 0., 0.);
      }
    } else {
      lead_type = fsi_outputs[out_idx].pdg;
      vLead_target = fsi_outputs[out_idx].p4;
    }
    ++out_idx;
  }
  if (rec_type == pCode || rec_type == nCode) {
    if (out_idx >= fsi_outputs.size() || !fsi_outputs[out_idx].survived) {
      fRecAbsorbed = true;
      if (!promote(1, rec_type, vRec_target)) {
        rec_type = 0;
        vRec_target = TLorentzVector(0., 0., 0., 0.);
      }
    } else {
      rec_type = fsi_outputs[out_idx].pdg;
      vRec_target = fsi_outputs[out_idx].p4;
    }
    ++out_idx;
  }

  fLastFSIEventStats.leadChargeExchange =
      ((lead_type_in == pCode || lead_type_in == nCode) &&
       (lead_type == pCode || lead_type == nCode) &&
       (lead_type != lead_type_in));
  fLastFSIEventStats.recoilChargeExchange =
      ((rec_type_in == pCode || rec_type_in == nCode) &&
       (rec_type == pCode || rec_type == nCode) &&
       (rec_type != rec_type_in));

  int nLeadNonNucleonSec = 0;
  int nRecNonNucleonSec = 0;
  for (const auto &sec : fLastFSISecondaries) {
    const bool isNucleon = (sec.pdg == pCode || sec.pdg == nCode);
    if (sec.parentRole == 0 && !isNucleon) nLeadNonNucleonSec++;
    if (sec.parentRole == 1 && !isNucleon) nRecNonNucleonSec++;
  }

  const double leadDeltaP = (vLead_target.Vect() - vLead_before.Vect()).Mag();
  const double recDeltaP  = (vRec_target.Vect()  - vRec_before.Vect()).Mag();
  const bool leadChanged = (leadDeltaP > 1e-6);
  const bool recChanged  = (recDeltaP  > 1e-6);

  // "Elastic-like": momentum changed (interaction occurred), no charge exchange,
  // and no non-nucleon secondaries (pions, etc.).  In GENIE's hN cascade an
  // elastic N-N scatter always ejects the struck nucleon as a secondary, so
  // requiring zero secondaries would make this condition unreachable.
  // Absorbed legs are excluded -- their post-FSI p4 is zero by convention,
  // which would otherwise look like a huge elastic kick.
  fLastFSIEventStats.leadElasticLike =
      !fLeadAbsorbed &&
      !fLastFSIEventStats.leadChargeExchange &&
      (nLeadNonNucleonSec == 0) &&
      leadChanged;
  fLastFSIEventStats.recoilElasticLike =
      !fRecAbsorbed &&
      !fLastFSIEventStats.recoilChargeExchange &&
      (nRecNonNucleonSec == 0) &&
      recChanged;

  fLastFSIEventStats.nSecondaries = static_cast<int>(fLastFSISecondaries.size());
  for (const auto &sec : fLastFSISecondaries) {
    if (sec.pdg == 211) {
      fLastFSIEventStats.nPiPlus++;
      fLastFSIEventStats.nPionsTotal++;
    } else if (sec.pdg == -211) {
      fLastFSIEventStats.nPiMinus++;
      fLastFSIEventStats.nPionsTotal++;
    } else if (sec.pdg == 111) {
      fLastFSIEventStats.nPiZero++;
      fLastFSIEventStats.nPionsTotal++;
    }
  }

  // No post-cascade Pauli veto: GENIE's cascade already Pauli-blocks each
  // collision internally (blocked outcomes are rerolled), so the surviving
  // final state needs no further momentum cut. Matches the 3N generator.
#endif
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
  fLastFSIEventStats = FSIEventStats();
  fLastFSISecondaries.clear();
  // Also clear the pre-FSI cache and absorbed flags: a kinematic rejection
  // below returns before ApplyFSI* runs, and consumers use
  // GetPreFSILead().E()==0 to detect that FSI never ran this event.
  fLeadPreFSI = TLorentzVector();
  fRecPreFSI  = TLorentzVector();
  fLeadAbsorbed = false;
  fRecAbsorbed  = false;
  weight = 1.;

  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  double mAm2 = get_mAm2(lead_type, rec_type, Estar);

  TVector3 v1, vRec;
  decay_function(weight, lead_type, rec_type, v1, vRec);

  if (weight <= 0.) return;

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

  if (doFSI && weight > 0.) {
    ApplyFSI(lead_type, rec_type, vLead_target, vRec_target, weight);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    if (rec_type == pCode) coulombCorrection(vRec_target, deltaECoul);
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
  fLastFSIEventStats = FSIEventStats();
  fLastFSISecondaries.clear();
  // Also clear the pre-FSI cache and absorbed flags: a kinematic rejection
  // below returns before ApplyFSI* runs, and consumers use
  // GetPreFSILead().E()==0 to detect that FSI never ran this event.
  fLeadPreFSI = TLorentzVector();
  fRecPreFSI  = TLorentzVector();
  fLeadAbsorbed = false;
  fRecAbsorbed  = false;
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

  if (doFSI && weight > 0.) {
    ApplyFSI_single(lead_type, vLead_target, weight);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    // rec_type is the no-recoil sentinel: nothing to correct.
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

  // Phase space (dphi=2pi, dcos=2, dp=kmax) * p^2 * n_local / (2 pi)^3.
  weight *= (2.*M_PI) * 2. * kmax * sq(p) * n_local / pow(2.*M_PI,3);

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
// ApplyFSI_single — transport ONE struck nucleon through the A-1 residual.
void QEGeneratorFSI::ApplyFSI_single(int &lead_type, TLorentzVector &vLead_target,
                                     double &weight)
{
  fLastFSIEventStats = FSIEventStats();
  fLastFSISecondaries.clear();
  fLeadPreFSI = vLead_target;
  fRecPreFSI  = TLorentzVector(0.,0.,0.,0.);
  fLeadAbsorbed = false;
  fRecAbsorbed  = false;  // no recoil leg in MF mode -- never absorbed
  if (!doFSI || weight <= 0.) return;

  const int lead_type_in = lead_type;

  // Residual nucleus after removing the single struck nucleon (A-1).
  const int removedP = (lead_type_in == pCode ? 1 : 0);
  const int removedN = (lead_type_in == nCode ? 1 : 0);
  const int Z_res = fZ - removedP;
  const int N_res = (fA - fZ) - removedN;
  const int A_res = Z_res + N_res;
  if (A_res < 3 || Z_res < 0 || N_res < 0) return;

#ifdef USE_FSI
  const TLorentzVector vLead_before = vLead_target;
  // Reuse the position sampled in sample_FG_momentum so the (r, p)
  // correlation set by the local Fermi gas is preserved in the cascade.
  // In global-FG mode the position is independent of |p| (still rho^1).
  const TLorentzVector x4 = fLastMFPosition;

  if (lead_type == pCode || lead_type == nCode) {
    if (!fsi::ApplyGenieFSIToNucleon(A_res, Z_res, lead_type, vLead_target, x4,
                                     myRand, 0, fLastFSISecondaries, fFSIModel)) {
      // Absorbed: keep the event with its GCF primary weight, mark the leg,
      // and (experiment-faithful) save the leading nucleon ejecta as the
      // measured nucleon; zero only if no nucleon ejecta exists.
      fLeadAbsorbed = true;
      int best = -1; double bestP = -1.;
      for (size_t i = 0; i < fLastFSISecondaries.size(); ++i) {
        if (fLastFSISecondaries[i].parentRole != 0) continue;
        if (fLastFSISecondaries[i].pdg != pCode &&
            fLastFSISecondaries[i].pdg != nCode) continue;
        if (fLastFSISecondaries[i].p4.P() > bestP) {
          bestP = fLastFSISecondaries[i].p4.P(); best = (int)i;
        }
      }
      if (best >= 0) {
        lead_type = fLastFSISecondaries[best].pdg;
        vLead_target = fLastFSISecondaries[best].p4;
        fLastFSISecondaries.erase(fLastFSISecondaries.begin() + best);
      } else {
        lead_type = 0;
        vLead_target = TLorentzVector(0., 0., 0., 0.);
      }
    }
  }

  fLastFSIEventStats.leadChargeExchange =
      ((lead_type_in == pCode || lead_type_in == nCode) &&
       (lead_type == pCode || lead_type == nCode) &&
       (lead_type != lead_type_in));

  int nLeadNonNucleonSec = 0;
  for (const auto &sec : fLastFSISecondaries) {
    const bool isNucleon = (sec.pdg == pCode || sec.pdg == nCode);
    if (sec.parentRole == 0 && !isNucleon) nLeadNonNucleonSec++;
  }
  const double leadDeltaP = (vLead_target.Vect() - vLead_before.Vect()).Mag();
  fLastFSIEventStats.leadElasticLike =
      !fLeadAbsorbed &&
      !fLastFSIEventStats.leadChargeExchange &&
      (nLeadNonNucleonSec == 0) &&
      (leadDeltaP > 1e-6);

  fLastFSIEventStats.nSecondaries = static_cast<int>(fLastFSISecondaries.size());
  for (const auto &sec : fLastFSISecondaries) {
    if (sec.pdg == 211)       { fLastFSIEventStats.nPiPlus++;  fLastFSIEventStats.nPionsTotal++; }
    else if (sec.pdg == -211) { fLastFSIEventStats.nPiMinus++; fLastFSIEventStats.nPionsTotal++; }
    else if (sec.pdg == 111)  { fLastFSIEventStats.nPiZero++;  fLastFSIEventStats.nPionsTotal++; }
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
  fLastFSIEventStats = FSIEventStats();
  fLastFSISecondaries.clear();
  // Also clear the pre-FSI cache and absorbed flags: a kinematic rejection
  // below returns before ApplyFSI* runs, and consumers use
  // GetPreFSILead().E()==0 to detect that FSI never ran this event.
  fLeadPreFSI = TLorentzVector();
  fRecPreFSI  = TLorentzVector();
  fLeadAbsorbed = false;
  fRecAbsorbed  = false;
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

  if (doFSI && weight > 0.) {
    ApplyFSI(lead_type, rec_type, vLead_target, vRec_target, weight);
  }

  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);
    if (lead_type == pCode) coulombCorrection(vLead_target, deltaECoul);
    if (rec_type == pCode) coulombCorrection(vRec_target, deltaECoul);
  }
}
