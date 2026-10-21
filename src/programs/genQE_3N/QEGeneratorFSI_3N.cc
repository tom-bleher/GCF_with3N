#include "QEGeneratorFSI_3N.hh"
#include "constants.hh"
#include <iostream>

QEGeneratorFSI_3N::QEGeneratorFSI_3N(double E, eNCrossSection *thisCS,
                                     int thisU, TRandom3 *thisRand,
                                     int A, int Z)
  : QEGenerator_3N(E, thisCS, thisU, thisRand),
    fRand(thisRand),
    fDoFSI(true),
    fFSIModel(kHN2018),
    fA(A),
    fZ(Z)
{
  fsi::ResolveGenieXMLPath();
}

QEGeneratorFSI_3N::~QEGeneratorFSI_3N() {}

void QEGeneratorFSI_3N::SetTargetNucleus(int A, int Z)
{
  QEGenerator_3N::SetTargetNucleus(A, Z);  // validates + sets kinematic mA/mAm3
  fA = A;
  fZ = Z;
}

void QEGeneratorFSI_3N::generate_event_with_FSI(double &weight,
                                                int &N1_type, int &N2_type, int &N3_type,
                                                TLorentzVector &vk_target,
                                                TLorentzVector &vLead_target,
                                                TLorentzVector &v2_target,
                                                TLorentzVector &v3_target,
                                                TLorentzVector &vAm3_target,
                                                bool use_CM)
{
  // Generate the bare 3N event using the parent class.
  QEGenerator_3N::generate_event(weight, N1_type, N2_type, N3_type,
                                 vk_target, vLead_target, v2_target,
                                 v3_target, vAm3_target, use_CM);

  fParticles.clear();
  if (weight <= 0.) return;

  // Default final state: the three generated nucleons unchanged. ApplyFSI
  // replaces this with GENIE's cascade output when FSI actually runs.
  auto seed = [&](int pdg, const TLorentzVector &p4, int role) {
    FSIParticle part;
    part.pdg = pdg; part.p4 = p4; part.rescatterCode = -1; part.parentRole = role;
    fParticles.push_back(part);
  };
  seed(N1_type, vLead_target, 0);
  seed(N2_type, v2_target, 1);
  seed(N3_type, v3_target, 2);

  if (fDoFSI) {
    ApplyFSI(N1_type, N2_type, N3_type, vLead_target, v2_target, v3_target);
  }
}

void QEGeneratorFSI_3N::ApplyFSI(int N1_type, int N2_type, int N3_type,
                                 const TLorentzVector &v1,
                                 const TLorentzVector &v2,
                                 const TLorentzVector &v3)
{
  // Residual nucleus: remove the three knocked-out nucleons.
  const int removedP = (N1_type == pCode ? 1 : 0)
                     + (N2_type == pCode ? 1 : 0)
                     + (N3_type == pCode ? 1 : 0);
  const int removedN = (N1_type == nCode ? 1 : 0)
                     + (N2_type == nCode ? 1 : 0)
                     + (N3_type == nCode ? 1 : 0);
  const int Z_res = fZ - removedP;
  const int N_res = (fA - fZ) - removedN;
  const int A_res = Z_res + N_res;

  // GENIE intranuke is not validated for A<3; skip FSI rather than risk NaN.
  // The particle list keeps the generated nucleons in that case.
  if (A_res < 3 || Z_res < 0 || N_res < 0) {
    return;
  }

  if (A_res < 7) {
    static bool warned = false;
    if (!warned) {
      std::cerr << "Warning: GENIE intranuke applied to light residual nucleus "
                << "(A_res=" << A_res << "). GENIE is validated for A>6; "
                << "results for A_res<7 should be treated with caution."
                << std::endl;
      warned = true;
    }
  }

#ifdef USE_FSI
  // Single shared position for the three nucleons of the SRC triplet:
  // rho^3-weighted (one factor of rho per nucleon at the same point).
  const TLorentzVector x4_src = fsi::SampleSRCPosition3N(fA, fRand);

  std::vector<FSIInputNucleon> fsi_inputs;
  fsi_inputs.reserve(3);
  auto push = [&](int pdg, const TLorentzVector &p4, int role) {
    if (pdg != pCode && pdg != nCode) return;
    FSIInputNucleon in;
    in.pdg = pdg; in.p4 = p4; in.x4 = x4_src; in.parentRole = role;
    fsi_inputs.push_back(in);
  };
  push(N1_type, v1, 0);
  push(N2_type, v2, 1);
  push(N3_type, v3, 2);
  if (fsi_inputs.empty()) return;

  std::vector<FSIParticle> particles;
  bool ran = false;
  if (fIndependentFSI) {
    // Independent mode: each nucleon gets its own GENIE record and sees the
    // full undepleted A-3 remnant, regardless of what happens to the others.
    for (const auto &in : fsi_inputs) {
      std::vector<FSIInputNucleon> one(1, in);
      ran |= fsi::ApplyGenieFSI(A_res, Z_res, one, particles, fRand, fFSIModel);
    }
  } else {
    // Shared mode: all three nucleons in ONE GENIE record so absorption of
    // one leg depletes the remnant seen by the others. GENIE's
    // TransportHadrons iterates the record in insertion order, so a fixed
    // order would systematically give the first nucleon the fullest medium;
    // all three leave the same vertex simultaneously and no order is
    // physical, so the transport order is always randomized.
    for (size_t i = fsi_inputs.size(); i > 1; --i) {
      const size_t j = static_cast<size_t>(fRand->Integer(static_cast<UInt_t>(i)));
      std::swap(fsi_inputs[i - 1], fsi_inputs[j]);
    }
    ran = fsi::ApplyGenieFSI(A_res, Z_res, fsi_inputs, particles,
                             fRand, fFSIModel);
  }

  // Replace the seeded pre-FSI nucleons with GENIE's final state. If the
  // cascade could not run at all, keep the generated nucleons.
  if (ran) fParticles = particles;

  // No post-cascade Pauli veto. Pauli blocking is GENIE's own, and it is
  // partial: absorption and pion production throw INukeException on a
  // blocked final state and the fate/kinematics are rerolled, but the NN
  // elastic / charge-exchange block in TwoBodyCollision is commented out
  // (INukeUtils2018.cxx:826-832), so NN rescattering is not Pauli-blocked
  // (#10).
#else
  (void)v1; (void)v2; (void)v3;
#endif
}
