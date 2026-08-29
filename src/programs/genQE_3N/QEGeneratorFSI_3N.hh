#ifndef __QE_GENERATOR_FSI_3N_H__
#define __QE_GENERATOR_FSI_3N_H__

#include "QEGenerator_3N.hh"
#include "fsi/GenieFSIHelpers.hh"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <vector>

// QEGeneratorFSI_3N — 3N quasi-elastic generator with GENIE FSI.
//
// Wraps QEGenerator_3N. After the parent samples the (N1, N2, N3) final state
// it transports the three outgoing nucleons through GENIE's intranuclear
// cascade (hN or hA mode), using a single shared SRC-triplet position sampled
// from rho^3(r). Pauli blocking is handled inside GENIE's cascade (blocked
// collisions are rerolled); no additional post-FSI momentum cut is applied.
//
// The FSI output is purely a particle list (GetFinalStateParticles()): GENIE's
// stable final state verbatim, with no lead/recoil designation and no
// absorption special-casing. The generated (pre-FSI) nucleons are returned
// through the generate_event_with_FSI out-parameters unchanged. Which
// final-state particle is "lead" or "recoil" is a detection-level choice made
// downstream from particle types and momenta alone.
//
// The transport medium is the residual nucleus (A_res = A - 3,
// Z_res = Z - n_protons_in_triplet).

class QEGeneratorFSI_3N : public QEGenerator_3N
{
public:
  QEGeneratorFSI_3N(double E, eNCrossSection *thisCS, int thisU,
                    TRandom3 *thisRand,
                    int A = 12, int Z = 6);
  ~QEGeneratorFSI_3N();

  // Generate one event. Identical signature to QEGenerator_3N::generate_event.
  // The out-parameters hold the GENERATED (pre-FSI) nucleons; the post-FSI
  // final state is available via GetFinalStateParticles(). With FSI disabled
  // (or skipped for an unphysical residual) the particle list is simply the
  // three generated nucleons.
  void generate_event_with_FSI(double &weight,
                               int &N1_type, int &N2_type, int &N3_type,
                               TLorentzVector &vk_target,
                               TLorentzVector &vLead_target,
                               TLorentzVector &v2_target,
                               TLorentzVector &v3_target,
                               TLorentzVector &vAm3_target,
                               bool use_CM);

  void EnableFSI(bool e = true)              { fDoFSI = e; }
  void SetFSIModel(FSIModel m)               { fFSIModel = m; }

  // Independent-FSI mode: transport each nucleon through its own GENIE
  // record with an undepleted A-3 remnant, instead of the default shared
  // depleting remnant in a single record.
  void SetIndependentFSI(bool i = true)      { fIndependentFSI = i; }
  bool IndependentFSI() const                { return fIndependentFSI; }

  // Override of QEGenerator_3N::SetTargetNucleus: updates both the
  // kinematic mA / mAm3 (via the parent) and the FSI residual nucleus
  // (A, Z used by ApplyFSI for GENIE transport).
  void SetTargetNucleus(int A, int Z);

  // Post-FSI final state of the last generated event: a flat particle list
  // (parentRole = 0/1/2 for descendants of N1/N2/N3, diagnostic only).
  const std::vector<FSIParticle> & GetFinalStateParticles() const { return fParticles; }

  bool   FSIEnabled() const { return fDoFSI; }
  FSIModel CurrentFSIModel() const { return fFSIModel; }

private:
  void ApplyFSI(int N1_type, int N2_type, int N3_type,
                const TLorentzVector &v1, const TLorentzVector &v2,
                const TLorentzVector &v3);

  TRandom3  *fRand;          // mirrors parent's myRand for use here
  bool       fDoFSI;
  bool       fIndependentFSI = false;
  FSIModel   fFSIModel;
  int        fA;
  int        fZ;

  std::vector<FSIParticle> fParticles;
};

#endif
