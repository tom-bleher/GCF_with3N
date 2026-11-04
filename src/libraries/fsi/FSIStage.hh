#ifndef __FSI_STAGE_H__
#define __FSI_STAGE_H__

#include "fsi/GenieFSIHelpers.hh"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <vector>

// FSI as a stand-alone step: "knocked-out nucleons in, final state out",
// independent of which generator produced the event. This is the step that
// QEGeneratorFSI (2N) and QEGeneratorFSI_3N do inline after generating an
// event; gcf_fsi applies it to every event of an event file.
//
// Conventions, identical to the generators:
//   - the transport medium is the residual nucleus: the target minus the
//     knocked-out protons and neutrons;
//   - the cascade is skipped (event passed through unchanged) when the
//     residual has A < 3, as GENIE INTRANUKE is not validated there;
//   - all nucleons start from one common vertex inside the FULL target;
//   - shared mode puts them in one GENIE record (absorption of one depletes
//     the remnant seen by the others) in a random transport order;
//     independent mode gives each its own record and an undepleted remnant.

namespace fsi {

enum TransportMode { kSharedRemnant, kIndependent };

struct StageInput {
  int targetA = 0;
  int targetZ = 0;
  // Particles leaving the primary vertex. Protons and neutrons are
  // transported; anything else is passed through to the final state as is.
  // x4 is ignored: all of them start from `vertex`.
  std::vector<FSIInputNucleon> particles;
  bool haveVertex = false;
  TLorentzVector vertex;     // fm, target rest frame
};

struct StageResult {
  bool ran = false;          // the GENIE cascade ran for this event
  bool vertexSampled = false;  // vertex was drawn here (input had none)
  TLorentzVector vertex;     // vertex used (fm); zero if none was needed
  int residualA = 0;         // transport medium
  int residualZ = 0;
  // Final state. When the cascade did not run this is the input particles
  // unchanged (rescatterCode -1, generation 0).
  std::vector<FSIParticle> particles;
};

// Radial weight of a sampled vertex, rho^n(r) r^2. kVertexAuto picks n from
// the number of transported nucleons: 1 -> rho (mean-field nucleon),
// 2 -> the configured SRC-pair power (SetSRCVertexRhoPower, default 2),
// >= 3 -> rho^3 (triplet).
const int kVertexAuto = 0;

void TransportEvent(const StageInput &in, TransportMode mode, FSIModel model,
                    int vertexPower, TRandom3 *rnd, StageResult &out);

// Coulomb correction of the particles leaving the nucleus, the rule of the
// generators' -C / --coulomb: a particle of charge q (in units of e) gains
// q * deltaE (gcfGenerator::coulombCorrection, fixed mass). A negative
// particle with kinetic energy T <= |q| deltaE cannot escape and is dropped;
// the return value is the number dropped. Neutral and unknown PDG codes are
// untouched. deltaE is Delta E_C of the nucleus (gcfGenerator::calcCoulombEnergy).
int CoulombCorrectFinalState(std::vector<FSIParticle> &particles, double deltaE);

// The inverse shift on the particles entering the cascade: q * deltaE is
// taken off again, so a list stored with the correction applied (gcf_generate
// --coulomb) is transported with its momenta inside the nucleus.
void CoulombUncorrect(std::vector<FSIInputNucleon> &particles, double deltaE);

} // namespace fsi

#endif
