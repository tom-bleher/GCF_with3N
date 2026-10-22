#ifndef __GENIE_FSI_HELPERS_H__
#define __GENIE_FSI_HELPERS_H__

#include "TLorentzVector.h"
#include "TRandom3.h"
#include <vector>

// FSI model selection for GENIE intranuclear cascade.
//   kHN2018: full hadron-nucleon cascade (HNIntranuke2018) - more accurate
//   kHA2018: effective hadron-nucleus cascade (HAIntranuke2018) - faster
enum FSIModel { kHN2018, kHA2018 };

// One stable final-state particle of the GENIE cascade. The cascade output is
// purely a particle list: no particle is designated lead/recoil here, and
// absorption of an input nucleon is not special-cased (the absorbed nucleon
// simply does not appear; its ejecta do). Any lead/recoil identification is a
// DETECTION-level choice made downstream from types and momenta alone.
struct FSIParticle {
  int pdg;
  TLorentzVector p4;
  // Fate of the INPUT nucleon this particle descends from, as stamped by
  // Intranuke (INukeHadroFates2018.h). Note Intranuke stamps the fate on the
  // incident particle, not on the stable final state, so this must be read
  // from the ancestor -- reading it off the final-state particle always
  // yields -1. A nucleon that left the nucleus without interacting has fate
  // 1 (kIHNFtNoInteraction / kIHAFtNoInteraction) and exactly its input
  // momentum; -1 means the cascade did not run at all.
  int rescatterCode;
  // Mother-chain depth from the transported input nucleon. Intranuke
  // re-emits even a non-interacting nucleon as a daughter, so this is >= 1
  // for every cascade output; 0 only when the cascade did not run. To tell
  // cascade products from nucleons that left unchanged use rescatterCode
  // (== 1 for the latter), not generation. parentRole cannot answer it
  // either: a knocked-out medium nucleon inherits the role of whichever
  // incident nucleon knocked it out.
  int generation;
  int parentRole;     // caller tag of the input nucleon this particle
                      // descends from (provenance diagnostic; -1 if unknown)
};

// Per-nucleon input for the FSI helper.
struct FSIInputNucleon {
  int pdg;                  // PDG code on entry (proton/neutron)
  TLorentzVector p4;        // 4-momentum on entry
  TLorentzVector x4;        // initial 4-position inside the nucleus
  int parentRole;           // provenance tag copied onto descendants
};

namespace fsi {

// Set GXMLPATH so GENIE can locate ModelConfiguration.xml /
// TuneGeneratorList.xml when the user has not sourced GENIE setup scripts.
// Idempotent; safe to call multiple times.
void ResolveGenieXMLPath();

// One-body nucleon number density rho(r) [fm^-3] at radius r [fm] for nucleus
// A, from GENIE's nuclear density model. Returns -1 if FSI/GENIE is not
// compiled in. Used to build the local-Fermi-gas momentum distribution.
double NuclearDensity(int A, double r_fm);

// Sample a position inside the nucleus weighted by rho^n(r) * r^2, where n
// is the configured SRC-vertex rho power (default 2: pair density ~ square
// of the one-body density; n=1 reproduces GENIE's one-body vertex, the
// choice of Wright et al., arXiv:2104.05090, doi:10.48550/arXiv.2104.05090
// -- REQUIRED for paper comparisons).
TLorentzVector SampleSRCPosition(int A, TRandom3 *rnd);

// Configure the SRC-pair vertex radial weight rho^n r^2 (n in [1,3]).
// Also settable via env FSI_SRC_VERTEX_RHO_POWER; the programmatic setter
// wins. Default 2.
void SetSRCVertexRhoPower(int n);
int  GetSRCVertexRhoPower();

// Sample a position inside the nucleus weighted by rho^3(r) * r^2.
// Models the SRC triplet location: triplet density goes as the cube of the
// single-nucleon density.
TLorentzVector SampleSRCPosition3N(int A, TRandom3 *rnd);

// Sample a position inside the nucleus weighted by rho^1(r) * r^2.
// Models a single mean-field nucleon location: it follows the one-body
// nucleon density directly (not its square as for an SRC pair).
TLorentzVector SampleMFPosition(int A, TRandom3 *rnd);

// Local Fermi momentum k_F^q(r) for nucleon species q at radius r, computed
// per GENIE's LocalFGM convention:
//   k_F^q(r) = ( 3 pi^2 * N_q * Density(r, A) )^{1/3}  [fm^-1]
// where N_q is the count of species-q nucleons (Z for protons, A-Z for
// neutrons), Density() is the GENIE one-body density (normalised to 1), and
// the result is converted to GeV/c with hbar c. Returns 0 if FSI is off
// or the density is non-positive.
double LocalFermiMomentum(int A, int N_q, double r_fm);

// Transport the input nucleons through GENIE's intranuclear cascade and
// APPEND every stable final-state particle to `particles`.
//
// All inputs are placed into a single GHepRecord as kIStHadronInTheNucleus
// particles; GENIE's Intranuke2018::TransportHadrons iterates over them
// sharing a single remnant nucleus (fRemnA/fRemnZ), so absorption during one
// nucleon's cascade depletes the medium seen by the others. For independent
// transport, call once per nucleon with a single-element input vector.
//
// The output is GENIE's final state verbatim: no surviving-nucleon selection,
// no daughter-chain walking, no absorption special-casing. Pseudo-particles
// and the remnant nucleus are excluded; everything else with status
// kIStStableFinalState is returned, tagged with the parentRole of the input
// nucleon it descends from.
//
// Returns true if the cascade ran; false if FSI is not compiled in or no
// valid nucleon input was given (in which case `particles` is not modified,
// and the caller should treat the inputs as the final state).
bool ApplyGenieFSI(int A_residual, int Z_residual,
                   const std::vector<FSIInputNucleon> &inputs,
                   std::vector<FSIParticle> &particles,
                   TRandom3 *rnd,
                   FSIModel model);

} // namespace fsi

#endif
