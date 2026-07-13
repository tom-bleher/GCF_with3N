#ifndef __GENIE_FSI_HELPERS_H__
#define __GENIE_FSI_HELPERS_H__

#include "TLorentzVector.h"
#include "TRandom3.h"
#include <vector>

// FSI model selection for GENIE intranuclear cascade.
//   kHN2018: full hadron-nucleon cascade (HNIntranuke2018) - more accurate
//   kHA2018: effective hadron-nucleus cascade (HAIntranuke2018) - faster
enum FSIModel { kHN2018, kHA2018 };

// Stable secondary produced during a single ApplyGenieFSIToNucleon call.
// The selected outgoing nucleon used to update lead/recoil state is NOT
// included in the secondaries list.
struct FSISecondary {
  int parentRole;     // caller-defined tag (e.g. 0=lead, 1=recoil, 2=third)
  int pdg;
  TLorentzVector p4;
  int rescatterCode;  // GENIE GHepParticle::RescatterCode()
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
// choice of Wright et al. -- REQUIRED for paper comparisons).
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

// Apply GENIE intranuclear cascade to a single nucleon.
// Inputs:
//   A_residual, Z_residual : transport medium (nucleus the nucleon traverses)
//   nucleon_type           : input PDG code (proton/neutron); updated in place
//                            to the post-FSI nucleon PDG (may flip on charge
//                            exchange).
//   p4                     : input 4-momentum; updated in place to post-FSI.
//   x4                     : initial 4-position inside the nucleus.
//   rnd                    : ROOT RNG used to seed GENIE for this call.
//   parentRole             : tag attached to all returned secondaries.
//   secondaries            : appended to with stable descendants other than
//                            the selected outgoing nucleon.
//   model                  : kHN2018 or kHA2018.
// Returns true on success (a stable nucleon was produced); false if the
// nucleon was absorbed or no nucleon descendant could be selected.
bool ApplyGenieFSIToNucleon(int A_residual, int Z_residual,
                            int &nucleon_type,
                            TLorentzVector &p4,
                            const TLorentzVector &x4,
                            TRandom3 *rnd,
                            int parentRole,
                            std::vector<FSISecondary> &secondaries,
                            FSIModel model);

// Per-nucleon input for the multi-hadron FSI helper.
struct FSIInputNucleon {
  int pdg;                  // PDG code on entry (proton/neutron)
  TLorentzVector p4;        // 4-momentum on entry
  TLorentzVector x4;        // initial 4-position inside the nucleus
  int parentRole;           // tag attached to secondaries from this nucleon
};

// Per-nucleon output from the multi-hadron FSI helper.
struct FSIOutputNucleon {
  bool survived;            // false if absorbed / no stable nucleon descendant
  int pdg;                  // post-FSI PDG (may flip on charge exchange)
  TLorentzVector p4;        // post-FSI 4-momentum
};

// Apply GENIE intranuclear cascade to SEVERAL nucleons SIMULTANEOUSLY.
//
// All input nucleons are placed into a single GHepRecord as
// kIStHadronInTheNucleus particles. GENIE's Intranuke2018::TransportHadrons
// then iterates over them sharing a single remnant nucleus (fRemnA/fRemnZ),
// so any absorption / 2p2h removal during the first nucleon's cascade
// properly depletes the medium seen by subsequent nucleons. This is the
// physically correct way to FSI a knocked-out SRC pair or triplet through
// GENIE Intranuke (still sequential within Intranuke's cascade, but with
// shared bookkeeping — no longer independent passes through an undepleted
// A-2 / A-3 remnant).
//
//   A_residual, Z_residual : transport medium AT ENTRY (before any cascade
//                            depletion).
//   inputs                 : one entry per outgoing nucleon.
//   outputs                : resized to inputs.size(); per-nucleon survival
//                            flag and post-FSI (pdg, p4).
//   secondaries            : appended to with stable descendants other than
//                            the selected outgoing nucleon, tagged with the
//                            originating input's parentRole.
//   model                  : kHN2018 or kHA2018.
// Returns true iff EVERY input produced a surviving nucleon. The caller
// can also inspect outputs[i].survived for per-nucleon results.
bool ApplyGenieFSIToNucleons(int A_residual, int Z_residual,
                              const std::vector<FSIInputNucleon> &inputs,
                              std::vector<FSIOutputNucleon> &outputs,
                              TRandom3 *rnd,
                              std::vector<FSISecondary> &secondaries,
                              FSIModel model);

} // namespace fsi

#endif
