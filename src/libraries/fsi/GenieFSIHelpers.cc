#include "GenieFSIHelpers.hh"
#include "constants.hh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <TMath.h>

#ifdef USE_FSI
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"
#include "Physics/HadronTransport/HNIntranuke2018.h"
#include "Physics/HadronTransport/HAIntranuke2018.h"
#include "Physics/NuclearState/NuclearUtils.h"

// NOTE on the rescatterCode diagnostic stored on each output particle: it is
// the raw fate integer stamped by {HA,HN}Intranuke2018 per
// INukeHadroFates2018.h. Do NOT compare it against the genie::kIHAFt*/kIHNFt*
// enum constants in analysis code compiled against GENIE headers: GENIE ships
// two headers (INukeHadroFates.h and INukeHadroFates2018.h) defining the same
// enum names under the SAME include guard but with DIFFERENT values, so which
// integers the names evaluate to depends on include order. Interpret the raw
// integers against INukeHadroFates2018.h directly.
#endif

namespace {
#ifdef USE_FSI
bool PathListContainsFile(const std::string &pathList, const std::string &filename)
{
  std::stringstream ss(pathList);
  std::string dir;
  while (std::getline(ss, dir, ':')) {
    if (dir.empty()) continue;
    std::ifstream f((dir + "/" + filename).c_str());
    if (f.good()) return true;
  }
  return false;
}

// Rejection-sample r from the radial weight w(r) = rho^n(r) * r^2 inside the
// nucleus; return a 4-position with that |r| and an isotropic direction.
// n=2 -> 2N (pair) sampling, n=3 -> 3N (triplet) sampling.
TLorentzVector SampleRadialPosition(int A, int n, TRandom3 *rnd,
                                    int &cached_A_io,
                                    double &cached_Rmax_io,
                                    double &cached_fmax_io)
{
  if (A != cached_A_io) {
    cached_A_io = A;
    cached_Rmax_io = 2.5 * 1.4 * TMath::Power(A, 1./3.); // fm
    cached_fmax_io = 0.;
    for (int i = 1; i <= 200; i++) {
      double r   = cached_Rmax_io * i / 200.;
      double rho = genie::utils::nuclear::Density(r, A);
      double f   = TMath::Power(rho, n) * r * r;
      if (f > cached_fmax_io) cached_fmax_io = f;
    }
    cached_fmax_io *= 1.01; // small safety margin
  }

  double r;
  do {
    r = cached_Rmax_io * rnd->Rndm();
    double rho = genie::utils::nuclear::Density(r, A);
    double f = TMath::Power(rho, n) * r * r;
    if (rnd->Rndm() * cached_fmax_io <= f) break;
  } while (true);

  double cosTheta = 2. * rnd->Rndm() - 1.;
  double sinTheta = TMath::Sqrt(1. - cosTheta * cosTheta);
  double phi = 2. * TMath::Pi() * rnd->Rndm();
  return TLorentzVector(r * sinTheta * TMath::Cos(phi),
                        r * sinTheta * TMath::Sin(phi),
                        r * cosTheta,
                        0.);
}
#endif // USE_FSI
} // namespace

namespace fsi {

void ResolveGenieXMLPath()
{
#ifdef USE_FSI
  const char *genie_home = std::getenv("GENIE");
  const char *gxml_env = std::getenv("GXMLPATH");
  if (!genie_home) return;

  const std::string config_dir = std::string(genie_home) + "/config";
  const std::string default_tune_dir = config_dir + "/G18_10a";

  // Project-side override directory, defined at compile time by the
  // CMakeLists of the program that links us. Empty if not provided.
  const std::string override_dir =
#ifdef GENIE_OVERRIDE_DIR
    GENIE_OVERRIDE_DIR;
#else
    "";
#endif

  bool has_model = false;
  bool has_tune_list = false;
  bool has_override = override_dir.empty();
  if (gxml_env && std::string(gxml_env).size() > 0) {
    const std::string gxml = gxml_env;
    has_model     = PathListContainsFile(gxml, "ModelConfiguration.xml");
    has_tune_list = PathListContainsFile(gxml, "TuneGeneratorList.xml");
    if (!override_dir.empty()) {
      // crude "is this path already on GXMLPATH" check via the override's
      // signature file
      has_override = PathListContainsFile(gxml, "HNIntranuke2018.xml")
                  && (gxml.find(override_dir) != std::string::npos);
    }
  }

  if (!has_model || !has_tune_list || !has_override) {
    // Build a path that PREPENDS the project override (so it wins over the
    // upstream GENIE configs) and keeps any pre-existing entries the user set.
    std::string fixed_gxml;
    if (!override_dir.empty()) fixed_gxml = override_dir;
    if (!fixed_gxml.empty()) fixed_gxml += ":";
    fixed_gxml += default_tune_dir + ":" + config_dir;
    if (gxml_env && std::string(gxml_env).size() > 0) {
      fixed_gxml += ":";
      fixed_gxml += gxml_env;
    }
    setenv("GXMLPATH", fixed_gxml.c_str(), 1);
  }
#endif
}

namespace {
// SRC-pair vertex radial weight rho^n(r) r^2. Default n=2 (pair density);
// n=1 reproduces GENIE's one-body vertex (the Wright et al. paper setup).
int gSRCVertexRhoPower = 0;   // 0 = unresolved: env or default on first use
int resolveSRCVertexRhoPower()
{
  if (gSRCVertexRhoPower >= 1) return gSRCVertexRhoPower;
  const char *e = std::getenv("FSI_SRC_VERTEX_RHO_POWER");
  int n = e ? std::atoi(e) : 2;
  if (n < 1 || n > 3) n = 2;
  if (e) fprintf(stderr, "SampleSRCPosition: rho power = %d (env override)\n", n);
  gSRCVertexRhoPower = n;
  return n;
}
}

void SetSRCVertexRhoPower(int n)
{
  if (n >= 1 && n <= 3) gSRCVertexRhoPower = n;
}

int GetSRCVertexRhoPower()
{
  return resolveSRCVertexRhoPower();
}

TLorentzVector SampleSRCPosition(int A, TRandom3 *rnd)
{
#ifdef USE_FSI
  const int rho_n = resolveSRCVertexRhoPower();
  static int    cached_n    = -1;
  static int    cached_A    = -1;
  static double cached_Rmax = 0.;
  static double cached_fmax = 0.;
  if (rho_n != cached_n) { cached_n = rho_n; cached_A = -1; }
  return SampleRadialPosition(A, rho_n, rnd, cached_A, cached_Rmax, cached_fmax);
#else
  (void)A; (void)rnd;
  return TLorentzVector(0., 0., 0., 0.);
#endif
}

TLorentzVector SampleSRCPosition3N(int A, TRandom3 *rnd)
{
#ifdef USE_FSI
  static int    cached_A    = -1;
  static double cached_Rmax = 0.;
  static double cached_fmax = 0.;
  return SampleRadialPosition(A, 3, rnd, cached_A, cached_Rmax, cached_fmax);
#else
  (void)A; (void)rnd;
  return TLorentzVector(0., 0., 0., 0.);
#endif
}

TLorentzVector SampleMFPosition(int A, TRandom3 *rnd)
{
#ifdef USE_FSI
  static int    cached_A    = -1;
  static double cached_Rmax = 0.;
  static double cached_fmax = 0.;
  return SampleRadialPosition(A, 1, rnd, cached_A, cached_Rmax, cached_fmax);
#else
  (void)A; (void)rnd;
  return TLorentzVector(0., 0., 0., 0.);
#endif
}

double LocalFermiMomentum(int A, int N_q, double r_fm)
{
#ifdef USE_FSI
  if (A <= 0 || N_q <= 0) return 0.;
  const double rho = genie::utils::nuclear::Density(r_fm, A);  // fm^-3, ∫ = 1
  if (rho <= 0.) return 0.;
  const double kF_fm = std::pow(3.0 * M_PI * M_PI * N_q * rho, 1.0/3.0);
  const double hbarc = 0.1973269804;  // GeV * fm
  return hbarc * kF_fm;               // GeV/c
#else
  (void)A; (void)N_q; (void)r_fm;
  return 0.;
#endif
}

double NuclearDensity(int A, double r_fm)
{
#ifdef USE_FSI
  return genie::utils::nuclear::Density(r_fm, A);
#else
  (void)A; (void)r_fm;
  return -1.;
#endif
}

bool ApplyGenieFSI(int A_residual, int Z_residual,
                   const std::vector<FSIInputNucleon> &inputs,
                   std::vector<FSIParticle> &particles,
                   TRandom3 *rnd,
                   FSIModel model)
{
#ifdef USE_FSI
  if (inputs.empty()) return false;

  // Only proton/neutron inputs can be transported; anything else is passed
  // through to the output list unchanged.
  size_t n_valid = 0;
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (inputs[i].pdg == pCode || inputs[i].pdg == nCode) {
      n_valid++;
    } else {
      FSIParticle part;
      part.pdg = inputs[i].pdg;
      part.p4 = inputs[i].p4;
      part.rescatterCode = -1;
      part.generation = 0;
      part.parentRole = inputs[i].parentRole;
      particles.push_back(part);
    }
  }
  if (n_valid == 0) return false;

  // 1..2^31-1: SetSeed(0) would mean "random UUID seed" and break replay.
  genie::RandomGen::Instance()->SetSeed(1 + rnd->Integer(0x7ffffffe));

  static genie::HNIntranuke2018 hN2018("genie::HNIntranuke2018");
  static genie::HAIntranuke2018 hA2018("genie::HAIntranuke2018");
  static bool initialized = false;
  if (!initialized) {
    hN2018.Configure("Default");
    hA2018.Configure("Default");
    initialized = true;
  }

  genie::GHepRecord evrec;
  const int pdg_tgt = genie::pdg::IonPdgCode(A_residual, Z_residual);
  // Physical residual-ion mass from GENIE's PDG library. Fallback to A*mN
  // if the ion is unknown (the ~8 MeV/nucleon binding otherwise leaks into
  // the remnant four-momentum handed to the cascade).
  double m_tgt = mN * static_cast<double>(A_residual);
  if (TParticlePDG *ion = genie::PDGLibrary::Instance()->Find(pdg_tgt)) {
    if (ion->Mass() > 0.) m_tgt = ion->Mass();
  }

  // Build event record for lepton-nucleus mode so GENIE skips
  // GenerateVertex() and trusts the pre-set SRC vertex on each hadron.
  //
  // Slot 0: fake electron probe  -> triggers kGMdLeptonNucleus.
  // Slot 1: target nucleus       -> used by SetTrackingRadius.
  // Slot 2: remnant nucleus      -> used by TransportHadrons for the shared
  //                                  fRemnA/fRemnZ/fRemnP4 bookkeeping. As
  //                                  GENIE absorbs a nucleon during one
  //                                  hadron's cascade the remnant depletes,
  //                                  so subsequent hadrons in the same record
  //                                  traverse the properly reduced medium.
  // Slots 3..3+N-1: the N outgoing nucleons in kIStHadronInTheNucleus.

  evrec.AddParticle(11, genie::kIStInitialState,
                    -1, -1, -1, -1,
                    0., 0., 100., 100.,
                    0., 0., 0., 0.);

  evrec.AddParticle(pdg_tgt, genie::kIStInitialState,
                    -1, -1, 2, 2,
                    0., 0., 0., m_tgt,
                    0., 0., 0., 0.);

  evrec.AddParticle(pdg_tgt, genie::kIStStableFinalState,
                    1, -1, -1, -1,
                    0., 0., 0., m_tgt,
                    0., 0., 0., 0.);

  // Slot -> parentRole map for the transported inputs. Slots are recorded at
  // insertion time (invalid inputs add no particle, so slots are contiguous
  // over the valid ones only).
  const int first_hadron_slot = evrec.GetEntries();
  std::vector<int> slot_role;
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (!(inputs[i].pdg == pCode || inputs[i].pdg == nCode)) continue;
    slot_role.push_back(inputs[i].parentRole);
    evrec.AddParticle(inputs[i].pdg, genie::kIStHadronInTheNucleus,
                      0, -1, -1, -1,
                      inputs[i].p4.Px(), inputs[i].p4.Py(),
                      inputs[i].p4.Pz(), inputs[i].p4.E(),
                      inputs[i].x4.X(), inputs[i].x4.Y(),
                      inputs[i].x4.Z(), inputs[i].x4.T());
  }
  const int last_hadron_slot = first_hadron_slot
                             + static_cast<int>(slot_role.size()) - 1;

  if (model == kHN2018) {
    hN2018.ProcessEventRecord(&evrec);
  } else {
    hA2018.ProcessEventRecord(&evrec);
  }

  // Collect GENIE's final state verbatim: every stable final-state particle
  // except the remnant nucleus / ions and GENIE pseudo-particles (e.g. the
  // hadronic blob hA uses for 4-momentum bookkeeping). Provenance is found
  // by walking mother links back to a transported input slot.
  for (int idx = 0; idx < evrec.GetEntries(); ++idx) {
    genie::GHepParticle* sp = evrec.Particle(idx);
    if (!sp || !sp->P4()) continue;
    if (sp->Status() != genie::kIStStableFinalState) continue;
    const int pdg = sp->Pdg();
    if (genie::pdg::IsIon(pdg) || genie::pdg::IsPseudoParticle(pdg)) continue;

    int role = -1;
    int fate = -1;
    int depth = 0;               // steps up the mother chain to an input slot
    int cur = idx;
    while (cur >= 0) {
      if (cur >= first_hadron_slot && cur <= last_hadron_slot) {
        role = slot_role[cur - first_hadron_slot];
        // Intranuke stamps the fate on the FirstMother of each product
        // (HAIntranuke2018.cxx:152, HNIntranuke2018.cxx:155), i.e. on the
        // INCIDENT particle -- never on the stable final state itself.
        // Reading sp->RescatterCode() therefore always returned -1. Take it
        // from the transported input slot this particle descends from.
        // A nucleon that left without interacting has fate 1 (no
        // interaction, set in Intranuke2018::TransportHadrons); -1 means the
        // input was never stepped.
        genie::GHepParticle* anc = evrec.Particle(cur);
        if (anc) fate = anc->RescatterCode();
        break;
      }
      genie::GHepParticle* cp = evrec.Particle(cur);
      if (!cp) break;
      cur = cp->FirstMother();
      ++depth;
    }

    FSIParticle part;
    part.pdg = pdg;
    part.p4 = *(sp->P4());
    part.rescatterCode = fate;
    // Mother-chain depth from the input slot. Intranuke rescatters a clone
    // whose FirstMother is the input slot, so every cascade output, an
    // unscattered nucleon included, has depth >= 1. The products of a
    // collision take the input slot as mother too (INukeUtils2018.cxx:848,
    // HNIntranuke2018.cxx:656): in hA they are final at depth 1, so a
    // rescattered or knocked-out nucleon is also at depth 1; in hN they are
    // stepped again and end up deeper. Use rescatterCode == 1, not
    // generation, to find nucleons that left unchanged; 0 occurs only for
    // the pass-through branch above.
    part.generation = depth;
    part.parentRole = role;
    particles.push_back(part);
  }

  return true;
#else
  (void)A_residual; (void)Z_residual; (void)inputs; (void)particles;
  (void)rnd; (void)model;
  return false;
#endif
}

} // namespace fsi
