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

// Intranuke fate codes as stamped by {HA,HN}Intranuke2018 on the interacting
// node's RescatterCode (INukeHadroFates2018.h).
// CAUTION: do NOT use the genie::kIHAFt*/kIHNFt* enum constants here. GENIE
// ships two headers (INukeHadroFates.h and INukeHadroFates2018.h) that define
// the same enum names under the SAME include guard (_INTRANUKE_FATES_H_) but
// with DIFFERENT values (the old header keeps kIHAFtElas, shifting Abs 4->5
// and Cmp 6->7). Which one a translation unit sees depends on include order,
// so the named constants silently evaluate to the wrong integers here while
// the Intranuke .cxx files stamp the 2018 values.
static const int kFateHA2018_Abs = 4;   // kIHAFtAbs  in INukeHadroFates2018.h
static const int kFateHA2018_Cmp = 6;   // kIHAFtCmp  in INukeHadroFates2018.h
static const int kFateHN2018_Abs = 5;   // kIHNFtAbs  in INukeHadroFates2018.h
static const int kFateHN2018_Cmp = 6;   // kIHNFtCmp  in INukeHadroFates2018.h

// True if this node's stamped fate means the incident hadron ceased to exist
// (absorption / compound-nucleus pre-equilibrium): its daughters are
// multi-nucleon ejecta, not a continuation of the incident particle.
static bool FateKillsIncident(int fate, FSIModel model)
{
  if (model == kHA2018) return fate == kFateHA2018_Abs || fate == kFateHA2018_Cmp;
  return fate == kFateHN2018_Abs || fate == kFateHN2018_Cmp;
}
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

#ifdef USE_FSI
namespace {

// Pick the surviving nucleon from input slot `transported_idx` and append all
// other stable descendants to `secondaries`. Returns the GHEP index of the
// surviving nucleon (>=0) or -1 if absorbed / no nucleon descendant found.
int SelectSurvivingNucleonAndCollectSecondaries(
    genie::GHepRecord &evrec,
    int transported_idx,
    int parentRole,
    std::vector<FSISecondary> &secondaries,
    FSIModel model)
{
  // Intranuke stamps each interaction's fate on the interacting node's
  // RescatterCode (in hA there is a single interaction, stamped on the input
  // slot; in hN every cascade generation is stamped on the node whose
  // interaction it was -- exactly the nodes the chain-walk below visits).
  // For the absorption/compound fates the incident nucleon is genuinely
  // gone: its "daughters" are multi-nucleon ejecta, NOT a continuation of
  // the incident particle. Without this check the FirstDaughter chain-walk
  // (or the highest-E fallback below) promotes one of those soft ejecta to
  // "the surviving nucleon", which (a) removes all absorption attenuation
  // and (b) fakes n->p conversions whenever the picked ejecta's charge
  // differs -- inflating the (e,e'pp) yield.
  //
  // Walk the cascade chain via FirstDaughter until a stable particle is
  // found, checking each node's fate for leg death along the way.
  int chain_idx = transported_idx;
  while (true) {
    genie::GHepParticle* cp = evrec.Particle(chain_idx);
    if (!cp) break;
    const int fate = cp->RescatterCode();
    if (std::getenv("FSI_DEBUG_FATE"))
      fprintf(stderr, "FATEDBG slot=%d pdg=%d fate=%d\n",
              chain_idx, cp->Pdg(), fate);
    if (FateKillsIncident(fate, model)) {
      // Absorbed at this generation: the leg dies. Everything the input
      // produced (ejecta of this and earlier generations) becomes ordinary
      // secondaries.
      std::vector<int>* ds = evrec.GetStableDescendants(transported_idx);
      if (ds) {
        for (int idx : *ds) {
          genie::GHepParticle* sp = evrec.Particle(idx);
          if (!sp || !sp->P4()) continue;
          if (sp->Status() != genie::kIStStableFinalState) continue;
          FSISecondary sec;
          sec.parentRole = parentRole;
          sec.pdg = sp->Pdg();
          sec.p4 = *(sp->P4());
          sec.rescatterCode = sp->RescatterCode();
          secondaries.push_back(sec);
        }
        delete ds;
      }
      return -1;   // absorbed: no surviving nucleon
    }
    if (cp->Status() == genie::kIStStableFinalState) break;
    int fd = cp->FirstDaughter();
    if (fd < 0) break;
    chain_idx = fd;
  }

  std::vector<int>* daughters = evrec.GetStableDescendants(transported_idx);
  std::vector<int> stable_descendant_indices;
  if (daughters) {
    for (int idx : *daughters) {
      genie::GHepParticle* sp = evrec.Particle(idx);
      if (!sp) continue;
      if (sp->Status() != genie::kIStStableFinalState) continue;
      stable_descendant_indices.push_back(idx);
    }
    delete daughters;
  }

  int best_idx = -1;
  {
    genie::GHepParticle* ep = evrec.Particle(chain_idx);
    if (ep && ep->Status() == genie::kIStStableFinalState &&
        (ep->Pdg() == pCode || ep->Pdg() == nCode)) {
      best_idx = chain_idx;
    }
  }

  // Fallback: highest-energy stable nucleon descendant.
  if (best_idx < 0) {
    double best_E = -1.;
    for (int idx : stable_descendant_indices) {
      genie::GHepParticle* sp = evrec.Particle(idx);
      if (!sp) continue;
      const int pdg = sp->Pdg();
      if (pdg != pCode && pdg != nCode) continue;
      if (sp->E() > best_E) {
        best_E = sp->E();
        best_idx = idx;
      }
    }
  }

  for (int idx : stable_descendant_indices) {
    if (idx == best_idx) continue;
    genie::GHepParticle* sp = evrec.Particle(idx);
    if (!sp || !sp->P4()) continue;

    FSISecondary sec;
    sec.parentRole = parentRole;
    sec.pdg = sp->Pdg();
    sec.p4 = *(sp->P4());
    sec.rescatterCode = sp->RescatterCode();
    secondaries.push_back(sec);
  }

  return best_idx;
}

} // anonymous namespace
#endif // USE_FSI

bool ApplyGenieFSIToNucleon(int A_residual, int Z_residual,
                            int &nucleon_type,
                            TLorentzVector &p4,
                            const TLorentzVector &x4,
                            TRandom3 *rnd,
                            int parentRole,
                            std::vector<FSISecondary> &secondaries,
                            FSIModel model)
{
#ifdef USE_FSI
  // Single-nucleon convenience wrapper around the multi-hadron path: keeps
  // the call site simple for the mean-field generator (no shared remnant
  // concern, only one outgoing nucleon).
  std::vector<FSIInputNucleon> inputs(1);
  inputs[0].pdg = nucleon_type;
  inputs[0].p4 = p4;
  inputs[0].x4 = x4;
  inputs[0].parentRole = parentRole;

  std::vector<FSIOutputNucleon> outputs;
  const bool ok = ApplyGenieFSIToNucleons(A_residual, Z_residual,
                                          inputs, outputs,
                                          rnd, secondaries, model);
  if (outputs.size() != 1) return false;
  if (!outputs[0].survived) return false;
  nucleon_type = outputs[0].pdg;
  p4 = outputs[0].p4;
  return ok;
#else
  (void)A_residual; (void)Z_residual; (void)nucleon_type; (void)p4;
  (void)x4; (void)rnd; (void)parentRole; (void)secondaries; (void)model;
  return false;
#endif
}

bool ApplyGenieFSIToNucleons(int A_residual, int Z_residual,
                              const std::vector<FSIInputNucleon> &inputs,
                              std::vector<FSIOutputNucleon> &outputs,
                              TRandom3 *rnd,
                              std::vector<FSISecondary> &secondaries,
                              FSIModel model)
{
  outputs.assign(inputs.size(), FSIOutputNucleon{false, 0, TLorentzVector()});
#ifdef USE_FSI
  if (inputs.empty()) return false;

  // Sanity-check inputs: only proton/neutron are supported.
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (!(inputs[i].pdg == pCode || inputs[i].pdg == nCode)) {
      // Echo the input through unchanged so the caller can see which slot
      // was invalid, but report failure.
      outputs[i].survived = false;
      outputs[i].pdg = inputs[i].pdg;
      outputs[i].p4 = inputs[i].p4;
    }
  }

  genie::RandomGen::Instance()->SetSeed(rnd->Integer(0x7fffffff));

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

  const int first_hadron_slot = 3;
  std::vector<int> input_slot(inputs.size(), -1);
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (!(inputs[i].pdg == pCode || inputs[i].pdg == nCode)) continue;
    input_slot[i] = first_hadron_slot + static_cast<int>(i);
    evrec.AddParticle(inputs[i].pdg, genie::kIStHadronInTheNucleus,
                      0, -1, -1, -1,
                      inputs[i].p4.Px(), inputs[i].p4.Py(),
                      inputs[i].p4.Pz(), inputs[i].p4.E(),
                      inputs[i].x4.X(), inputs[i].x4.Y(),
                      inputs[i].x4.Z(), inputs[i].x4.T());
  }

  if (model == kHN2018) {
    hN2018.ProcessEventRecord(&evrec);
  } else {
    hA2018.ProcessEventRecord(&evrec);
  }

  bool all_survived = true;
  for (size_t i = 0; i < inputs.size(); ++i) {
    if (input_slot[i] < 0) {
      all_survived = false;
      continue;
    }
    const int best_idx = SelectSurvivingNucleonAndCollectSecondaries(
        evrec, input_slot[i], inputs[i].parentRole, secondaries, model);
    if (best_idx < 0) {
      outputs[i].survived = false;
      all_survived = false;
      continue;
    }
    genie::GHepParticle* best = evrec.Particle(best_idx);
    if (!best || !best->P4()) {
      outputs[i].survived = false;
      all_survived = false;
      continue;
    }
    outputs[i].survived = true;
    outputs[i].pdg = best->Pdg();
    outputs[i].p4 = *(best->P4());
  }

  return all_survived;
#else
  (void)A_residual; (void)Z_residual; (void)inputs;
  (void)rnd; (void)secondaries; (void)model;
  return false;
#endif
}

} // namespace fsi
