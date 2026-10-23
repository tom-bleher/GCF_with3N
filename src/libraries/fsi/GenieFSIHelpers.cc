#include "GenieFSIHelpers.hh"
#include "constants.hh"
#include <cmath>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TMath.h>
#include <TSystem.h>

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

#ifdef USE_FSI
namespace {
// The entries of a GXMLPATH as GENIE reads them (utils::xml::GetXMLFilePath):
// split on ':', ';' and ',', and $VARS and ~ expanded. Without the expansion
// an entry such as a literal '$GENIE/config' would not be recognised and
// would stay ahead of the project override, where GENIE expands it.
std::vector<std::string> SplitPathList(const std::string &pathList)
{
  std::vector<std::string> dirs;
  std::string dir;
  for (size_t i = 0; i <= pathList.size(); ++i) {
    if (i < pathList.size() && std::string(":;,").find(pathList[i]) == std::string::npos) {
      dir += pathList[i];
      continue;
    }
    if (!dir.empty()) {
      char *e = gSystem->ExpandPathName(dir.c_str());
      dirs.push_back(e ? e : dir);
      delete[] e;
    }
    dir.clear();
  }
  return dirs;
}

// Canonical spelling of a directory for comparisons (the spelling itself is
// kept on the path); the input unchanged if it does not exist.
std::string CanonicalDir(const std::string &dir)
{
  char buf[PATH_MAX];
  if (realpath(dir.c_str(), buf)) return buf;
  return dir;
}

// The NuclearModel of the "Default" parameter set of an Intranuke XML file,
// ignoring commented-out lines; empty if not found.
std::string IntranukeNuclearModel(const std::string &file)
{
  std::ifstream in(file.c_str());
  if (!in) return "";
  std::stringstream buf;
  buf << in.rdbuf();
  std::string s = buf.str(), text;
  for (size_t pos = 0; pos < s.size();) {
    const size_t c = s.find("<!--", pos);
    text += s.substr(pos, c == std::string::npos ? std::string::npos : c - pos);
    if (c == std::string::npos) break;
    const size_t e = s.find("-->", c + 4);
    if (e == std::string::npos) break;
    pos = e + 3;
  }
  const size_t set = text.find("<param_set name=\"Default\">");
  if (set == std::string::npos) return "";
  const size_t setEnd = text.find("</param_set>", set);
  const size_t key = text.find("name=\"NuclearModel\"", set);
  if (key == std::string::npos || key > setEnd) return "";
  const size_t a = text.find('>', key), b = text.find("</param>", key);
  if (a == std::string::npos || b == std::string::npos || b < a) return "";
  std::string v = text.substr(a + 1, b - a - 1);
  const size_t f = v.find_first_not_of(" \t\n"), l = v.find_last_not_of(" \t\n");
  return f == std::string::npos ? "" : v.substr(f, l - f + 1);
}
} // namespace
#endif

std::string GenieOverrideDir()
{
  const char *env = std::getenv("GENIE_OVERRIDE_DIR");
  if (env && *env) return env;
#ifdef GENIE_OVERRIDE_DIR
  return GENIE_OVERRIDE_DIR;
#else
  return "";
#endif
}

std::string ResolveGenieXMLPath()
{
#ifdef USE_FSI
  const char *genie_home = std::getenv("GENIE");
  const char *gxml_env = std::getenv("GXMLPATH");
  if (!genie_home) return gxml_env ? gxml_env : "";

  const std::string config_dir = std::string(genie_home) + "/config";
  const std::string default_tune_dir = config_dir + "/G18_10a";
  const std::string override_dir = GenieOverrideDir();
  // The compile-time default override. When $GENIE_OVERRIDE_DIR selects a
  // different directory, an entry naming the default (build/gcf_env.sh and
  // the gcf front end export one) is dropped rather than kept ahead of the
  // selection, so the variable always takes effect.
  std::string default_override;
#ifdef GENIE_OVERRIDE_DIR
  default_override = GENIE_OVERRIDE_DIR;
#endif

  // A selected override without the Intranuke files would silently hand
  // those files to GENIE's own config further down the path.
  if (!override_dir.empty()) {
    std::ifstream probe((override_dir + "/HNIntranuke2018.xml").c_str());
    static bool warned = false;
    if (!probe && !warned) {
      std::cerr << "Warning: GENIE override directory " << override_dir
                << " has no HNIntranuke2018.xml (GENIE_OVERRIDE_DIR?); GENIE "
                   "reads its own defaults for the files it lacks\n";
      warned = true;
    }
  }

  // GENIE reads each XML file from the first GXMLPATH entry that has it.
  // Order: (1) the user's own entries, in their order (e.g. a gcf_fsi
  // --channel directory, or a variant such as genie_override_srctail),
  // (2) the selected project override, (3) the G18_10a tune, (4) GENIE's
  // config. Entries naming (2)-(4) are recognised by canonical path and
  // moved to their slot, so any spelling of them works and a user
  // directory is never shadowed by the project override. An entry naming
  // the default override when another one is selected is dropped.
  std::vector<std::string> fixed;
  std::vector<std::string> seen;
  auto add = [&](const std::string &dir) {
    const std::string canon = CanonicalDir(dir);
    for (const std::string &s : seen) if (s == canon) return;
    seen.push_back(canon);
    fixed.push_back(dir);
  };
  std::vector<std::string> tail;
  if (!override_dir.empty()) tail.push_back(override_dir);
  tail.push_back(default_tune_dir);
  tail.push_back(config_dir);
  std::vector<std::string> tailCanon;
  for (const std::string &d : tail) tailCanon.push_back(CanonicalDir(d));
  if (!default_override.empty()) tailCanon.push_back(CanonicalDir(default_override));
  if (gxml_env)
    for (const std::string &dir : SplitPathList(gxml_env)) {
      const std::string canon = CanonicalDir(dir);
      bool isTail = false;
      for (const std::string &t : tailCanon) if (t == canon) isTail = true;
      if (!isTail) add(dir);
    }
  for (const std::string &d : tail) add(d);

  std::string path;
  for (const std::string &d : fixed) path += (path.empty() ? "" : ":") + d;
  setenv("GXMLPATH", path.c_str(), 1);
  return path;
#else
  const char *gxml_env = std::getenv("GXMLPATH");
  return gxml_env ? gxml_env : "";
#endif
}

std::string FindOnGenieXMLPath(const std::string &file)
{
  const char *gxml_env = std::getenv("GXMLPATH");
  if (!gxml_env) return "";
  std::stringstream ss(gxml_env);
  std::string dir;
  while (std::getline(ss, dir, ':')) {
    if (dir.empty()) continue;
    std::ifstream f((dir + "/" + file).c_str());
    if (f.good()) return dir;
  }
  return "";
}

void PrintGenieXMLConfig(std::ostream &os)
{
#ifdef USE_FSI
  const char *gxml_env = std::getenv("GXMLPATH");
  os << "GENIE XML path: " << (gxml_env ? gxml_env : "(unset)") << "\n";
  for (const char *base : {"HNIntranuke2018.xml", "HAIntranuke2018.xml", "Messenger.xml"}) {
    const std::string dir = FindOnGenieXMLPath(base);
    os << "  " << base << " from " << (dir.empty() ? "(not found)" : dir);
    if (!dir.empty() && std::string(base) != "Messenger.xml") {
      const std::string m = IntranukeNuclearModel(dir + "/" + base);
      if (!m.empty()) os << " (NuclearModel " << m << ")";
    }
    os << "\n";
  }
#else
  (void)os;
#endif
}

namespace {
// SRC-pair vertex radial weight rho^n(r) r^2. Default n=2 (pair density);
// n=1 reproduces GENIE's one-body vertex (the Wright et al. paper setup,
// arXiv:2104.05090, doi:10.48550/arXiv.2104.05090).
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
    int depth = 0;               // 0 = this particle IS the transported input
    int cur = idx;
    while (cur >= 0) {
      if (cur >= first_hadron_slot && cur <= last_hadron_slot) {
        role = slot_role[cur - first_hadron_slot];
        // Intranuke stamps the fate on the FirstMother of each product
        // (HAIntranuke2018.cxx:152, HNIntranuke2018.cxx:155), i.e. on the
        // INCIDENT particle -- never on the stable final state itself.
        // Reading sp->RescatterCode() therefore always returned -1. Take it
        // from the transported input slot this particle descends from:
        // -1 there means that nucleon never interacted.
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
    // 0 = this particle IS the transported input nucleon, unchanged by the
    // cascade; >0 = it is a descendant, i.e. an FSI product.
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
