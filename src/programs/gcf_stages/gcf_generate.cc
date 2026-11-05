// gcf_generate -- stage 1 of the GCF pipeline: quasi-elastic (e,e'N...) events
// from the Generalized Contact Formalism in PWIA, i.e. WITHOUT final-state
// interactions, written in the GCF event-file format (docs/PIPELINE.md).
// gcf_fsi adds FSI to such a file as a separate step.
//
// Models:
//   3N  struck nucleon from an SRC triplet (QEGenerator_3N)
//   2N  struck nucleon from an SRC pair (QEGeneratorFSI, FSI off)
//   MF  single mean-field nucleon from a Fermi gas (QEGeneratorFSI, FSI off)

#include "QEGenerator_3N.hh"
#include "QEGeneratorFSI.hh"
#include "constants.hh"
#include "fsi/GenieFSIHelpers.hh"
#include "io/GCFEventFile.hh"
#include "io/Provenance.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TTree.h"

#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

using namespace std;

namespace {

void Usage()
{
  cerr <<
    "Usage: gcf_generate --model 3N|2N|MF --beam E --events N --output FILE [options]\n"
    "\n"
    "Generates N accepted (weight > 0) GCF quasi-elastic events in PWIA (no FSI).\n"
    "Add FSI afterwards with gcf_fsi. An option of another model is an error.\n"
    "\n"
    "Common options:\n"
    "  -m, --model 3N|2N|MF     SRC triplet, SRC pair, or single mean-field nucleon\n"
    "  -b, --beam E             beam energy (GeV)\n"
    "  -n, --events N           number of accepted events to write\n"
    "  -o, --output FILE        output ROOT file\n"
    "  -A, --target-A A         target mass number (default 12)\n"
    "  -Z, --target-Z Z         target charge (default 6)\n"
    "  -S, --seed S             random seed, positive (default: drawn and printed)\n"
    "      --ff MODEL           nucleon form factors: kelly (default), ye, dipole\n"
    "  -c, --cs MODEL           eN cross section: onshell, cc1 (default), cc2\n"
    "      --q2 MIN:MAX         Q^2 window (GeV^2); 3N default open, 2N/MF default 1:10\n"
    "      --print-params       print the nuclear inputs of this setup as one line,\n"
    "                           PARAMS {json}, then exit without generating\n"
    "      --check              validate the options and set the generator up\n"
    "                           (nucleus tables, ...) and try up to 10^6 events: fails\n"
    "                           if none has a nonzero weight. Writes nothing\n"
    "  -v, --verbose            progress output\n"
    "  -h, --help               this message\n"
    "\n"
    "3N only:\n"
    "  -u, --wave-function W    3N wave function: AV18, AV18+UIX, AV18+TM, CDBonn,\n"
    "                           CDBonn+TM (3He Faddeev, default AV18), or the same\n"
    "                           with -pD for the deuteron channel only (AV18-pD, ...);\n"
    "                           the table numbers 1-10 are accepted too\n"
    "      --theta-e MIN:MAX    electron polar angle window (deg, default 5:40)\n"
    "      --kF K               internal-momentum cutoff k > kF (GeV/c, default 0.25);\n"
    "                           also sets sigma-cm = sqrt(3/5) kF unless given\n"
    "      --no-cm              no CM motion of the triplet\n"
    "\n"
    "3N and 2N:\n"
    "      --sigma-cm S         CM momentum width of the SRC cluster per component\n"
    "                           (GeV/c; default 0.15 for 3N, the nucleus table for 2N)\n"
    "      --estar E            excitation energy E* of the residual (GeV, default 0);\n"
    "                           not for a target with no residual\n"
    "\n"
    "2N only:\n"
    "      --nn MODEL           NN interaction for the pair contacts (default AV18)\n"
    "      --sigma-e S          draw E* per event from a Gaussian of this width around\n"
    "                           --estar, keeping E* > 0 (GeV)\n"
    "      --cpp0 C, --cnn0 C, --cpn0 C, --cpn1 C\n"
    "                           override the pair contacts of the nucleus table, in the\n"
    "                           table's units (--print-params shows them); only change\n"
    "                           the weights. Not for the deuteron\n"
    "      --randomize S        vary the nuclear inputs once, within their uncertainties,\n"
    "                           with a random stream of its own seeded by S: sigma_CM and\n"
    "                           the contacts (Gaussian), E* (uniform up to the table's\n"
    "                           maximum). The options above then override the result.\n"
    "                           Give every job of a run the same S\n"
    "      --light-cone         use the light-cone cross section; excludes --prel-is and\n"
    "                           --onset local\n"
    "      --kcut K             relative-momentum cutoff (GeV/c, default 0.25)\n"
    "      --onset fixed|local  SRC onset k_rel > kcut (default) or > k_F(r)\n"
    "      --prel-is            importance-sample |p_rel| (lower variance)\n"
    "\n"
    "MF only (the A-1 residual is left in its ground state):\n"
    "      --fermi-gas global|local  Fermi-gas model (default global)\n"
    "      --kF-fg K            global Fermi-gas momentum (GeV/c, default 0.25)\n"
    "\n"
    "2N and MF:\n"
    "      --rad                peaking-approximation radiation\n"
    "      --coulomb            Coulomb correction (applied before any FSI stage)\n"
    "      --phase-space FILE   phase-space limits file (gcfGenerator format)\n";
}

bool ParseRange(const char *arg, double &lo, double &hi)
{
  const char *sep = strchr(arg, ':');
  if (!sep) return false;
  lo = atof(arg);
  hi = atof(sep + 1);
  return hi > lo;
}

string AbsolutePath(const string &path)
{
  if (!path.empty() && path[0] == '/') return path;
  char cwd[PATH_MAX];
  if (!getcwd(cwd, sizeof cwd)) return path;
  return string(cwd) + "/" + path;
}

// The models each model-specific option applies to; an option outside its
// models is an error, not a silent no-op. The options not listed (model, beam,
// events, output, target, seed, cs, q2, check, verbose) apply to every model.
// tools/gcf/schema.py carries the same table (Opt.models), and
// test_model_options.py checks that the binary and the schema agree.
struct ModelOption {
  const char *flag;
  const char *models;   // '/'-separated
  const char *what;     // why it does not apply elsewhere
};
const ModelOption kModelOptions[] = {
  {"--sigma-cm", "3N/2N", "CM motion of the SRC cluster; MF has no cluster"},
  {"--wave-function", "3N", "three-body density of the triplet; 2N uses --nn"},
  {"--theta-e", "3N", "electron-angle window of the 3N generator; the 2N/MF generator has none"},
  {"--kF", "3N", "internal-momentum cutoff of the triplet; 2N uses --kcut, MF --kF-fg"},
  {"--no-cm", "3N", "CM motion of the triplet; 2N sets its width with --sigma-cm"},
  {"--nn", "2N", "NN interaction of the pair contacts; 3N uses --wave-function, MF has no contacts"},
  {"--estar", "3N/2N", "excitation of the cluster's residual; MF leaves the A-1 residual in its ground state"},
  {"--sigma-e", "2N", "per-event E* distribution of the pair's A-2 residual"},
  {"--cpp0", "2N", "pair contact; 3N has no contacts yet (#73), MF has no pair"},
  {"--cnn0", "2N", "pair contact; 3N has no contacts yet (#73), MF has no pair"},
  {"--cpn0", "2N", "pair contact; 3N has no contacts yet (#73), MF has no pair"},
  {"--cpn1", "2N", "pair contact; 3N has no contacts yet (#73), MF has no pair"},
  {"--randomize", "2N", "the table of pair inputs and their uncertainties exists for 2N only"},
  {"--light-cone", "2N", "the light-cone cross section exists for the pair generator only"},
  {"--kcut", "2N", "relative-momentum cutoff of the SRC pair"},
  {"--onset", "2N", "SRC onset of the pair's relative momentum"},
  {"--prel-is", "2N", "importance sampling of the pair's |p_rel|"},
  {"--fermi-gas", "MF", "Fermi-gas model of the mean-field nucleon"},
  {"--kF-fg", "MF", "Fermi momentum of the mean-field Fermi gas"},
  {"--rad", "2N/MF", "radiation exists in the 2N/MF generator only"},
  {"--coulomb", "2N/MF", "Coulomb correction exists in the 2N/MF generator only"},
  {"--phase-space", "2N/MF", "limits file of the 2N/MF generator; 3N uses --q2 and --theta-e"},
};

bool AppliesTo(const char *models, const string &model)
{
  const string list = string("/") + models + "/";
  return list.find("/" + model + "/") != string::npos;
}

} // namespace

int main(int argc, char **argv)
{
  string model, output, csName = "cc1", ffName = "kelly", nnName = "AV18", onset = "fixed",
         fermiGas = "global", phaseSpace;
  double beamE = 0., sigCM = 0., q2min = 0., q2max = 0., thetaMin = 5.,
         thetaMax = 40., kF = 0., estar = 0., kcut = 0., kFfg = 0., sigmaE = 0.,
         cpp0 = -1., cnn0 = -1., cpn0 = -1., cpn1 = -1.;
  long long nEvents = 0;
  unsigned long seed = 0, randomizeSeed = 0;
  int A = 12, Z = 6, u = 1;
  bool lightCone = false, printParams = false, sigmaESet = false, randomize = false;
  bool sigCMSet = false, q2Set = false, kFSet = false, estarSet = false,
       kcutSet = false, kFfgSet = false, useCM = true, doRad = false,
       doCoul = false, prelIS = false, verbose = false, checkOnly = false;
  bool contactGiven = false, badContact = false;
  vector<string> given;   // model-specific options on the command line

  enum { kSigmaCM = 1000, kQ2, kThetaE, kKF, kNoCM, kNN, kEstar, kKcut,
         kOnset, kPrelIS, kFermiGas, kKFfg, kRad, kCoulomb, kPhaseSpace, kCheck, kFF,
         kSigmaE, kCpp0, kCnn0, kCpn0, kCpn1, kRandomize, kLightCone, kPrintParams };
  static const option longOpts[] = {
    {"model", required_argument, nullptr, 'm'},
    {"beam", required_argument, nullptr, 'b'},
    {"events", required_argument, nullptr, 'n'},
    {"output", required_argument, nullptr, 'o'},
    {"target-A", required_argument, nullptr, 'A'},
    {"target-Z", required_argument, nullptr, 'Z'},
    {"seed", required_argument, nullptr, 'S'},
    {"cs", required_argument, nullptr, 'c'},
    {"wave-function", required_argument, nullptr, 'u'},
    {"interaction", required_argument, nullptr, 'u'},   // old name
    {"verbose", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {"sigma-cm", required_argument, nullptr, kSigmaCM},
    {"q2", required_argument, nullptr, kQ2},
    {"theta-e", required_argument, nullptr, kThetaE},
    {"kF", required_argument, nullptr, kKF},
    {"no-cm", no_argument, nullptr, kNoCM},
    {"nn", required_argument, nullptr, kNN},
    {"estar", required_argument, nullptr, kEstar},
    {"kcut", required_argument, nullptr, kKcut},
    {"onset", required_argument, nullptr, kOnset},
    {"prel-is", no_argument, nullptr, kPrelIS},
    {"fermi-gas", required_argument, nullptr, kFermiGas},
    {"kF-fg", required_argument, nullptr, kKFfg},
    {"rad", no_argument, nullptr, kRad},
    {"coulomb", no_argument, nullptr, kCoulomb},
    {"phase-space", required_argument, nullptr, kPhaseSpace},
    {"check", no_argument, nullptr, kCheck},
    {"ff", required_argument, nullptr, kFF},
    {"sigma-e", required_argument, nullptr, kSigmaE},
    {"cpp0", required_argument, nullptr, kCpp0},
    {"cnn0", required_argument, nullptr, kCnn0},
    {"cpn0", required_argument, nullptr, kCpn0},
    {"cpn1", required_argument, nullptr, kCpn1},
    {"randomize", required_argument, nullptr, kRandomize},
    {"light-cone", no_argument, nullptr, kLightCone},
    {"print-params", no_argument, nullptr, kPrintParams},
    {nullptr, 0, nullptr, 0}};

  int c;
  int optIndex = -1;
  while ((c = getopt_long(argc, argv, "m:b:n:o:A:Z:S:c:u:vh", longOpts, &optIndex)) != -1) {
    string flag;
    if (c == 'u') flag = "--wave-function";   // also -u and --interaction
    else if (c >= kSigmaCM && optIndex >= 0) flag = string("--") + longOpts[optIndex].name;
    if (!flag.empty() && find(given.begin(), given.end(), flag) == given.end())
      given.push_back(flag);
    optIndex = -1;
    switch (c) {
      case 'm': model = optarg; break;
      case 'b': beamE = atof(optarg); break;
      case 'n': nEvents = atoll(optarg); break;
      case 'o': output = optarg; break;
      case 'A': A = atoi(optarg); break;
      case 'Z': Z = atoi(optarg); break;
      case 'S': seed = strtoul(optarg, nullptr, 10); break;
      case 'c': csName = optarg; break;
      case 'u': {
        static const char *names[] = {"AV18", "AV18+UIX", "AV18+TM", "CDBonn", "CDBonn+TM",
                                      "AV18-pD", "AV18+UIX-pD", "AV18+TM-pD", "CDBonn-pD", "CDBonn+TM-pD"};
        u = 0;
        for (int i = 0; i < 10; ++i) if (strcmp(optarg, names[i]) == 0) u = i + 1;
        if (u == 0) u = atoi(optarg);
        break;
      }
      case 'v': verbose = true; break;
      case 'h': Usage(); return 0;
      case kSigmaCM: sigCM = atof(optarg); sigCMSet = true; break;
      case kQ2:
        if (!ParseRange(optarg, q2min, q2max)) {
          cerr << "Invalid --q2 (expect MIN:MAX with MAX > MIN).\n"; return 1;
        }
        q2Set = true; break;
      case kThetaE:
        if (!ParseRange(optarg, thetaMin, thetaMax)) {
          cerr << "Invalid --theta-e (expect MIN:MAX in deg with MAX > MIN).\n"; return 1;
        }
        break;
      case kKF: kF = atof(optarg); kFSet = true; break;
      case kNoCM: useCM = false; break;
      case kNN: nnName = optarg; break;
      case kEstar: estar = atof(optarg); estarSet = true; break;
      case kKcut: kcut = atof(optarg); kcutSet = true; break;
      case kOnset: onset = optarg; break;
      case kPrelIS: prelIS = true; break;
      case kFermiGas: fermiGas = optarg; break;
      case kKFfg: kFfg = atof(optarg); kFfgSet = true; break;
      case kRad: doRad = true; break;
      case kCoulomb: doCoul = true; break;
      case kPhaseSpace: phaseSpace = optarg; break;
      case kCheck: checkOnly = true; break;
      case kFF: ffName = optarg; break;
      case kSigmaE: sigmaE = atof(optarg); sigmaESet = true; break;
      case kCpp0: cpp0 = atof(optarg); contactGiven = true; if (cpp0 < 0.) badContact = true; break;
      case kCnn0: cnn0 = atof(optarg); contactGiven = true; if (cnn0 < 0.) badContact = true; break;
      case kCpn0: cpn0 = atof(optarg); contactGiven = true; if (cpn0 < 0.) badContact = true; break;
      case kCpn1: cpn1 = atof(optarg); contactGiven = true; if (cpn1 < 0.) badContact = true; break;
      case kRandomize: randomizeSeed = strtoul(optarg, nullptr, 10); randomize = true; break;
      case kLightCone: lightCone = true; break;
      case kPrintParams: printParams = true; break;
      default: Usage(); return 1;   // getopt already printed the reason
    }
  }

  // ---- Validate ----------------------------------------------------------
  bool ok = true;
  auto fail = [&ok](const string &msg) { cerr << "Error: " << msg << "\n"; ok = false; };
  if (model != "3N" && model != "2N" && model != "MF") fail("--model must be 3N, 2N or MF");
  if (beamE <= 0.) fail("--beam must be a positive energy in GeV");
  const bool noOutput = checkOnly || printParams;
  if (nEvents <= 0 && !noOutput) fail("--events must be positive");
  if (output.empty() && !noOutput) fail("--output is required");
  if (A < 1 || Z < 0 || Z > A) fail("invalid target (A, Z)");
  csMethod csMeth = cc1;
  if (csName == "onshell") csMeth = onshell;
  else if (csName == "cc1" || csName == "1") csMeth = cc1;
  else if (csName == "cc2" || csName == "2") csMeth = cc2;
  else fail("--cs must be onshell, cc1 or cc2");
  ffModel ffMod = kelly;
  if (!parse_ffModel(ffName.c_str(), ffMod)) fail("--ff must be kelly, ye or dipole");
  if (model == "3N" && (u < 1 || u > 10)) fail("--wave-function must be one of AV18, AV18+UIX, AV18+TM, CDBonn, CDBonn+TM (or with -pD), or 1-10");
  if (kFSet && (kF <= 0. || kF > 0.5)) fail("--kF expects k_F in GeV/c, e.g. 0.22");
  if (onset != "fixed" && onset != "local") fail("--onset must be fixed or local");
  if (fermiGas != "global" && fermiGas != "local") fail("--fermi-gas must be global or local");
  if (model == "3N" || model == "2N" || model == "MF") {
    for (const string &flag : given) {
      for (const ModelOption &o : kModelOptions) {
        if (flag != o.flag || AppliesTo(o.models, model)) continue;
        fail(flag + " applies to " + o.models + ", not " + model + " (" + o.what + ")");
        break;
      }
    }
  }
  // A negative contact would silently act like 0: only events with a positive
  // weight are kept.
  if (badContact) fail("the contacts --cpp0, --cnn0, --cpn0 and --cpn1 must be >= 0");
  if (estarSet && estar < 0.) fail("--estar must be >= 0");
  if (sigmaESet && sigmaE <= 0.) fail("--sigma-e must be positive");
  if (randomize && randomizeSeed == 0) fail("--randomize needs a positive seed");
  if (lightCone && prelIS) fail("--light-cone excludes --prel-is: it does not sample |p_rel| that way");
  if (lightCone && onset == "local") fail("--light-cone excludes --onset local: it has no local onset");
  if (model == "2N" && A == 2) {
    for (const char *flag : {"--cpp0", "--cnn0", "--cpn0", "--cpn1", "--randomize"})
      if (find(given.begin(), given.end(), flag) != given.end())
        fail(string(flag) + " does not apply to the deuteron: it uses a momentum distribution, not contacts");
  }
  if (model == "3N" && A == 3 && estarSet) fail("--estar does not apply to A = 3: a triplet leaves no residual");
  if (!ok) { cerr << "(gcf_generate --help for usage)\n"; return 1; }
  if (!output.empty()) output = AbsolutePath(output);

  // Seed 0 draws a fresh one; print it either way so the run can be replayed.
  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  if (!noOutput) cout << "Random seed: " << seed << "\n";
  TRandom3 rand(seed);
  eNCrossSection cs(csMeth, ffMod);

  // ---- Generator setup (before opening the output: a setup error must not
  // leave a truncated file behind) ----------------------------------------
  unique_ptr<QEGenerator_3N> gen3;
  unique_ptr<gcfNucleus> nucleus;
  unique_ptr<QEGeneratorFSI> gen2;
  if (model == "3N") {
    // QEGenerator_3N reads its density tables from array/ relative to the
    // working directory; construct it from the genQE_3N source directory.
    char cwd[PATH_MAX];
    const bool haveCwd = getcwd(cwd, sizeof cwd) != nullptr;
    if (chdir(GCF_3N_SOURCE_DIR) != 0)
      cerr << "Warning: cannot enter " << GCF_3N_SOURCE_DIR
           << "; looking for the 3N tables in the current directory\n";
    gen3.reset(new QEGenerator_3N(beamE, &cs, u, &rand));
    if (haveCwd && chdir(cwd) != 0) { cerr << "Error: cannot return to " << cwd << "\n"; return 1; }
    gen3->SetTargetNucleus(A, Z);
    if (kFSet) {
      gen3->SetInternalMomentumMin(kF);
      if (!sigCMSet) gen3->SetSigCM(sqrt(3. / 5.) * kF);
    }
    if (sigCMSet) gen3->SetSigCM(sigCM);
    gen3->set_theta_k_maxmin(thetaMin, thetaMax);
    if (q2Set) gen3->SetQ2Range(q2min, q2max);
  } else {
    const bool deuteron = (A == 2 && Z == 1);
    string uType = (deuteron && nnName == "AV18") ? "AV18_deut" : nnName;
    nucleus.reset(new gcfNucleus(Z, A - Z, const_cast<char *>(uType.c_str())));
    // One draw per run: its own stream, so the events of the jobs still differ
    // while all of them use the same nuclear inputs.
    if (randomize) {
      TRandom3 randomizer(randomizeSeed);
      nucleus->randomize(&randomizer);
    }
    if (sigCMSet) nucleus->set_sigmaCM(sigCM);
    if (estarSet) nucleus->set_Estar(estar);
    if (sigmaESet) nucleus->set_sigmaE(sigmaE);
    gcfSRC *src = nucleus->get_SRC();
    if (cpp0 >= 0.) src->set_Cpp0(cpp0);
    if (cnn0 >= 0.) src->set_Cnn0(cnn0);
    if (cpn0 >= 0.) src->set_Cpn0(cpn0);
    if (cpn1 >= 0.) src->set_Cpn1(cpn1);
    if (src->get_Cpp0() + src->get_Cnn0() + src->get_Cpn0() + src->get_Cpn1() <= 0.) {
      cerr << "Error: all contacts are zero, so every event has zero weight\n"; return 1;
    }
    gen2.reset(new QEGeneratorFSI(beamE, nucleus.get(), &cs, &rand));
    gen2->EnableFSI(false);
    if (deuteron) gen2->set_deuteron();
    if (!phaseSpace.empty() &&
        !gen2->parse_phase_space_file(const_cast<char *>(phaseSpace.c_str()))) {
      cerr << "Error: could not read phase-space file " << phaseSpace << "\n"; return 1;
    }
    if (kcutSet) gen2->set_pRel_cut(kcut);
    if (q2Set) gen2->set_QSq_range(q2min, q2max);
    if (doRad) gen2->set_doRad(true);
    if (doCoul) gen2->set_doCoul(true);
    if (prelIS) gen2->set_pRel_importance_sampling(true);
    if (onset == "local") gen2->SetSrcOnset(QEGeneratorFSI::kLocalKF);
    gen2->SetFermiGasMode(fermiGas == "local" ? QEGeneratorFSI::kLocalFG
                                              : QEGeneratorFSI::kGlobalFG);
    if (kFfgSet) gen2->SetFermiMomentumFG(kFfg);
    // Without the A-1 masses every MF event gets weight 0 and the run would
    // spin to the attempt limit.
    if (model == "MF" && (nucleus->get_mAm1p() <= 0. || nucleus->get_mAm1n() <= 0.)) {
      cerr << "Error: --model MF needs the A-1 residual masses, which the nucleus tables "
              "have only for 12C (not A=" << A << ", Z=" << Z << ")\n";
      return 1;
    }
  }

  // What the nuclear inputs of this run are, for the provenance and the front end.
  ostringstream inputs;
  inputs << setprecision(8);
  if (gen3) {
    inputs << "sigma_cm=" << (useCM ? gen3->GetSigCM() : 0.) << " estar=" << estar;
  } else if (nucleus) {
    gcfSRC *src = nucleus->get_SRC();
    inputs << "sigma_cm=" << nucleus->get_sigmaCM() << " estar=" << nucleus->get_Estar()
           << " sigma_e=" << nucleus->get_sigmaE();
    if (model == "2N")
      inputs << " cpp0=" << src->get_Cpp0() << " cnn0=" << src->get_Cnn0()
             << " cpn0=" << src->get_Cpn0() << " cpn1=" << src->get_Cpn1();
  }

  if (printParams) {
    ostringstream j;
    j << setprecision(8) << "{\"model\": \"" << model << "\", \"A\": " << A << ", \"Z\": " << Z;
    if (gen3) {
      j << ", \"sigma_cm\": [" << gen3->GetSigCM() << ", null], \"estar\": " << estar
        << ", \"kF\": " << gen3->GetInternalMomentumMin()
        << ", \"mA\": " << gen3->GetTargetMass() << ", \"mRes\": {\"residual\": " << gen3->GetResidualMass() << "}";
    } else {
      gcfSRC *src = nucleus->get_SRC();
      j << ", \"sigma_cm\": [" << nucleus->get_sigmaCM() << ", " << nucleus->get_d_sigmaCM() << "]"
        << ", \"estar\": " << nucleus->get_Estar() << ", \"sigma_e\": " << nucleus->get_sigmaE()
        << ", \"estar_max\": " << nucleus->get_Estar_max() << ", \"mA\": " << nucleus->get_mA();
      if (model == "2N") {
        j << ", \"nn\": \"" << nnName << "\", \"contact_source\": \"" << src->get_ContactSource() << "\""
          << ", \"contacts\": {\"cpp0\": [" << src->get_Cpp0() << ", " << src->get_dCpp0() << "]"
          << ", \"cnn0\": [" << src->get_Cnn0() << ", " << src->get_dCnn0() << "]"
          << ", \"cpn0\": [" << src->get_Cpn0() << ", " << src->get_dCpn0() << "]"
          << ", \"cpn1\": [" << src->get_Cpn1() << ", " << src->get_dCpn1() << "]}"
          << ", \"mRes\": {\"pp\": " << nucleus->get_mAmpp() << ", \"pn\": " << nucleus->get_mAmpn()
          << ", \"nn\": " << nucleus->get_mAmnn() << "}";
      } else {
        j << ", \"mRes\": {\"p\": " << nucleus->get_mAm1p() << ", \"n\": " << nucleus->get_mAm1n() << "}";
      }
    }
    j << "}";
    cout << "PARAMS " << j.str() << "\n";
    return 0;
  }

  if (checkOnly) {
    // A setup can be valid and still never give a nonzero weight (an empty
    // phase space): try events until the first one is accepted.
    const long long kCheckAttempts = 1000000LL;
    long long tried = 0;
    double w = 0.;
    while (w <= 0. && tried < kCheckAttempts) {
      int t1, t2, t3;
      TLorentzVector vk, v1, v2, v3, vRes;
      if (gen3) { double e3 = estar; gen3->generate_event(w, t1, t2, t3, vk, v1, v2, v3, vRes, useCM, e3); }
      else if (model == "2N" && lightCone) gen2->generate_event_lightcone(w, t1, t2, vk, v1, v2, vRes);
      else if (model == "2N") gen2->generate_event(w, t1, t2, vk, v1, v2, vRes);
      else gen2->generate_event_MF(w, t1, t2, vk, v1, v2, vRes);
      ++tried;
    }
    if (w <= 0.) {
      cerr << "Error: no event with a nonzero weight in " << tried << " attempts: the "
           << "requested phase space looks empty. Check the beam energy, target and windows.\n";
      return 1;
    }
    cout << "OK: " << model << " on A=" << A << " Z=" << Z << " at " << beamE
         << " GeV is set up (first accepted event after " << tried
         << " attempts); nothing written (--check)\n";
    return 0;
  }

  // ---- Output ------------------------------------------------------------
  TFile outfile(output.c_str(), "RECREATE");
  if (outfile.IsZombie()) { cerr << "Error: cannot create " << output << "\n"; return 1; }
  TTree *tree = new TTree(gcfio::kTreeName, "GCF events (PWIA)");
  gcfio::EventBuffers ev;
  ev.Book(tree);
  ev.targetA = A; ev.targetZ = Z; ev.beamE = beamE;

  // ---- Event loop: until nEvents events with weight > 0 -------------------
  const long long kMaxAttemptsWithoutSuccess = 100000000LL;
  long long attempts = 0, filled = 0;
  const long long progressEvery = std::max(1LL, nEvents / 100);
  double sumW = 0.;
  while (filled < nEvents) {
    double weight = 0.;
    TLorentzVector vk, vRes;
    ev.ClearParticles();
    TLorentzVector vertex;
    bool haveVertex = false;

    if (gen3) {
      int t1, t2, t3;
      TLorentzVector v1, v2, v3;
      double e3 = estar;   // 0 unless --estar: gives the same events as without it
      gen3->generate_event(weight, t1, t2, t3, vk, v1, v2, v3, vRes, useCM, e3);
      if (weight > 0.) {
        ev.AddParticle(t1, 0, v1);
        ev.AddParticle(t2, 1, v2);
        ev.AddParticle(t3, 2, v3);
        // Triplet vertex, rho^3-weighted, as the 3N FSI generator draws it.
        vertex = fsi::SampleSRCPosition3N(A, &rand);
        haveVertex = true;
      }
    } else {
      int leadType, recType;
      TLorentzVector vLead, vRec;
      if (model == "2N") {
        if (lightCone) gen2->generate_event_lightcone(weight, leadType, recType, vk, vLead, vRec, vRes);
        else gen2->generate_event(weight, leadType, recType, vk, vLead, vRec, vRes);
        if (weight > 0.) haveVertex = gen2->GetLastSrcVertex(vertex);
      } else {
        gen2->generate_event_MF(weight, leadType, recType, vk, vLead, vRec, vRes);
        if (weight > 0.) { vertex = gen2->GetLastMFPosition(); haveVertex = true; }
      }
      if (weight > 0.) {
        ev.AddParticle(leadType, 0, vLead);
        if (recType == pCode || recType == nCode) ev.AddParticle(recType, 1, vRec);
      }
    }
    ++attempts;

    if (weight > 0.) {
      ev.weight = weight;
      ev.SetElectron(vk);
      ev.SetResidual(vRes);
      ev.SetVertex(haveVertex, vertex);
      tree->Fill();
      ++filled;
      sumW += weight;
      if (filled % progressEvery == 0 || filled == nEvents)
        cout << "Progress: " << filled << "/" << nEvents << endl;
    }
    if (verbose && attempts % 1000000 == 0)
      cout << "Attempts " << attempts << ", accepted " << filled << "/" << nEvents << "\n";
    if (filled == 0 && attempts >= kMaxAttemptsWithoutSuccess) {
      cerr << "No accepted events after " << attempts << " attempts: the requested "
           << "phase space looks empty. Check the beam energy, target and windows.\n";
      break;
    }
  }

  // ---- Bookkeeping -------------------------------------------------------
  outfile.cd();
  tree->Write();
  TH1D hAttempts(gcfio::kAttemptsName, "generation attempts", 1, 0., 1.);
  hAttempts.SetBinContent(1, static_cast<double>(attempts));
  hAttempts.Write();
  string prov = gcfio::Provenance("gcf_generate", argc, argv, seed);
  prov += string("form factors: ") + ffModel_name(ffMod) + "\n";
  prov += "model: " + model + "\n";
  prov += "nuclear inputs: " + inputs.str() + "\n";
  TNamed("gcf_generate", prov.c_str()).Write();
  outfile.Close();

  cout << "Wrote " << filled << " events (" << attempts << " attempts, efficiency "
       << (attempts ? 100. * filled / attempts : 0.) << "%) to " << output << "\n"
       << "Sum of weights / attempts = " << (attempts ? sumW / attempts : 0.)
       << " (the per-attempt normalization)\n";
  return filled == nEvents ? 0 : 2;
}
