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
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>

using namespace std;

namespace {

void Usage()
{
  cerr <<
    "Usage: gcf_generate --model 3N|2N|MF --beam E --events N --output FILE [options]\n"
    "\n"
    "Generates N accepted (weight > 0) GCF quasi-elastic events in PWIA (no FSI).\n"
    "Add FSI afterwards with gcf_fsi.\n"
    "\n"
    "Common options:\n"
    "  -m, --model 3N|2N|MF     SRC triplet, SRC pair, or single mean-field nucleon\n"
    "  -b, --beam E             beam energy (GeV)\n"
    "  -n, --events N           number of accepted events to write\n"
    "  -o, --output FILE        output ROOT file\n"
    "  -A, --target-A A         target mass number (default 12)\n"
    "  -Z, --target-Z Z         target charge (default 6)\n"
    "  -S, --seed S             random seed, positive (default: drawn and printed)\n"
    "  -c, --cs MODEL           eN cross section: onshell, cc1 (default), cc2\n"
    "      --sigma-cm S         3N/2N: CM momentum width per component (GeV/c)\n"
    "      --q2 MIN:MAX         Q^2 window (GeV^2); 3N default open, 2N/MF default 1:10\n"
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
    "2N / MF only:\n"
    "      --nn MODEL           2N: NN interaction for the pair contacts (default AV18)\n"
    "      --estar E            2N: residual excitation energy E* (GeV); MF leaves the\n"
    "                           A-1 residual in its ground state\n"
    "      --kcut K             2N: relative-momentum cutoff (GeV/c, default 0.25)\n"
    "      --onset fixed|local  2N: SRC onset k_rel > kcut (default) or > k_F(r)\n"
    "      --prel-is            2N: importance-sample |p_rel| (lower variance)\n"
    "      --fermi-gas global|local  MF: Fermi-gas model (default global)\n"
    "      --kF-fg K            MF: global Fermi-gas momentum (GeV/c, default 0.25)\n"
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

} // namespace

int main(int argc, char **argv)
{
  string model, output, csName = "cc1", nnName = "AV18", onset = "fixed",
         fermiGas = "global", phaseSpace;
  double beamE = 0., sigCM = 0., q2min = 0., q2max = 0., thetaMin = 5.,
         thetaMax = 40., kF = 0., estar = 0., kcut = 0., kFfg = 0.;
  long long nEvents = 0;
  unsigned long seed = 0;
  int A = 12, Z = 6, u = 1;
  bool sigCMSet = false, q2Set = false, kFSet = false, estarSet = false,
       kcutSet = false, kFfgSet = false, useCM = true, doRad = false,
       doCoul = false, prelIS = false, verbose = false, checkOnly = false,
       nnSet = false, onsetSet = false;

  enum { kSigmaCM = 1000, kQ2, kThetaE, kKF, kNoCM, kNN, kEstar, kKcut,
         kOnset, kPrelIS, kFermiGas, kKFfg, kRad, kCoulomb, kPhaseSpace, kCheck };
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
    {nullptr, 0, nullptr, 0}};

  int c;
  while ((c = getopt_long(argc, argv, "m:b:n:o:A:Z:S:c:u:vh", longOpts, nullptr)) != -1) {
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
      case kNN: nnName = optarg; nnSet = true; break;
      case kEstar: estar = atof(optarg); estarSet = true; break;
      case kKcut: kcut = atof(optarg); kcutSet = true; break;
      case kOnset: onset = optarg; onsetSet = true; break;
      case kPrelIS: prelIS = true; break;
      case kFermiGas: fermiGas = optarg; break;
      case kKFfg: kFfg = atof(optarg); kFfgSet = true; break;
      case kRad: doRad = true; break;
      case kCoulomb: doCoul = true; break;
      case kPhaseSpace: phaseSpace = optarg; break;
      case kCheck: checkOnly = true; break;
      default: Usage(); return 1;   // getopt already printed the reason
    }
  }

  // ---- Validate ----------------------------------------------------------
  bool ok = true;
  auto fail = [&ok](const string &msg) { cerr << "Error: " << msg << "\n"; ok = false; };
  if (model != "3N" && model != "2N" && model != "MF") fail("--model must be 3N, 2N or MF");
  if (beamE <= 0.) fail("--beam must be a positive energy in GeV");
  if (nEvents <= 0 && !checkOnly) fail("--events must be positive");
  if (output.empty() && !checkOnly) fail("--output is required");
  if (A < 1 || Z < 0 || Z > A) fail("invalid target (A, Z)");
  csMethod csMeth = cc1;
  if (csName == "onshell") csMeth = onshell;
  else if (csName == "cc1" || csName == "1") csMeth = cc1;
  else if (csName == "cc2" || csName == "2") csMeth = cc2;
  else fail("--cs must be onshell, cc1 or cc2");
  if (model == "3N" && (u < 1 || u > 10)) fail("--wave-function must be one of AV18, AV18+UIX, AV18+TM, CDBonn, CDBonn+TM (or with -pD), or 1-10");
  if (kFSet && (kF <= 0. || kF > 0.5)) fail("--kF expects k_F in GeV/c, e.g. 0.22");
  if (onset != "fixed" && onset != "local") fail("--onset must be fixed or local");
  if (fermiGas != "global" && fermiGas != "local") fail("--fermi-gas must be global or local");
  if (model == "3N" && (!nnName.empty() && nnName != "AV18"))
    fail("--nn applies to 2N; the 3N wave function is chosen with --wave-function");
  // The mean-field model has no pair: no contacts (--nn), no pair CM motion
  // (--sigma-cm), no relative momentum (--kcut, --onset, --prel-is), and the
  // A-1 residual is left in its ground state (--estar).
  if (model == "MF") {
    if (nnSet) fail("--nn applies to 2N, not MF (the mean-field weight has no NN contacts)");
    if (sigCMSet) fail("--sigma-cm applies to 3N/2N, not MF (no SRC-cluster CM motion)");
    if (estarSet) fail("--estar applies to 2N, not MF (the A-1 residual is left in its ground state)");
    if (kcutSet) fail("--kcut applies to 2N, not MF (no SRC pair, so no relative-momentum cutoff)");
    if (onsetSet) fail("--onset applies to 2N, not MF (the SRC onset is a property of the pair)");
    if (prelIS) fail("--prel-is applies to 2N, not MF (no |p_rel| to importance-sample)");
  }
  if (!ok) { cerr << "(gcf_generate --help for usage)\n"; return 1; }
  if (!output.empty()) output = AbsolutePath(output);

  // Seed 0 draws a fresh one; print it either way so the run can be replayed.
  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  if (!checkOnly) cout << "Random seed: " << seed << "\n";
  TRandom3 rand(seed);
  eNCrossSection cs(csMeth, kelly);

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
    if (sigCMSet) nucleus->set_sigmaCM(sigCM);
    if (estarSet) nucleus->set_Estar(estar);
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

  if (checkOnly) {
    // A setup can be valid and still never give a nonzero weight (an empty
    // phase space): try events until the first one is accepted.
    const long long kCheckAttempts = 1000000LL;
    long long tried = 0;
    double w = 0.;
    while (w <= 0. && tried < kCheckAttempts) {
      int t1, t2, t3;
      TLorentzVector vk, v1, v2, v3, vRes;
      if (gen3) gen3->generate_event(w, t1, t2, t3, vk, v1, v2, v3, vRes, useCM);
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
      gen3->generate_event(weight, t1, t2, t3, vk, v1, v2, v3, vRes, useCM);
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
        gen2->generate_event(weight, leadType, recType, vk, vLead, vRec, vRes);
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
  prov += "model: " + model + "\n";
  TNamed("gcf_generate", prov.c_str()).Write();
  outfile.Close();

  cout << "Wrote " << filled << " events (" << attempts << " attempts, efficiency "
       << (attempts ? 100. * filled / attempts : 0.) << "%) to " << output << "\n"
       << "Sum of weights / attempts = " << (attempts ? sumW / attempts : 0.)
       << " (the per-attempt normalization)\n";
  return filled == nEvents ? 0 : 2;
}
