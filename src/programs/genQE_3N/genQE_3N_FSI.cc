#include <cstring>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "QEGeneratorFSI_3N.hh"
#include "fsi/GenieFSIHelpers.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
bool useCM   = true;   // -C 0 to disable CM smearing
TRandom3 * myRand;
eNCrossSection * myCS;
QEGeneratorFSI_3N * myGen;
TTree * outtree;

// Tree variables
Double_t pe[3], p1_gen[3], p2_gen[3], p3_gen[3], pAm3[3];
Double_t weight;
Int_t N1_type, N2_type, N3_type;
Int_t doFSI_flag;
Int_t fsiModelInt;       // 0 = hN, 1 = hA
// Post-FSI final state: purely a flat particle list (GENIE's cascade output
// verbatim; with FSI off, the three generated nucleons). No particle is
// designated lead/recoil here -- the analysis builds the measured event
// (detection-based) from types and momenta alone.
Int_t nPart;
std::vector<int>    part_pdg, part_parentRole, part_rescatterCode;
std::vector<double> part_px, part_py, part_pz, part_E;
long long totalAttempts = 0;

void Usage()
{
  cerr << "Usage: ./code <Beam energy (GeV)> [Interaction #] <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-c: Specify eN cross section model (default cc1)\n"
       << "-f: FSI model (hN or hA, default hN)\n"
       << "-A: Target mass number (default 12)\n"
       << "-Z: Target charge number (default 6)\n"
       << "-s: CM Gaussian width sigCM in GeV/c per component (default 0.15)\n"
       << "-C: Toggle CM smearing: 1=on (default), 0=off\n"
       << "-n: Disable FSI (PWIA only)\n"
       << "-I: Independent FSI: transport each nucleon through its own GENIE\n"
       << "    record with an undepleted A-3 remnant (default: shared depleting\n"
       << "    remnant in a single record)\n"
       << "-t: Electron polar-angle range in deg, format MIN:MAX (default 5:40)\n"
       << "-q: Q^2 window in GeV^2, format MIN:MAX (default open). Events outside\n"
       << "    are rejected at generation time (not filled).\n"
       << "-K: Fermi momentum kF in GeV/c. Sets the internal-momentum lower cutoff\n"
       << "    (k > kF on the initial-state wavefunction-frame momenta; default 0.25)\n"
       << "    AND, unless -s is given, sigCM = sqrt(3/5)*kF (3N triplet CM width\n"
       << "    from three uncorrelated mean-field nucleons).\n"
       << "-h: Print this message and exit\n\n\n";
}

bool init(int argc, char ** argv)
{
  int numargs = 5;

  if (argc < numargs)
    {
      Usage();
      return false;
    }

  double Ebeam = atof(argv[1]);
  int u = atoi(argv[2]);
  nEvents = atoi(argv[4]);

  ffModel ffMod = kelly;
  csMethod csMeth = cc1;
  FSIModel fsiModel = kHN2018;
  int A = 12;
  int Z = 6;
  double sigCM_GeV = 0.15;
  bool   sigCM_set = false;
  bool enableFSI = true;
  bool independentFSI = false;
  double thetaKmin_deg = 5.0;
  double thetaKmax_deg = 40.0;
  double Q2min = 0.;
  double Q2max = 1.e9;
  double kF_GeV = 0.;
  bool   kF_set = false;

  int c;
  while ((c = getopt(argc-numargs+1, &argv[numargs-1], "vc:f:A:Z:s:C:nIt:q:K:h")) != -1)
    switch (c) {
      case 'v':
        verbose = true; break;
      case 'c':
        if (strcmp(optarg,"onshell") == 0) csMeth = onshell;
        else if ((strcmp(optarg,"cc1") == 0) || (atoi(optarg) == 1)) csMeth = cc1;
        else if ((strcmp(optarg,"cc2") == 0) || (atoi(optarg) == 2)) csMeth = cc2;
        else { cerr << "Invalid cross section. Allowed: onshell, cc1, cc2.\n"; return false; }
        break;
      case 'f':
        if (strcmp(optarg,"hN") == 0 || strcmp(optarg,"HN") == 0) fsiModel = kHN2018;
        else if (strcmp(optarg,"hA") == 0 || strcmp(optarg,"HA") == 0) fsiModel = kHA2018;
        else { cerr << "Invalid FSI model. Allowed: hN, hA.\n"; return false; }
        break;
      case 'A':
        A = atoi(optarg); break;
      case 'Z':
        Z = atoi(optarg); break;
      case 's':
        sigCM_GeV = atof(optarg); sigCM_set = true; break;
      case 'C':
        useCM = (atoi(optarg) != 0); break;
      case 'n':
        enableFSI = false; break;
      case 'I':
        independentFSI = true; break;
      case 't': {
        // parse MIN:MAX in degrees
        char *sep = strchr(optarg, ':');
        if (!sep) { cerr << "Invalid -t format. Use MIN:MAX (deg).\n"; return false; }
        *sep = '\0';
        thetaKmin_deg = atof(optarg);
        thetaKmax_deg = atof(sep + 1);
        if (thetaKmax_deg <= thetaKmin_deg) {
          cerr << "Invalid -t range: MAX must exceed MIN.\n"; return false;
        }
        break;
      }
      case 'q': {
        char *sep = strchr(optarg, ':');
        if (!sep) { cerr << "Invalid -q format. Use MIN:MAX (GeV^2).\n"; return false; }
        *sep = '\0';
        Q2min = atof(optarg);
        Q2max = atof(sep + 1);
        if (Q2max <= Q2min) {
          cerr << "Invalid -q range: MAX must exceed MIN.\n"; return false;
        }
        break;
      }
      case 'K':
        kF_GeV = atof(optarg);
        kF_set = true;
        if (kF_GeV <= 0. || kF_GeV > 0.5) {
          cerr << "Invalid -K value (expect kF in GeV/c, e.g. 0.220).\n"; return false;
        }
        break;
      case 'h':
        Usage(); return false;
      default:
        // getopt already printed "invalid option"
        Usage(); return false;
    }

  myRand = new TRandom3(0);
  myCS   = new eNCrossSection(csMeth, ffMod);

  myGen = new QEGeneratorFSI_3N(Ebeam, myCS, u, myRand, A, Z);
  myGen->SetTargetNucleus(A, Z);  // updates both kinematic mA/mAm3 and FSI fA/fZ
  if (kF_set) {
    myGen->SetInternalMomentumMin(kF_GeV);
    if (!sigCM_set) myGen->SetSigCM(sqrt(3./5.) * kF_GeV);
  }
  if (sigCM_set) myGen->SetSigCM(sigCM_GeV);
  myGen->set_theta_k_maxmin(thetaKmin_deg, thetaKmax_deg);
  myGen->SetQ2Range(Q2min, Q2max);
  myGen->EnableFSI(enableFSI);
  myGen->SetFSIModel(fsiModel);
  myGen->SetIndependentFSI(independentFSI);

  doFSI_flag  = enableFSI ? 1 : 0;
  fsiModelInt = (fsiModel == kHN2018) ? 0 : 1;

  // Open the output only after the generator is set up: exit() from a
  // setup error with an open TFile crashes in ROOT's teardown and buries
  // the real error message.
  outfile = new TFile(argv[3],"RECREATE");
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("N1_type",  &N1_type,  "N1_type/I");
  outtree->Branch("N2_type",  &N2_type,  "N2_type/I");
  outtree->Branch("N3_type",  &N3_type,  "N3_type/I");
  outtree->Branch("pe",        pe,       "pe[3]/D");
  // Generated (pre-FSI) nucleon momenta: p1_gen is the struck nucleon,
  // p2_gen/p3_gen the spectator pair. Initial-state truth, not the measured
  // final state.
  outtree->Branch("p1_gen",    p1_gen,   "p1_gen[3]/D");
  outtree->Branch("p2_gen",    p2_gen,   "p2_gen[3]/D");
  outtree->Branch("p3_gen",    p3_gen,   "p3_gen[3]/D");
  outtree->Branch("pAm3",      pAm3,     "pAm3[3]/D");
  outtree->Branch("doFSI",    &doFSI_flag,  "doFSI/I");
  outtree->Branch("fsiModel", &fsiModelInt, "fsiModel/I");
  outtree->Branch("weight",   &weight,      "weight/D");
  outtree->Branch("nPart",    &nPart,       "nPart/I");
  outtree->Branch("part_pdg",           &part_pdg);
  outtree->Branch("part_parentRole",    &part_parentRole);
  outtree->Branch("part_rescatterCode", &part_rescatterCode);
  outtree->Branch("part_px",            &part_px);
  outtree->Branch("part_py",            &part_py);
  outtree->Branch("part_pz",            &part_pz);
  outtree->Branch("part_E",             &part_E);

  return true;
}

void evnt(int /*event*/)
{
  TLorentzVector vk, vLead, v2, v3, vAm3;

  myGen->generate_event_with_FSI(weight, N1_type, N2_type, N3_type,
                                 vk, vLead, v2, v3, vAm3, useCM);

  pe[0] = vk.X();    pe[1] = vk.Y();    pe[2] = vk.Z();
  p1_gen[0] = vLead.X(); p1_gen[1] = vLead.Y(); p1_gen[2] = vLead.Z();
  p2_gen[0] = v2.X();    p2_gen[1] = v2.Y();    p2_gen[2] = v2.Z();
  p3_gen[0] = v3.X();    p3_gen[1] = v3.Y();    p3_gen[2] = v3.Z();
  pAm3[0] = vAm3.X(); pAm3[1] = vAm3.Y(); pAm3[2] = vAm3.Z();

  part_pdg.clear(); part_parentRole.clear(); part_rescatterCode.clear();
  part_px.clear();  part_py.clear();  part_pz.clear();  part_E.clear();
  const std::vector<FSIParticle> &parts = myGen->GetFinalStateParticles();
  for (const auto &p : parts) {
    part_pdg.push_back(p.pdg);
    part_parentRole.push_back(p.parentRole);
    part_rescatterCode.push_back(p.rescatterCode);
    part_px.push_back(p.p4.X()); part_py.push_back(p.p4.Y());
    part_pz.push_back(p.p4.Z()); part_E.push_back(p.p4.T());
  }
  nPart = static_cast<int>(part_pdg.size());

  if (weight > 0.) {
    outtree->Fill();
  }
}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
  // Total generation attempts (including weight = 0 rejections, which are
  // not stored in the tree). Kept as a 1-bin histogram so hadd sums it
  // across parallel chunks; needed for per-attempt weight normalization
  // when generation-time cuts (e.g. -q) reject different fractions in
  // different samples.
  TH1D hAttempts("hAttempts", "generation attempts", 1, 0., 1.);
  hAttempts.SetBinContent(1, (double)totalAttempts);
  hAttempts.Write();
  outfile->Delete("genTbuffer;*");
  outfile->Close();
}

int main(int argc, char ** argv)
{
  if (not init(argc, argv)) return -1;

  // Count successful (weight > 0) events. Keep generating until the target
  // number of successful events have been written to the tree.
  const long long kMaxAttemptsWithoutSuccess = 100000000LL;
  long long attempts = 0;
  long long filled   = 0;
  while (filled < (long long)nEvents) {
    if ((attempts % 1000000 == 0) && verbose)
      cout << "Attempts " << attempts << ", filled " << filled
           << "/" << nEvents << "\n";
    evnt((int)attempts);
    if (weight > 0.) filled++;
    attempts++;
    // An impossible -t/-q/target choice has zero acceptance and would
    // otherwise loop forever; give up once it is clear nothing gets through.
    if (filled == 0 && attempts >= kMaxAttemptsWithoutSuccess) {
      cerr << "No accepted events after " << attempts << " attempts: the "
           << "requested phase space looks empty. Check the beam energy and "
           << "-t/-q/-A/-Z settings. Aborting.\n";
      totalAttempts = attempts;
      fini();
      return 2;
    }
  }
  totalAttempts = attempts;
  if (verbose) {
    cout << "Done: " << filled << " filled out of " << attempts
         << " attempts (efficiency " << 100.0 * filled / attempts
         << "%)\n";
  }
  fini();
  return 0;
}
