#include <iostream>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TParameter.h"
#include "QEGeneratorFSI.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

string nucleusName(int Z, int N)
{
  int A = Z + N;
  if (Z==1  && N==1)   return "Deuterium (2H)";
  if (Z==2  && N==1)   return "Helium-3 (3He)";
  if (Z==1  && N==2)   return "Tritium (3H)";
  if (Z==2  && N==2)   return "Helium-4 (4He)";
  if (Z==6  && N==6)   return "Carbon-12 (12C)";
  if (Z==13 && N==14)  return "Aluminum-27 (27Al)";
  if (Z==18 && N==22)  return "Argon-40 (40Ar)";
  if (Z==20 && N==20)  return "Calcium-40 (40Ca)";
  if (Z==26 && N==30)  return "Iron-56 (56Fe)";
  if (Z==82 && N==126) return "Lead-208 (208Pb)";
  return "Z=" + to_string(Z) + " N=" + to_string(N) + " A=" + to_string(A);
}

string nnModelName(NNModel m)
{
  switch(m) {
    case AV18:      return "AV18";
    case AV4Pc:     return "AV4'";
    case N2LO_10:   return "N2LO (1.0 fm)";
    case N2LO_12:   return "N2LO (1.2 fm)";
    case N3LO_600:  return "N3LO (600 MeV)";
    case NV2_1a:    return "NV2-1a";
    case AV18_deut: return "AV18 (deuteron)";
    default:        return "unknown";
  }
}

string csMethodName(csMethod m)
{
  switch(m) {
    case onshell: return "on-shell";
    case cc1:     return "CC1";
    case cc2:     return "CC2";
    default:      return "unknown";
  }
}

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
QEGeneratorFSI * myGen;
TTree * outtree;
bool doLC;
bool doCoul = false;
double deltaECoul = 0;

// Tree variables. pLead/pRec are the GENERATED (pre-FSI) nucleons; the
// post-FSI final state is purely the flat particle list below (GENIE's
// cascade output verbatim -- no lead/recoil designation, no absorption
// special-casing).
Double_t pe[3], pLead[3], pRec[3], pAm2[3], pRel[3], q[3];
Double_t weight;
Int_t lead_type, rec_type;
Int_t doFSI_flag;
UInt_t seed = 0;
bool doMF = false;      // -F: single-nucleon mean-field knockout instead of an SRC pair
double beamE = 0.;      // for q in MF mode
Int_t nPart;
std::vector<int>    part_pdg, part_parentRole, part_rescatterCode;
std::vector<double> part_px, part_py, part_pz, part_E;

void Usage()
{
  cerr << "Usage: ./genQE_FSI <Z> <N> <Beam energy (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-u: Specify NN interaction (default AV18)\n"
       << "-c: Specify eN cross section model (default cc1)\n"
       << "-s: Specify sigma_CM [GeV/c]\n"
       << "-E: Specify E* [GeV]\n"
       << "-M: Use randomized E* according to Barack's values\n"
       << "-k: Specify pRel hard cutoff [GeV/c]\n"
       << "-O: Turn on peaking radiation\n"
       << "-C: Turn on coulomb correction\n"
       << "-r: Randomize nuclear properties\n"
       << "-l: Use Lightcone cross section\n"
       << "-n: Disable FSI (PWIA, same generator)\n"
       << "-F: Mean-field single-nucleon knockout from a global Fermi gas, kF in GeV/c (e.g. 0.25)\n"
       << "-S: Random seed (positive integer; default 0 = random, printed at start)\n"
       << "-h: Print this message and exit\n\n\n";
}

bool init(int argc, char ** argv)
{
  int numargs = 6;
 
  if (argc < numargs)
    {
      Usage();
      return false;
    }

  // Read in the arguments
  int Z = atoi(argv[1]);
  int N = atoi(argv[2]);
  double Ebeam = atof(argv[3]);
  nEvents = atoi(argv[5]);

  ffModel ffMod=kelly;
  
  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  char * uType = (char *)"AV18";
  csMethod csMeth=cc1;
  double sigCM = 0.;
  bool do_sigCM = false;
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  double kCut = 0.25;
  bool do_kCut = false;
  bool doRad = false;
  bool rand_flag = false;
  bool noFSI = false;
  double mfKF = 0.;
  doLC = false;
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:c:s:E:Mk:OCrlnF:S:h")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
	break;
      case 'u':
	uType = optarg;
	break;
      case 'c':
	if (strcmp(optarg,"onshell")==0)
	  csMeth=onshell;
	else if ((strcmp(optarg,"cc1")==0) or (atoi(optarg)==1))
	  csMeth=cc1;
	else if ((strcmp(optarg,"cc2")==0) or (atoi(optarg)==2))
	  csMeth=cc2;
	else
	  {
	    cerr << "Invalid cross section designation. Allowed values are onshell, cc1 and cc2. Aborting...\n";
	    return false;
	  }
	break;
      case 's':
	do_sigCM = true;
	sigCM = atof(optarg);
	break;
      case 'E':
	do_Estar = true;
	Estar = atof(optarg);
	break;
      case 'M':
        do_Estar = true;
	Estar = 0.01732;
        do_sigmaE = true;
	sigmaE = 0.009571;
	break;
      case 'k':
	do_kCut = true;
	kCut = atof(optarg);
	break;
      case 'O':
	doRad = true;
	break;
      case 'C':
	doCoul = true;
	break;
      case 'r':
        rand_flag = true;
	break;
      case 'l':
	doLC = true;
	break;
      case 'n':
	noFSI = true;
	break;
      case 'F':
	doMF = true;
	mfKF = atof(optarg);
	if (mfKF <= 0. || mfKF > 0.5) { cerr << "Invalid -F kF (GeV/c, e.g. 0.25).\n"; return false; }
	break;
      case 'S':
	seed = strtoul(optarg, nullptr, 10);
	break;
      case 'h':
	Usage();
	return false;
      default:
	// getopt already printed "invalid option"
	Usage();
	return false;
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,uType);
  // Seed 0 (the default) draws a fresh seed; print and store it either way. GENIE's
  // generator is seeded from myRand (GenieFSIHelpers), so FSI replays too.
  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  cout << "Random seed: " << seed << "\n";
  myRand = new TRandom3(seed);
  myCS = new eNCrossSection(csMeth,ffMod);

  if (rand_flag)
    myInfo->randomize(myRand);
  if (do_sigCM)
    myInfo->set_sigmaCM(sigCM);
  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);
  
  // Initialize FSI-enabled generator
  myGen = new QEGeneratorFSI(Ebeam + deltaECoul, myInfo, myCS, myRand);
  myGen->EnableFSI(!noFSI);
  doFSI_flag = noFSI ? 0 : 1;
  beamE = Ebeam + deltaECoul;
  if (doMF) {
    myGen->SetFermiMomentumFG(mfKF);
    myGen->SetFermiGasMode(QEGeneratorFSI::kGlobalFG);
  }
  
  if ((Z == 1) and (N == 1))
    myGen->set_deuteron();
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  if (do_kCut)
    myGen->set_pRel_cut(kCut);
  if (doRad)
    myGen->set_doRad(true);
  if (doCoul)
    myGen->set_doCoul(true);
  
  // Set up the tree
  // Open the output only after the generator is set up: exit() from a
  // setup error with an open TFile crashes in ROOT's teardown and buries
  // the real error message.
  outfile = new TFile(argv[4],"RECREATE");
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("pRel",pRel,"pRel[3]/D");
  outtree->Branch("q",q,"q[3]/D");
  outtree->Branch("doFSI",&doFSI_flag,"doFSI/I");
  outtree->Branch("weight",&weight,"weight/D");
  outtree->Branch("nPart",    &nPart,       "nPart/I");
  outtree->Branch("part_pdg",           &part_pdg);
  outtree->Branch("part_parentRole",    &part_parentRole);
  outtree->Branch("part_rescatterCode", &part_rescatterCode);
  outtree->Branch("part_px",            &part_px);
  outtree->Branch("part_py",            &part_py);
  outtree->Branch("part_pz",            &part_pz);
  outtree->Branch("part_E",             &part_E);

  // Print run configuration
  cout << "\n"
       << "========================================\n"
       << "  GCF Event Generator Configuration\n"
       << "========================================\n"
       << "  Nucleus:        " << nucleusName(Z,N) << "\n"
       << "  Beam energy:    " << Ebeam << " GeV\n"
       << "  NN interaction: " << nnModelName(myInfo->get_SRC()->get_InteractionType()) << "\n"
       << "  eN cross sec:   " << csMethodName(csMeth) << "\n"
       << "  Events:         " << nEvents << "\n"
       << "  Output:         " << argv[4] << "\n"
       << "----------------------------------------\n"
       << "  Q2 range:       [" << myGen->get_QSqmin() << ", " << myGen->get_QSqmax() << "] (GeV/c)^2\n";
  if (myGen->get_xBmin() != 0. || myGen->get_xBmax() != 0.)
    cout << "  xB range:       [" << myGen->get_xBmin() << ", " << myGen->get_xBmax() << "]\n";
  else
    cout << "  xB range:       (no cut)\n";
  cout << "  sigma_CM:       " << myInfo->get_sigmaCM() << " GeV/c\n"
       << "  pRel cutoff:    " << kCut << " GeV/c\n"
       << "  FSI:            " << (noFSI ? "disabled (PWIA)" : "enabled (hN2018)") << "\n";
  if (doMF)    cout << "  Mean field:     single nucleon, global Fermi gas, kF=" << mfKF << " GeV/c\n";
  if (doLC)    cout << "  Lightcone:      ON\n";
  if (doCoul)  cout << "  Coulomb corr:   ON\n";
  if (doRad)   cout << "  Radiation:      ON\n";
  if (rand_flag) cout << "  Randomization:  ON\n";
  if (!noFSI)  fsi::PrintGenieXMLConfig(cout);
  cout << "========================================\n\n";

  return true;
  
}

void evnt(int event)
{

  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector vRec;
  TLorentzVector vAm2;
  TVector3 vRel;
  TLorentzVector vq;
  double Estar;

  if (doMF) {
    myGen->generate_event_MF(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);
    vq = TLorentzVector(0., 0., beamE, beamE) - vk;
    vRel.SetXYZ(0., 0., 0.);
  }
  else if (doLC)
    myGen->generate_event_lightcone(weight, lead_type, rec_type, vk, vLead, vRec, vAm2);
  else
    myGen->generate_event(weight, lead_type, rec_type, vk, vLead, vRec, vAm2, vRel, vq, Estar);
  
  pe[0] = vk.X();
  pe[1] = vk.Y();
  pe[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();
  pRel[0] = vRel.X();
  pRel[1] = vRel.Y();
  pRel[2] = vRel.Z();
  q[0] = vq.X();
  q[1] = vq.Y();
  q[2] = vq.Z();

  // Post-FSI final state: flat particle list.
  part_pdg.clear(); part_parentRole.clear(); part_rescatterCode.clear();
  part_px.clear();  part_py.clear();  part_pz.clear();  part_E.clear();
  for (const auto &p : myGen->GetFinalStateParticles()) {
    part_pdg.push_back(p.pdg);
    part_parentRole.push_back(p.parentRole);
    part_rescatterCode.push_back(p.rescatterCode);
    part_px.push_back(p.p4.X()); part_py.push_back(p.p4.Y());
    part_pz.push_back(p.p4.Z()); part_E.push_back(p.p4.T());
  }
  nPart = static_cast<int>(part_pdg.size());

  if (weight > 0.)
    outtree->Fill();
  
}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
  // Attempts (exactly nEvents; only weight > 0 is stored) for per-attempt normalization, and the seed.
  TH1D hAttempts("hAttempts", "generation attempts", 1, 0., 1.);
  hAttempts.SetBinContent(1, (double)nEvents);
  hAttempts.Write();
  TParameter<Long64_t>("seed", seed).Write();
  outfile->Delete("genTbuffer;*");
  outfile->Close();
}

int main(int argc, char ** argv)
{

  if (not init(argc,argv))
    return -1;

  for (int event=0; event < nEvents; event++)
    {
      if ((event %100000==0) && (verbose))
	cout << "Working on event " << event << "\n";

      evnt(event);
    }

  fini();
  
  return 0;
  
}
