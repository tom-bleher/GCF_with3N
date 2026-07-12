#include <iostream>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
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

// Tree variables
Double_t pe[3], pLead[3], pRec[3], pAm2[3], pRel[3], q[3];
Double_t pLead_pre[3], pRec_pre[3];   // pre-FSI lead/recoil 3-momenta
Double_t weight;
Int_t lead_type, rec_type;
Int_t doFSI_flag;
Int_t lead_absorbed, rec_absorbed;

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
  outfile = new TFile(argv[4],"RECREATE");
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
  doLC = false;
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:c:s:E:Mk:OCrlh")) != -1)
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
	    return -1;
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
      case 'h':
	Usage();
	return false;
      default:
	abort();
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,uType);
  myRand = new TRandom3(0);
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
  myGen->EnableFSI(true);  // FSI is enabled by default, but explicit for clarity
  doFSI_flag = 1;
  
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
  outtree->Branch("pLead_pre",pLead_pre,"pLead_pre[3]/D");
  outtree->Branch("pRec_pre",pRec_pre,"pRec_pre[3]/D");
  outtree->Branch("doFSI",&doFSI_flag,"doFSI/I");
  // Per-leg absorption flags (diagnostics only). True (1) when the ORIGINAL
  // nucleon ceased to exist in the cascade (typically NN absorption). The
  // event is NOT killed; the GCF primary weight is kept. The leg's saved
  // 4-momentum is then its leading nucleon ejectum, promoted from the
  // secondaries (what a detector would reconstruct as the leg); it is zero
  // only if the absorption produced no nucleon ejecta at all. Topology
  // selection must use the FINAL state (particle types and momenta), never
  // these flags -- an experiment cannot veto on absorption.
  outtree->Branch("lead_absorbed",&lead_absorbed,"lead_absorbed/I");
  outtree->Branch("rec_absorbed", &rec_absorbed, "rec_absorbed/I");
  outtree->Branch("weight",&weight,"weight/D");

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
       << "  FSI:            enabled (hN2018)\n";
  if (doLC)    cout << "  Lightcone:      ON\n";
  if (doCoul)  cout << "  Coulomb corr:   ON\n";
  if (doRad)   cout << "  Radiation:      ON\n";
  if (rand_flag) cout << "  Randomization:  ON\n";
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

  if (doLC)
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

  // Pre-FSI lead/recoil (saved before the GENIE intranuclear cascade
  // modifies them in ApplyFSI). When FSI is disabled, fall back to post.
  if (doFSI_flag) {
    const TLorentzVector & preL = myGen->GetPreFSILead();
    const TLorentzVector & preR = myGen->GetPreFSIRecoil();
    pLead_pre[0] = preL.X(); pLead_pre[1] = preL.Y(); pLead_pre[2] = preL.Z();
    pRec_pre[0]  = preR.X(); pRec_pre[1]  = preR.Y(); pRec_pre[2]  = preR.Z();
    lead_absorbed = myGen->IsLeadAbsorbed() ? 1 : 0;
    rec_absorbed  = myGen->IsRecAbsorbed()  ? 1 : 0;
  } else {
    pLead_pre[0] = pLead[0]; pLead_pre[1] = pLead[1]; pLead_pre[2] = pLead[2];
    pRec_pre[0]  = pRec[0];  pRec_pre[1]  = pRec[1];  pRec_pre[2]  = pRec[2];
    lead_absorbed = 0;
    rec_absorbed  = 0;
  }

  if (weight > 0.)
    outtree->Fill();
  
}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
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
