#include <iostream>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "QEGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
QEGenerator * myGen;
TTree * outtree;
bool doLC;
bool doCoul = false;
double deltaECoul = 0;

// Tree variables
Double_t pe[3], pLead[3], pRec[3], pAm2[3], pRel[3], q[3];
Double_t weight;
Int_t lead_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genQE <Z> <N> <Beam energy (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-u: Specify NN interaction (default AV18)\n"
       << "-c: Specify eN cross section model (default cc1)\n"
       << "-G: Nucleon form factors: kelly (default), ye or dipole\n"
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
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:c:s:E:Mk:OCrlG:h")) != -1)
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
      case 'G':
	if (!parse_ffModel(optarg,ffMod))
	  {
	    cerr << "Invalid form-factor model. Allowed values are kelly, ye and dipole. Aborting...\n";
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
  
  // Initialize generator
  myGen = new QEGenerator(Ebeam + deltaECoul, myInfo, myCS, myRand);
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
  outtree->Branch("weight",&weight,"weight/D");
  
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
