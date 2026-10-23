#include <iostream>
#include <unistd.h>
#include "TFile.h"
#include "TTree.h"
#include "CCQEGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
CCQEGenerator * myGen;
TTree * outtree;

// Tree variables
Double_t pk[3], pLead[3], pRec[3], pAm2[3];
Double_t weight, Eneutrino;
Int_t rec_type;

void Usage()
{
  cerr << "Usage: ./genCCQE <Z> <N> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-G: Nucleon form factors: kelly (default), ye or dipole\n"
       << "-M: Use randomized E* according to Barack's values\n"
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

  // Read in the arguments
  int Z = atoi(argv[1]);
  int N = atoi(argv[2]);
  outfile = new TFile(argv[3],"RECREATE");
  nEvents = atoi(argv[4]);

  csMethod csMeth=cc1;
  ffModel ffMod=kelly;
  
  // Optional flags
  bool custom_ps = false;
  char * phase_space;
  char * uType = (char *)"AV18";
  double Estar = 0.;
  bool do_Estar = false;
  double sigmaE = 0.;
  bool do_sigmaE = false;
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";

  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:MG:h")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
	break;
      case 'G':
	if (!parse_ffModel(optarg,ffMod))
	  {
	    cerr << "Invalid form-factor model. Allowed values are kelly, ye and dipole. Aborting...\n";
	    return false;
	  }
	break;
      case 'M':
        do_Estar = true;
	Estar = 0.01732;
        do_sigmaE = true;
	sigmaE = 0.009571;
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
  
  if (do_Estar)
    myInfo->set_Estar(Estar);
  if (do_sigmaE)
    myInfo->set_sigmaE(sigmaE);
  
  // Initialize generator
  myGen = new CCQEGenerator(myInfo, myCS, myRand);
  if ((Z == 1) and (N == 1))
    myGen->set_deuteron();
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  
  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("Eneutrino",&Eneutrino,"Eneutrino/D");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("pk",pk,"pk[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{

  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, Eneutrino, rec_type, vk, vLead, vRec, vAm2);
  
  pk[0] = vk.X();
  pk[1] = vk.Y();
  pk[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();

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
