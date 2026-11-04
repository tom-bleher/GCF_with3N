#include <iostream>
#include <unistd.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
long long totalAttempts = 0;
bool verbose = false;
bool useCM   = true;   // -C 0 to disable CM smearing
gcfNucleus * myInfo;
TRandom3 * myRand;
eNCrossSection * myCS;
QEGenerator_3N * myGen;
TTree * outtree;

// Tree variables
Double_t pe[3], pLead[3], p2[3], p3[3], pAm3[3];
Double_t weight;
Int_t N1_type, N2_type, N3_type;

void Usage()
{
  cerr << "Usage: ./code <Beam energy (GeV)> [Interaction #] <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-c: Specify eN cross section model (default cc1)\n"
       << "-G: Nucleon form factors: kelly (default), ye or dipole\n"
       << "-A: Target mass number (default 12)\n"
       << "-Z: Target charge number (default 6)\n"
       << "-s: CM Gaussian width sigCM in GeV/c per component (default 0.15)\n"
       << "-C: Toggle CM smearing: 1=on (default), 0=off\n"
       << "-S: Random seed (positive integer; default 0 = random, printed at start)\n"
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
  double Ebeam = atof(argv[1]);
  double u = atoi(argv[2]);
  nEvents = atoi(argv[4]);

  ffModel ffMod=kelly;
  csMethod csMeth=cc1;
  int targetA = 12;
  int targetZ = 6;
  double sigCM_GeV = 0.15;
  bool   sigCM_set = false;
  UInt_t seed = 0;

  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vc:A:Z:s:C:S:G:h")) != -1)
    switch(c)
      {

      case 'v':
	verbose = true;
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
      case 'A':
	targetA = atoi(optarg); break;
      case 'Z':
	targetZ = atoi(optarg); break;
      case 's':
	sigCM_GeV = atof(optarg); sigCM_set = true; break;
      case 'C':
	useCM = (atoi(optarg) != 0); break;
      case 'S':
	seed = strtoul(optarg, nullptr, 10); break;
      case 'h':
	Usage();
	return false;
      default:
	// getopt already printed "invalid option"
	Usage();
	return false;
      }

  // Initialize objects
  // Seed 0 (the default) draws a fresh seed; print it either way so any run
  // (e.g. a crashed parallel worker) can be replayed exactly with -S.
  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  cout << "Random seed: " << seed << "\n";
  myRand = new TRandom3(seed);
  myCS = new eNCrossSection(csMeth,ffMod);
  
  // Initialize generator
  myGen = new QEGenerator_3N(Ebeam, myCS, u, myRand);
  myGen->SetTargetNucleus(targetA, targetZ);
  if (sigCM_set) myGen->SetSigCM(sigCM_GeV);
  myGen->set_theta_k_maxmin(5,40);
  // Set up the tree
  // Open the output only after the generator is set up: exit() from a
  // setup error with an open TFile crashes in ROOT's teardown and buries
  // the real error message.
  outfile = new TFile(argv[3],"RECREATE");
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("N1_type",&N1_type,"N1_type/I");
  outtree->Branch("N2_type",&N2_type,"N2_type/I");
  outtree->Branch("N3_type",&N3_type,"N3_type/I");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pLead",pLead,"pLead[3]/D");
  outtree->Branch("p2",p2,"p2[3]/D");
  outtree->Branch("p3",p3,"p3[3]/D");
  outtree->Branch("pAm3",pAm3,"pAm3[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{

  TLorentzVector vk;
  TLorentzVector vLead;
  TLorentzVector v2;
  TLorentzVector v3;
  TLorentzVector vAm3;

  myGen->generate_event(weight, N1_type, N2_type, N3_type, vk, vLead, v2, v3, vAm3, useCM);

  pe[0] = vk.X();
  pe[1] = vk.Y();
  pe[2] = vk.Z();
  pLead[0] = vLead.X();
  pLead[1] = vLead.Y();
  pLead[2] = vLead.Z();
  p2[0] = v2.X();
  p2[1] = v2.Y();
  p2[2] = v2.Z();
  p3[0] = v3.X();
  p3[1] = v3.Y();
  p3[2] = v3.Z();
  pAm3[0] = vAm3.X();
  pAm3[1] = vAm3.Y();
  pAm3[2] = vAm3.Z();

  if(weight>0.){
    outtree->Fill();
  }
}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
  // Total generation attempts (including weight = 0 rejections, which are
  // not stored in the tree), as in genQE_3N_FSI: a 1-bin histogram so hadd
  // sums it across parallel chunks, for per-attempt weight normalization.
  TH1D hAttempts("hAttempts", "generation attempts", 1, 0., 1.);
  hAttempts.SetBinContent(1, (double)totalAttempts);
  hAttempts.Write();
  outfile->Delete("genTbuffer;*");
  outfile->Close();
}

int main(int argc, char ** argv)
{

  if (not init(argc,argv))
    return -1;

  // Count successful (weight > 0) events. Keep generating until the target
  // number of successful events has been written to the tree.
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
