#include <iostream>
#include <unistd.h>
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "ionGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
UInt_t seed = 0;
TFile * outfile;
bool verbose = false;
bool inv_kin = false;
TVector3 boost_vector;
gcfNucleus * myInfo;
TRandom3 * myRand;
NNCrossSection * myCS;
ionGenerator * myGen;
TTree * outtree;

// Tree variables
Double_t p3[3], p4[3], pRec[3], pAm2[3];
Double_t weight;
Int_t lead_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genIon <Z> <N> <Beam energy per nucleon (GeV)> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-I: Boost to inverse kinematics. In this case, the given beam energy is the kinetic energy per nucleon"
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
  outfile = new TFile(argv[4],"RECREATE");
  nEvents = atoi(argv[5]);

  // Optional flags
  csParam myParam = Panin;
  bool custom_ps = false;
  char * phase_space;
  char * uType = (char *)"AV18";
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";
  
  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:IS:h")) != -1)
    switch(c)
      {
	
      case 'v':
	verbose = true;
	break;
      case 'P':
	custom_ps = true;
	phase_space = optarg;
	break;
      case 'S':
	seed = strtoul(optarg, nullptr, 10);
	break;
      case 'h':
	Usage();
	return false;
      case 'I':
	inv_kin = true;
	break;
      default:
	abort();
	
      }

  // Initialize objects
  myInfo = new gcfNucleus(Z,N,uType);
  // Seed 0 (the default) draws a fresh seed; print and store it either way so
  // any run can be reproduced.
  if (seed == 0) { TRandom3 seeder(0); seed = 1 + seeder.Integer(0xfffffffe); }
  cout << "Random seed: " << seed << "\n";
  myRand = new TRandom3(seed);
  myCS = new NNCrossSection(myParam);
  
  if (inv_kin)
    {
      double mA = myInfo->get_mA();
      int Anum = Z + N;
      double EA = mA + Ebeam*Anum;
      double pA = sqrt(sq(EA) - sq(mA));
      double p1 = mN*pA/mA;
      double E1 = sqrt(sq(mN) + sq(p1));
      TLorentzVector v1_target(0.,0.,p1,E1);
      boost_vector = -v1_target.BoostVector();
      Ebeam = E1 - mN;
    }
  
  // Initialize generator
  myGen = new ionGenerator(Ebeam, myInfo, myCS, myRand);
  if ((Z == 1) and (N == 1))
    myGen->set_deuteron();
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);

  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("lead_type",&lead_type,"lead_type/I");
  outtree->Branch("rec_type",&rec_type,"rec_type/I");
  outtree->Branch("p3",p3,"p3[3]/D");
  outtree->Branch("p4",p4,"p4[3]/D");
  outtree->Branch("pRec",pRec,"pRec[3]/D");
  outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
  outtree->Branch("weight",&weight,"weight/D");
  
  return true;
  
}

void evnt(int event)
{
  
  TLorentzVector v3;
  TLorentzVector v4;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, lead_type, rec_type, v3, v4, vRec, vAm2);

  if (inv_kin)
    {
      v3.Boost(boost_vector);
      v4.Boost(boost_vector);
      vRec.Boost(boost_vector);
      vAm2.Boost(boost_vector);
      v3.RotateX(M_PI);
      v4.RotateX(M_PI);
      vRec.RotateX(M_PI);
      vAm2.RotateX(M_PI);
    }
  
  p3[0] = v3.X();
  p3[1] = v3.Y();
  p3[2] = v3.Z();
  p4[0] = v4.X();
  p4[1] = v4.Y();
  p4[2] = v4.Z();
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
