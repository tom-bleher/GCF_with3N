#include <iostream>
#include <fstream>
#include <unistd.h>
#include "TFile.h"
#include "TTree.h"
#include "photoGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

int nEvents;
TFile * outfile;
bool verbose = false;
bool do_ascii = false;
gcfNucleus * myInfo;
TRandom3 * myRand;
photoCrossSection * myCS;
photoGenerator * myGen;
TTree * outtree;
ofstream * asciiWriter = nullptr;

int numOut = 0;

// Tree variables
Double_t pMeson[3], pBaryon[3], pRec[3], pAm2[3];
Double_t weight, Ephoton, mMeson, mBaryon;
Int_t meson_type, baryon_type, rec_type;

void Usage()
{
  cerr << "Usage: ./genPhoto <Z> <N> <path/to/output.root> <# of events>\n\n"
       << "Optional flags:\n"
       << "-v: Verbose\n"
       << "-P: Use text file to specify phase space\n"
       << "-u: Specify NN interaction (default AV18)\n"
       << "-k: Specify pRel hard cutoff [GeV/c]\n"
       << "-R: Specify the reaction channel, default pim+proton. (pim, rho0, rhom, omega, phi, phin, delta++ [Delta++ pi-], delta+pi- [Delta+ pi-])\n"
       << "-A: Use the energy spectrum of the amorphous GlueX radiator rather than the diamond radiator\n"
       << "-B: Specify a fixed beam energy (default is HallD spectrum)\n"
       << "-l: Use the lightcone cross section to generate events\n"
       << "-h: Print this message and exit\n\n\n";
}

void asciiWriteProton(int num, TLorentzVector v)
{
  (*asciiWriter) << "   " << num << " " << pCode_geant << " " << mN << endl;
  (*asciiWriter) << "   " << 1 << " " << v.X() << " " << v.Y() << " " << v.Z() << " " << v.T() << endl;
}

void asciiWriteNeutron(int num, TLorentzVector v)
{
  (*asciiWriter) << "   " << num << " " << nCode_geant << " " << mN << endl;
  (*asciiWriter) << "   " << 0 << " " << v.X() << " " << v.Y() << " " << v.Z() << " " << v.T() << endl;
}

void asciiWritePiPlus(int num, TLorentzVector v)
{
  (*asciiWriter) << "   " << num << " " << pipCode_geant << " " << mpip << endl;
  (*asciiWriter) << "   " << 1 << " " << v.X() << " " << v.Y() << " " << v.Z() << " " << v.T() << endl;
}

void asciiWritePiMinus(int num, TLorentzVector v)
{
  (*asciiWriter) << "   " << num << " " << pimCode_geant << " " << mpip << endl;
  (*asciiWriter) << "   " << -1 << " " << v.X() << " " << v.Y() << " " << v.Z() << " " << v.T() << endl;
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

  // Optional flags
  char * asciiFile;
  bool custom_ps = false;
  char * phase_space;
  char * uType = (char *)"AV18";
  double kCut = 0.25;
  bool do_kCut = false;
  char * react;
  spectrum mySpectrum = diamond;
  reaction myReaction = proton_piMinus;
  bool usefixedE=false;
  double fixedE=0;
  bool doLC = false;
  if ((Z == 1) and (N == 1))
    uType = (char *)"AV18_deut";

  int c;
  while ((c = getopt (argc-numargs+1, &argv[numargs-1], "vP:u:k:R:AB:lh")) != -1)
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
      case 'k':
	do_kCut = true;
	kCut = atof(optarg);
	break;
      case 'R':
	react = optarg;
	if ((strcmp(react, "pim")==0) or (strcmp(react, "proton_piMinus")==0))
	  myReaction=proton_piMinus;
	else if ((strcmp(react, "rho0")==0) or (strcmp(react, "proton_rho0")==0))
	  myReaction=proton_rho0;
	else if ((strcmp(react, "rhom")==0) or (strcmp(react, "rho-")==0) or (strcmp(react, "proton_rhoMinus")==0))
	  myReaction=proton_rhoMinus;
	else if ((strcmp(react, "omega")==0) or (strcmp(react, "proton_omega")==0))
	  myReaction=proton_omega;
	else if ((strcmp(react, "phi")==0) or (strcmp(react, "proton_phi")==0))
	  myReaction=proton_phi;
	else if ((strcmp(react, "phin")==0) or (strcmp(react, "phiN")==0) or strcmp(react, "neutron_phi")==0)
	  myReaction=neutron_phi;
	else if ((strcmp(react, "deltapp")==0) or (strcmp(react, "delta++")==0) or strcmp(react, "DeltaPlusPlus_piMinus")==0)
	  myReaction=DeltaPlusPlus_piMinus;
	else if ((strcmp(react, "deltappim")==0) or (strcmp(react, "delta+pi-")==0) or strcmp(react, "DeltaPlus_piMinus")==0)
	  myReaction=DeltaPlus_piMinus;
	else
	  {
	    cerr << "This reaction is not yet implemented.\n";
	    exit(-1);
	  }
	break;
      case 'A':
	mySpectrum = amorphous;
        break;
      case 'B':
	usefixedE=true;
	fixedE=atof(optarg);
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
  myCS = new photoCrossSection();

  // Initialize generator
  myGen = new photoGenerator(myInfo, myCS, myRand, mySpectrum, myReaction);
  if ((Z == 1) and (N == 1))
    myGen->set_deuteron();
  if (custom_ps)
    myGen->parse_phase_space_file(phase_space);
  if (do_kCut)
    myGen->set_pRel_cut(kCut);
  if (usefixedE)
    myGen->setfixedE(fixedE);
  if (verbose)
    myGen->print_beam_info();
  if (doLC)
    myGen->setDoLC(true);

  // Set up the tree
  outfile->cd();
  outtree = new TTree("genTbuffer","Generator Tree");
  if (not do_ascii)
    {
      outtree->Branch("Ephoton",&Ephoton,"Ephoton/D");
      outtree->Branch("meson_type",&meson_type,"meson_type/I");
      outtree->Branch("mMeson",&mMeson,"mMeson/D");
      outtree->Branch("baryon_type",&baryon_type,"baryon_type/I");
      outtree->Branch("mBaryon",&mBaryon,"mBaryon/D");
      outtree->Branch("rec_type",&rec_type,"rec_type/I");
      outtree->Branch("pMeson",pMeson,"pMeson[3]/D");
      outtree->Branch("pBaryon",pBaryon,"pBaryon[3]/D");
      outtree->Branch("pRec",pRec,"pRec[3]/D");
      outtree->Branch("pAm2",pAm2,"pAm2[3]/D");
    }
  outtree->Branch("weight",&weight,"weight/D");

  if (do_ascii)
    asciiWriter = new ofstream(asciiFile);

  return true;

}

void evnt(int event)
{

  TLorentzVector vMeson;
  TLorentzVector vBaryon;
  TLorentzVector vRec;
  TLorentzVector vAm2;

  myGen->generate_event(weight, Ephoton, meson_type, mMeson, baryon_type, mBaryon, rec_type, vMeson, vBaryon, vRec, vAm2);

  pMeson[0] = vMeson.X();
  pMeson[1] = vMeson.Y();
  pMeson[2] = vMeson.Z();
  pBaryon[0] = vBaryon.X();
  pBaryon[1] = vBaryon.Y();
  pBaryon[2] = vBaryon.Z();
  pRec[0] = vRec.X();
  pRec[1] = vRec.Y();
  pRec[2] = vRec.Z();
  pAm2[0] = vAm2.X();
  pAm2[1] = vAm2.Y();
  pAm2[2] = vAm2.Z();

  if (weight > 0.)
    {
      outtree->Fill();
      if (do_ascii)
	{
	  (*asciiWriter) << "1 " << numOut << " 3 "<< Ephoton << endl;
	  switch(meson_type)
	    {
	    case pipCode:
	      asciiWritePiPlus(1,vMeson);
	      break;
	    case pimCode:
	      asciiWritePiMinus(1,vMeson);
	      break;
	    default:
	      cout << "Cannot write out this meson to ASCII format. Aborting...";
	      abort();
	    }
	  switch(baryon_type)
	    {
	    case pCode:
	      asciiWriteProton(2,vBaryon);
	      break;
	    case nCode:
	      asciiWriteNeutron(2,vBaryon);
	      break;
	    default:
	      cout << "Cannot write out this baryon to ASCII format. Aborting...";
	      abort();
	    }
	  switch(rec_type)
	    {
	    case pCode:
	      asciiWriteProton(3,vRec);
	      break;
	    case nCode:
	      asciiWriteNeutron(3,vRec);
	      break;
	    default:
	      cout << "Cannot write out this recoil to ASCII format. Aborting...";
	      abort();
	    }

	}
      numOut++;
    }

}

void fini()
{
  outtree->SetName("genT");
  outtree->Write();
  outfile->Delete("genTbuffer;*");
  outfile->Close();

  if ( asciiWriter ) delete asciiWriter;
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
