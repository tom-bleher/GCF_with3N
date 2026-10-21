#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "constants.hh"
#include "helpers.hh"
#include "QEGenerator_3N.hh"

using namespace std;
eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

void Usage()
{
  cerr << "Usage: ./code [potential] <output file>\n\n";
}

int main(int argc, char ** argv)
{

  if(argc != 3)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  // Initialize objects
  myRand = new TRandom3(0);
  ffModel ffMod=kelly;
  csMethod csMeth=cc1;
  myCS = new eNCrossSection(csMeth,ffMod);
  // Initialize generator
  myGen = new QEGenerator_3N(6.0, myCS, atoi(argv[1]), myRand);
  char directory[100];
  strcpy(directory,argv[2]);
  char temp_name[100];

  const int array_bins = 500;
  const double X_min = -1/sqrt(3);
  const double X_max =  1/sqrt(3);
  const double Y_min =  0;
  const double Y_max =  1;       
  const double Z_tot_min = 0;
  const double Z_tot_max = 15;
  const double Z_miss_min = 0;
  const double Z_miss_max = 5;
  const double f1s[7] = {0.490,0.490,0.490,0.255,0.255,0.020,0.33333};
  const double f2s[7] = {0.490,0.255,0.020,0.255,0.490,0.490,0.33333};

  /////////////////////////////////////////////////////////////
  //2D arrays as a function of p Total
  /////////////////////////////////////////////////////////////
  for(int S = 0; S <= 50; S+=1){    
    double p_tot = (double)S * 3.0/50.0;
    sprintf(temp_name,"%s/Shape_3N_pp_%i_ptot.txt",directory,S);
    ofstream ppFile(temp_name);
    for(int i = 0; i < array_bins+1; i++){
      for(int j = 0; j < array_bins+1; j++){	
	double X = (((double)i/array_bins)*(X_max-X_min)) + X_min;
	double Y = (((double)j/array_bins)*(Y_max-Y_min)) + Y_min;
	double f1 = (1 - Y) / 2.0;
	double f2 = (sqrt(3)*X + Y + 1.0) / 4.0;
	double f3 = 1 - f1 - f2;	
	if((f1>0.5) || (f2>0.5) || (f3>0.5) || (f1<0.0) || (f2<0.0) || (f3<0.0)){
	  ppFile << 0 << " ";
	  continue;
	}
	double k_tot = p_tot / GeVfm;
	double k_1 = f1 * k_tot;
	double k_2 = f2 * k_tot;
	double k_3 = f3 * k_tot;
	double p_1 = k_1 * GeVfm;
	double p_2 = k_2 * GeVfm;
	double p_3 = k_3 * GeVfm;
	ppFile << myGen->get_rho_ptot_f1f2f3(k_1,k_2,k_3) << " ";		
      }
      ppFile << "\n";
    }
    ppFile.close();
  }

  /////////////////////////////////////////////////////////////
  //2D arrays as a function of p Miss
  /////////////////////////////////////////////////////////////
  for(int S = 0; S <= 50; S+=1){    
    double p_miss = (double)S / 50.0;
    sprintf(temp_name,"%s/Shape_3N_pp_%i_pmiss.txt",directory,S);
    ofstream ppFile(temp_name);
    for(int i = 0; i < array_bins+1; i++){
      for(int j = 0; j < array_bins+1; j++){	
	double X = (((double)i/array_bins)*(X_max-X_min)) + X_min;
	double Y = (((double)j/array_bins)*(Y_max-Y_min)) + Y_min;
	double f1 = (1 - Y) / 2.0;
	double f2 = (sqrt(3)*X + Y + 1.0) / 4.0;
	double f3 = 1 - f1 - f2;	
	if((f1>0.5) || (f2>0.5) || (f3>0.5) || (f1<=0.0) || (f2<0.0) || (f3<0.0)){
	  ppFile << 0 << " ";
	  continue;
	}
	double k_miss = p_miss / GeVfm;
	double k_1 = k_miss;
	double k_tot = k_1 / f1;
	double k_2 = f2 * k_tot;
	double k_3 = f3 * k_tot;
	double p_1 = k_1 * GeVfm;
	double p_2 = k_2 * GeVfm;
	double p_3 = k_3 * GeVfm;
	// density per (p_miss, f1, f2): p_tot = p_miss/f1, so dp_tot = dp_miss/f1
	ppFile << myGen->get_rho_ptot_f1f2f3(k_1,k_2,k_3) / f1 << " ";		
      }
      ppFile << "\n";
    }
    ppFile.close();
  }


  /////////////////////////////////////////////////////////////
  //1D arrays as a function of p Total
  /////////////////////////////////////////////////////////////
  for(int S = 0; S < 7; S++){    
    double f1 = f1s[S];
    double f2 = f2s[S];
    double f3 = 1 - f1 - f2;
    sprintf(temp_name,"%s/OneDim_3N_pp_%i_ptot.txt",directory,S+1);
    ofstream ppFile(temp_name);
    for(int i = 0; i < array_bins+1; i++){
      double Z = (((double)i/array_bins)*(Z_tot_max-Z_tot_min)) + Z_tot_min;
      double k_tot = Z;
      double k_1 = f1 * k_tot;
      double k_2 = f2 * k_tot;
      double k_3 = f3 * k_tot;
      double p_1 = k_1 * GeVfm;
      double p_2 = k_2 * GeVfm;
      double p_3 = k_3 * GeVfm;
      double p_tot = k_tot * GeVfm;
      ppFile << p_tot << " " << myGen->get_rho_ptot_f1f2f3(k_1,k_2,k_3) << "\n";			
    }
    ppFile.close();
  }  

  /////////////////////////////////////////////////////////////
  //1D arrays as a function of p Miss
  /////////////////////////////////////////////////////////////
  for(int S = 0; S < 7; S++){    
    double f1 = f1s[S];
    double f2 = f2s[S];
    double f3 = 1 - f1 - f2;
    sprintf(temp_name,"%s/OneDim_3N_pp_%i_pmiss.txt",directory,S+1);
    ofstream ppFile(temp_name);
    for(int i = 0; i < array_bins+1; i++){
      double Z = (((double)i/array_bins)*(Z_miss_max-Z_miss_min)) + Z_miss_min;
      double k_miss = Z;
      double k_1 = k_miss;
      double k_tot = k_1/f1;
      double k_2 = f2 * k_tot;
      double k_3 = f3 * k_tot;
      double p_1 = k_1 * GeVfm;
      double p_2 = k_2 * GeVfm;
      double p_3 = k_3 * GeVfm;
      double p_tot = k_tot * GeVfm;
      double p_miss = k_miss * GeVfm;
      // density per (p_miss, f1, f2): p_tot = p_miss/f1, so dp_tot = dp_miss/f1
      ppFile << p_miss << " " << myGen->get_rho_ptot_f1f2f3(k_1,k_2,k_3) / f1 << "\n";			
    }
    ppFile.close();
  }  
  
  return 0;  
}
  
