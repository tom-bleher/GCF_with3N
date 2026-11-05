#ifndef __GCF_NUCLEUS_H__
#define __GCF_NUCLEUS_H__

#include "TRandom3.h"
#include "gcfSRC.hh"

class TRandom3;

class gcfNucleus
{
 public:
  gcfNucleus(int thisZ, int thisN, char* uType);
  gcfNucleus(int thisZ, int thisN, NNModel uType);
  gcfNucleus(int thisZ, int thisN, gcfSRC * thisSRC);
  ~gcfNucleus();
  double get_S(double krel, int l_type, int r_type);
  int get_Z();
  int get_N();
  int get_A();
  double get_mA();
  double get_mbar();
  double get_mAmpp();
  double get_mAmpn();
  double get_mAmnn();
  // Mass of the A-1 residual after removing a SINGLE nucleon (mean-field
  // knockout). _p = proton removed, _n = neutron removed. Returns 0 for
  // nuclei whose A-1 masses are not tabulated (only carbon is, for now).
  double get_mAm1p();
  double get_mAm1n();
  double get_mAmpp_random(double &Estar, TRandom3* myRand);
  double get_mAmpn_random(double &Estar, TRandom3* myRand);
  double get_mAmnn_random(double &Estar, TRandom3* myRand);
  double get_sigmaCM();
  double get_d_sigmaCM() const { return d_sigmaCM; }
  double get_Estar_max() const { return Estar_max; }
  double get_Estar();
  double get_sigmaE();
  double get_Estar_random(TRandom3* myRand);
  gcfSRC * get_SRC();
  bool get_Estar_randomization();
  void randomize(TRandom3* myRand);
  void randomize_sigmaCM(TRandom3* myRand);
  void randomize_Estar(TRandom3* myRand);
  
  void set_Nucleus(int thisZ, int thisN);
  void setCustomValues(double newSigma, double newEstar, double newCpp0, double Cnn0, double newCpn0, double newCpn1);
  void set_Interaction(NNModel thisNNType);
  void set_Interaction(char* thisNNType);
  void set_sigmaCM(double newSigma);
  void set_Estar(double newEstar);
  void set_sigmaE(double newSigE);
  
 private:
  int Z;
  int N;
  int A;
  double mA;
  double mAmpp;
  double mAmpn;
  double mAmnn;
  double mAm1p;   // A-1 residual mass, proton removed (single-nucleon MF)
  double mAm1n;   // A-1 residual mass, neutron removed (single-nucleon MF)
  double sigmaCM;
  double d_sigmaCM;
  double Estar;
  double Estar_max;
  double sigmaE;
  bool random_Estar;

  gcfSRC * mySRC;
  
};

#endif
