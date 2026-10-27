#ifndef __ION_GENERATOR_H__
#define __ION_GENERATOR_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/NNCrossSection.hh"

class ionGenerator: public gcfGenerator
{
 public:
  ionGenerator(double E, gcfNucleus * thisInfo, NNCrossSection * thisCS, TRandom3 * thisRand);
  ~ionGenerator();
  void generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector &v3_target, TLorentzVector &v4_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target);
  double get_s_nominal();
  void print_pp_summary();
  
 private:
  NNCrossSection * myCS;
    
  double Ebeam;
  double E1;
  double p1;
  TVector3 v1;
  TLorentzVector v1_target;

  // Events and weight for which the pp s lies below or above the pp data
  long nCalls_pp, nLow_pp, nHigh_pp;
  double sumW_pp, sumWLow_pp, sumWHigh_pp;
  
};

#endif
