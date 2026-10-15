#ifndef __QE_GENERATOR_3N_H__
#define __QE_GENERATOR_3N_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "generator/gcfGenerator.hh"
#include "nucleus/gcfNucleus.hh"
#include "cross_sections/eNCrossSection.hh"
#include <fstream>
#include <iostream>
#include <string>

class QEGenerator_3N
{
 public:
  QEGenerator_3N(double E, eNCrossSection * thisCS, int thisU, TRandom3 * thisRand);
  ~QEGenerator_3N();
  void generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &v2_target, TLorentzVector &v3_target, TLorentzVector &vAm3_target, bool use_CM);
  void generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &v2_target, TLorentzVector &v3_target, TLorentzVector &vAm3_target, bool use_CM, double &Estar);
  void generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &v2_target, TLorentzVector &v3_target, TLorentzVector &vAm3_target, bool use_CM, double &Estar, double &rho_out, double &delta_jacobian_out);
  void set_theta_k_maxmin(double min, double max);
  double get_rho(double N2_Type, double N3_Type, double theta_ab, double p_a, double p_b);
  double get_rho_ptot_f1f2f3(double k_1, double k_2, double k_3);
  void create_wavefunction_heatmap();

  // Target nucleus configuration (kinematic side). Validates that the
  // mass of (A, Z) and the residual (A-3, Z-2) after ejecting a ppn
  // triplet are both known to the lookup table; exits with an error
  // message if not. Default is 12C.
  void SetTargetNucleus(int A, int Z);
  int  GetTargetA() const { return fTargetA; }
  int  GetTargetZ() const { return fTargetZ; }

  // Per-component CM Gaussian width (GeV/c). Whether CM smearing is
  // applied at all is still controlled by the `use_CM` flag passed to
  // generate_event(). Default 0.15 (12C 2N SRC pair literature value).
  void SetSigCM(double sigCM_GeV) { sigCM = sigCM_GeV; }
  double GetSigCM() const { return sigCM; }

  // Lower cutoff (GeV/c) on the internal (wavefunction-frame) momenta of
  // the three nucleons, i.e. the SRC condition k > kF on the initial
  // state. Default 0.25; set to the target's Fermi momentum for a
  // nucleus-dependent cut.
  void SetInternalMomentumMin(double kmin_GeV) { fkMinInternal = kmin_GeV; }
  double GetInternalMomentumMin() const { return fkMinInternal; }

  // Q^2 window (GeV^2) applied at generation time. Q^2 is not a sampled
  // variable (the scattered-electron energy is fixed by energy-momentum
  // conservation given the hadronic configuration), so this is a
  // rejection: events outside the window get weight = 0 and are not
  // filled. Default is fully open.
  void SetQ2Range(double q2min, double q2max) { fQ2min = q2min; fQ2max = q2max; }

  // Mass lookup for (A, Z). Returns -1 if the nucleus is not in the
  // hard-coded table (see QEGenerator_3N.cc). Handles A=0 (no
  // residual), A=1 (single nucleon) as special cases.
  static double LookupNuclearMass(int A, int Z);
  
 private:
  eNCrossSection * myCS;
    
  double Ebeam;
  TVector3 vbeam;
  TLorentzVector vbeam_target;
  TRandom3 * myRand;

  int u;
  std::string uType;      // path of the 3N density table loaded by fill_array()

  double sigCM;
  double fkMinInternal;   // GeV/c, internal-momentum lower cutoff (default 0.25)
  double fQ2min;          // GeV^2, generation-time Q^2 window (default open)
  double fQ2max;
  int    fTargetA;        // mass number of host nucleus (default 12 = 12C)
  int    fTargetZ;        // charge number of host nucleus (default 6)
  double fTargetMass;     // = LookupNuclearMass(fTargetA, fTargetZ)
  double fResidualMass;   // = LookupNuclearMass(fTargetA-3, fTargetZ-2)
  double phi_a_max;
  double phi_a_min;
  double cos_theta_a_max;
  double cos_theta_a_min;
  double phi_baz_max;
  double phi_baz_min;
  double mom_a_max;
  double mom_a_min;
  double mom_b_max;
  double mom_b_min;
  double theta_ab_max;
  double theta_ab_min;
  double cos_theta_k_max;
  double cos_theta_k_min;
  double phi_k_max;
  double phi_k_min;

  //Functions to get the 3D wavefunction
  static const int theta_bins = 61; // 61
  const double theta_min = 0;
  const double theta_max = 180;
  const double theta_width = (theta_max - theta_min)/((double)theta_bins-1);

  static const int k_cm_bins = 101; // 101
  const double k_cm_min = 0;
  const double k_cm_max = 10;
  const double k_cm_width = (k_cm_max - k_cm_min)/((double)k_cm_bins-1);

  static const int k_rel_bins = 101; // 101
  const double k_rel_min = 0;
  const double k_rel_max = 10;
  const double k_rel_width = (k_rel_max - k_rel_min)/((double)k_rel_bins-1);

  double _last_rho = 0;
  double _last_delta_jacobian = 0;

  double density_matrix_pp[theta_bins][k_cm_bins][k_rel_bins];
  //double density_matrix_pn[theta_bins][k_cm_bins][k_rel_bins];

  void fill_array();
  double interpolate(double matrix[k_rel_bins], double k_rel);
  double interpolate(double matrix[k_cm_bins][k_rel_bins], double k_cm, double k_rel);
  double interpolate(double matrix[theta_bins][k_cm_bins][k_rel_bins], double theta, double k_cm, double k_rel);
  
};

#endif
