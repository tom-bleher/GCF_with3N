#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "helpers.hh"
#include <iostream>
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TFile.h"

using namespace std;

QEGenerator_3N::QEGenerator_3N(double E, eNCrossSection * thisCS, int thisU, TRandom3 * thisRand)
{

  myCS = thisCS;
  myRand = thisRand;
  u = thisU;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  sigCM = 0.15;  // 12C 2N SRC pair CM width, literature value (Cohen et al.)
  fkMinInternal = 0.25;  // GeV/c, lower cutoff on the internal (wavefunction-
                         // frame) momenta; settable to the nuclear kF
  fQ2min = 0.;
  fQ2max = 1.e9;
  // sigCM = 10.0; // sanity check (very high CM momentum -> theta12~theta23~0)
  phi_a_max = M_PI;
  phi_a_min =-M_PI;
  cos_theta_a_max = 1;
  cos_theta_a_min =-1;
  phi_baz_max = M_PI;
  phi_baz_min =-M_PI;
  mom_a_max = 5.0; // 1.2?
  mom_a_min = 0.00;
  mom_b_max = 5.0;
  mom_b_min = 0.00;
  theta_ab_max = M_PI;
  theta_ab_min = 0.0001;
  cos_theta_k_max = 1;
  cos_theta_k_min = sqrt(2.)/2.;
  phi_k_max = M_PI;
  phi_k_min =-M_PI;

  // Default target nucleus: 12C. Sets fTargetMass / fResidualMass.
  SetTargetNucleus(12, 6);

  fill_array();
}

QEGenerator_3N::~QEGenerator_3N()
{
}

// ---------------------------------------------------------------------------
// LookupNuclearMass: hard-coded table of nuclear masses keyed by (A, Z).
// Includes special cases for A=0 (no residual) and A=1 (single nucleon).
// Returns -1 if the nucleus is not known; callers should error out and
// instruct the user to add the mass to constants/nuclearMasses.hh.
double QEGenerator_3N::LookupNuclearMass(int A, int Z)
{
  if (A == 0 && Z == 0) return 0.0;          // no residual
  if (A == 1 && Z == 0) return mN;           // free neutron (mN ~ mn ~ mp)
  if (A == 1 && Z == 1) return m_1H;         // proton
  if (A == 2 && Z == 1) return m_2H;
  if (A == 3 && Z == 1) return m_3H;
  if (A == 3 && Z == 2) return m_3He;
  if (A == 4 && Z == 2) return m_4He;
  if (A == 6 && Z == 3) return m_6Li;
  if (A == 8 && Z == 4) return m_8Be;
  if (A == 9 && Z == 4) return m_9Be;
  if (A == 10 && Z == 4) return m_10Be;
  if (A == 10 && Z == 5) return m_10B;
  if (A == 11 && Z == 5) return m_11B;
  if (A == 10 && Z == 6) return m_10C;
  if (A == 12 && Z == 6) return m_12C;
  if (A == 14 && Z == 7) return m_14N;
  if (A == 16 && Z == 8) return m_16O;
  if (A == 25 && Z == 11) return m_25Na;
  if (A == 25 && Z == 12) return m_25Mg;
  if (A == 25 && Z == 13) return m_25Al;
  if (A == 27 && Z == 13) return m_27Al;
  if (A == 38 && Z == 16) return m_38S;
  if (A == 38 && Z == 17) return m_38Cl;
  if (A == 38 && Z == 18) return m_38Ar;
  if (A == 38 && Z == 19) return m_38K;
  if (A == 38 && Z == 20) return m_38Ca;
  if (A == 40 && Z == 18) return m_40Ar;
  if (A == 40 && Z == 20) return m_40Ca;
  if (A == 54 && Z == 24) return m_54Cr;
  if (A == 54 && Z == 25) return m_54Mn;
  if (A == 54 && Z == 26) return m_54Fe;
  if (A == 56 && Z == 26) return m_56Fe;
  if (A == 206 && Z == 80) return m_206Hg;
  if (A == 206 && Z == 81) return m_206Tl;
  if (A == 206 && Z == 82) return m_206Pb;
  if (A == 208 && Z == 82) return m_208Pb;
  return -1.0;
}

void QEGenerator_3N::SetTargetNucleus(int A, int Z)
{
  if (A < 3 || Z < 2 || (A - Z) < 1) {
    std::cerr << "QEGenerator_3N::SetTargetNucleus: cannot eject a ppn triplet "
              << "from (A=" << A << ", Z=" << Z << "); need A>=3, Z>=2, N>=1."
              << std::endl;
    std::exit(1);
  }
  const double mA_new = LookupNuclearMass(A, Z);
  if (mA_new < 0) {
    std::cerr << "QEGenerator_3N::SetTargetNucleus: mass for target (A="
              << A << ", Z=" << Z << ") not in the hard-coded table. "
              << "Add it to constants/nuclearMasses.hh and extend "
              << "QEGenerator_3N::LookupNuclearMass." << std::endl;
    std::exit(1);
  }
  const int A_res = A - 3;
  const int Z_res = Z - 2;
  const double mAm3_new = LookupNuclearMass(A_res, Z_res);
  if (mAm3_new < 0) {
    std::cerr << "QEGenerator_3N::SetTargetNucleus: mass for residual (A="
              << A_res << ", Z=" << Z_res << ") after ejecting a ppn triplet "
              << "from (A=" << A << ", Z=" << Z << ") not in the table. "
              << "Add it to constants/nuclearMasses.hh and extend "
              << "QEGenerator_3N::LookupNuclearMass." << std::endl;
    std::exit(1);
  }
  fTargetA      = A;
  fTargetZ      = Z;
  fTargetMass   = mA_new;
  fResidualMass = mAm3_new;
}

void QEGenerator_3N::set_theta_k_maxmin(double min, double max)
{
  double c1=cos(min*M_PI/180);
  double c2=cos(max*M_PI/180);
  
  if(c1>c2){
    cos_theta_k_max = c1;
    cos_theta_k_min = c2;  
  }
  else{
    cos_theta_k_max = c2;
    cos_theta_k_min = c1;  
  }

}

void QEGenerator_3N::generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& v_k_target, TLorentzVector &v_Lead_target, TLorentzVector &v_2_target, TLorentzVector &v_3_target, TLorentzVector &v_Am3_target, bool use_CM)
{
  double Estar=0;
  generate_event(weight,N1_type,N2_type,N3_type,v_k_target,v_Lead_target,v_2_target,v_3_target,v_Am3_target,use_CM, Estar);
}

void QEGenerator_3N::generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& v_k_target, TLorentzVector &v_Lead_target, TLorentzVector &v_2_target, TLorentzVector &v_3_target, TLorentzVector &v_Am3_target, bool use_CM, double &Estar)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  //We only want ppn triplets
  N1_type = pCode;
  N2_type = pCode;
  N3_type = pCode;
  double ntypes = myRand->Rndm();
  ((ntypes<0.333333)?N1_type:(ntypes<0.666667)?N2_type:N3_type)=nCode;
  weight*=3.;
  
  //Combinitorial factor
  if(N1_type==pCode){
    //weight*=2.;
  }

  // Target nucleus and residual mass are configured via SetTargetNucleus
  // (default 12C -> 9Be).
  double mA   = fTargetMass;
  double mAm3 = fResidualMass + Estar;

  //1st, sample the center of mass momentum
  
 TVector3 v_cm(0,0,0);
  if (use_CM) {v_cm.SetXYZ(myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM));}
  TVector3 v_cm_component = v_cm;
  if(v_cm.Mag() != 0) v_cm_component.SetMag(v_cm_component.Mag()/3.0);
  weight*=1.;

  //2nd, sample the Euler angles
  double phi_a = phi_a_min + (phi_a_max-phi_a_min)*myRand->Rndm();
  double cos_theta_a = cos_theta_a_min + (cos_theta_a_max-cos_theta_a_min)*myRand->Rndm();
  double theta_a = acos(cos_theta_a);
  double phi_baz = phi_baz_min + (phi_baz_max-phi_baz_min)*myRand->Rndm();
  weight*=(phi_a_max-phi_a_min)*(cos_theta_a_max-cos_theta_a_min)*(phi_baz_max-phi_baz_min);
  
  //3rd, sample the NNN variables
  double mom_a = mom_a_min + (mom_a_max - mom_a_min)*myRand->Rndm();
  double mom_b = mom_b_min + (mom_b_max - mom_b_min)*myRand->Rndm();
  double theta_ab = theta_ab_min + (theta_ab_max - theta_ab_min)*myRand->Rndm();
  weight*=(mom_a_max - mom_a_min)*(mom_b_max - mom_b_min)*(theta_ab_max - theta_ab_min);

  //Now define all of our vectors
  TVector3 v_a, v_b, v_1, v_2, v_3, v_Am3;
  //Start with va and vb
  v_a.SetXYZ(0,0,mom_a);
  v_b.SetMagThetaPhi(mom_b,theta_ab,phi_baz);
  //Rotate to get the correct angles
  v_a.RotateY(theta_a);
  v_b.RotateY(theta_a);
  v_a.RotateZ(phi_a);
  v_b.RotateZ(phi_a);
  //Now get v1, v2, v3
  v_1 = v_a;
  v_2 = (- 0.5 * v_a) +  v_b;
  v_3 = (- 0.5 * v_a) -  v_b;

  // Momentum cuts on the INTERNAL (wavefunction-frame) momenta, i.e. the
  // ones get_rho evaluates the 3N density at. Cutting after the CM boost
  // let configurations with tiny internal momenta through, where the
  // 1/(p1*p2*p3) Jacobian in get_rho_ptot_f1f2f3 diverges and produced
  // rare huge-weight events.
  if(v_1.Mag() < fkMinInternal || v_2.Mag() < fkMinInternal || v_3.Mag() < fkMinInternal){
    weight = 0;
    return;
  }
  if(v_1.Mag() > 5 || v_2.Mag() > 5 || v_3.Mag() > 5){
    weight = 0;
    return;
  }

  //Now add the center of mass motion to all
  v_1 += v_cm_component;
  v_2 += v_cm_component;
  v_3 += v_cm_component;
  v_Am3 = - v_cm;

  if (weight <= 0.)
    return;
  
  //Define the 4 vectors for all
  double E_Am3;
  if (fTargetA == 3 && fTargetZ == 2) {E_Am3 = 0;}  // 3He: no residual
  else {E_Am3 = sqrt(sq(mAm3) + v_Am3.Mag2());}
  
  v_Am3_target.SetVect(v_Am3);
  v_Am3_target.SetT(E_Am3);
  
  double E_2 = sqrt(sq(mN) + v_2.Mag2());  
  v_2_target.SetVect(v_2);
  v_2_target.SetT(E_2);

  double E_3 = sqrt(sq(mN) + v_3.Mag2());  
  v_3_target.SetVect(v_3);
  v_3_target.SetT(E_3);
 
  double E_1 = mA - E_Am3 - E_2 - E_3;
  TLorentzVector v_1_target(v_1,E_1);

  double phi_k = phi_k_min + (phi_k_max - phi_k_min)*myRand->Rndm();
  double cos_theta_k = cos_theta_k_min + (cos_theta_k_max-cos_theta_k_min)*myRand->Rndm();
  double theta_k = acos(cos_theta_k);

  weight *= (cos_theta_k_max - cos_theta_k_min) * (phi_k_max - phi_k_min);
  //Define some variables for solving E_k
  double Z = mA + Ebeam - E_2 - E_3 - E_Am3;
  TVector3 v_X = v_1 + vbeam;
  TVector3 hat_k;
  hat_k.SetMagThetaPhi(1,theta_k,phi_k);
  double E_k = (Z*Z - (mN*mN) - v_X.Mag2())/(2*Z - 2*v_X.Dot(hat_k));

  if (E_k <= 0.)
    weight=0;

  if (Z-E_k <= 0.)
    weight=0;

  if (weight <= 0.)
    return;

  
  //Outgoing electron, q, and outgoing lead
  TVector3 v_k = hat_k;
  v_k.SetMag(E_k);
  v_k_target.SetVect(v_k);
  v_k_target.SetT(E_k);
  TLorentzVector v_q_target = vbeam_target - v_k_target;

  const double Q2 = -v_q_target.M2();
  if (Q2 < fQ2min || Q2 > fQ2max) { weight = 0; return; }

  v_Lead_target = v_1_target + v_q_target;
  TVector3 v_Lead(v_Lead_target.X(),v_Lead_target.Y(),v_Lead_target.Z());
  double E_Lead = v_Lead_target.T();

  // Jacobian for the E_k delta function: |E'_1 / (E'_1 - p'_1 . k_hat)|.
  // Differs deliberately from Eq. 18 of 3N_Generator.pdf, which mis-applies the
  // chain rule (dots p'_1 into the full argument instead of the inner derivative -k_hat).
  double J_delta = E_Lead - v_Lead.Dot(hat_k);
  weight *= fabs(E_Lead/J_delta);

  //Add in the actual cross section of the weight
  //Fudge some numbers for C
  double C = 0.12;
  double t = 2;
  // cout << "theta_ab" << theta_ab * 180. / M_PI << endl;
  
  double sigma = myCS->sigma_eN(Ebeam, v_k, v_Lead, (N1_type==pCode));
  double rho_val = get_rho(N2_type, N3_type, theta_ab, mom_a, mom_b);
  double delta_jac = fabs(E_Lead/J_delta);

  weight *= 0.5 * pow(2 * M_PI,-7) * C * t * sigma;
  weight *= rho_val;
  // weight *= myCS->sigma_eN(Ebeam, v_k, v_Lead, (N1_type==pCode));
  if(weight<0){weight=0;}

  // Store components for the extended overload
  _last_rho = rho_val;
  _last_delta_jacobian = delta_jac;
}

void QEGenerator_3N::generate_event(double &weight, int &N1_type, int &N2_type, int &N3_type, TLorentzVector& v_k_target, TLorentzVector &v_Lead_target, TLorentzVector &v_2_target, TLorentzVector &v_3_target, TLorentzVector &v_Am3_target, bool use_CM, double &Estar, double &rho_out, double &delta_jacobian_out)
{
  _last_rho = 0;
  _last_delta_jacobian = 0;
  generate_event(weight, N1_type, N2_type, N3_type, v_k_target, v_Lead_target, v_2_target, v_3_target, v_Am3_target, use_CM, Estar);
  rho_out = _last_rho;
  delta_jacobian_out = _last_delta_jacobian;
}

double QEGenerator_3N::get_rho(double N2_type, double N3_type, double theta_ab, double p_a, double p_b){
  ////////////////////////////////
  //Define Some Vectors to get pn
  ////////////////////////////////
  //Now define all of our vectors
  TVector3 v_a, v_b, v_1, v_2, v_3;
  //Start with va and vb
  v_a.SetXYZ(0,0,p_a);
  v_b.SetMagThetaPhi(p_b,theta_ab,0);
  //Now get v1, v2, v3
  v_1 = v_a;
  v_2 = (- 0.5 * v_a) +  v_b;
  v_3 = (- 0.5 * v_a) -  v_b;
  //Now change them to get 23 -> 12
  TVector3 v_a_prime, v_b_prime, v_1_prime, v_2_prime, v_3_prime;
  if((N2_type==pCode) && (N3_type==pCode)){
    v_1_prime = v_1;
    v_2_prime = v_2;
    v_3_prime = v_3;
  }
  else if((N2_type==pCode) && (N3_type==nCode)){
    v_1_prime = v_3;
    v_2_prime = v_2;
    v_3_prime = v_1;
  }
  else if((N2_type==nCode) && (N3_type==pCode)){
    v_1_prime = v_2;
    v_2_prime = v_1;
    v_3_prime = v_3;
  }
  else{
    return 0;
  }
  v_a_prime = v_1_prime;
  v_b_prime = 0.5 * (v_2_prime - v_3_prime);
  double p_a_prime = v_a_prime.Mag();
  double p_b_prime = v_b_prime.Mag();
  double theta_ab_prime = v_a_prime.Angle(v_b_prime);
  
  //The 3N interaction is given in terms of pp-SRCs: the table's k_cm axis is
  //the neutron momentum and k_rel is the pp half-relative momentum, so the
  //relabeling above (v_a->v_a_prime etc.) is a pointwise re-evaluation of the
  //invariant density -- it needs NO Jacobian (the Cartesian map between the
  //two Jacobi sets has |det| = 1).
  //
  //The table stores the raw density |phi(pa',pb',theta')|^2 per unit
  //d^3pa d^3pb (nonzero at pa=pb=0; integrates to ~2 against the full
  //measure). Since the generator samples pa, pb, theta_ab uniformly, the
  //weight needs the phase-space measure of the SAMPLED variables:
  //d^3pa d^3pb = pa^2 pb^2 sin(theta_ab) dpa dpb dtheta_ab dOmega_a dphi_b.
  //Constant unit factors (hbar*c powers) are absorbed into the contact C.
  double J = p_a*p_a*p_b*p_b*sin(theta_ab);

  ////////////////////////////////
  //Define Some Vectors to get pn
  ////////////////////////////////    
  double theta = theta_ab_prime * 180 / M_PI;
  double k_cm = p_a_prime / GeVfm;
  double k_rel = p_b_prime / GeVfm;
   if((k_cm>10) || (k_rel>10)){
    return 0;
  }

  return J*interpolate(density_matrix_pp,theta,k_cm,k_rel);
}

double QEGenerator_3N::get_rho_ptot_f1f2f3(double k_1, double k_2, double k_3){
  //We want to get the wave function in variables of
  // f1, f2, and p_tot. To do this, we need a different
  //Jacobian. So that is the purpose of some of these
  //calculations. 
  
  // Check for zero or negative momenta
  if (k_1 <= 0 || k_2 <= 0 || k_3 <= 0) {
    return 0;
  }
  
  double p_1 = k_1 * GeVfm;
  double p_2 = k_2 * GeVfm;
  double p_3 = k_3 * GeVfm;
  
  // Check for valid triangle inequality and acos argument
  double acos_arg = (sq(k_1) + sq(k_2) - sq(k_3))/(2*k_1*k_2);
  if (acos_arg < -1.0 || acos_arg > 1.0) {
    return 0;  // Invalid triangle
  }
  
  double theta_12 = M_PI - acos(acos_arg);
  
  
  TVector3 v_k_1(0.0,0.0,k_1);
  TVector3 v_k_2;
  v_k_2.SetMagThetaPhi(k_2,theta_12,0.0);
  TVector3 v_k_3 = - v_k_1 - v_k_2;
  if(fabs(v_k_3.Mag()-k_3)>0.000001){
    cout<<"Geometry Problem!"<<endl;
  }
  double sin_theta_12 = sin(theta_12);
  
  TVector3 v_k_rel = (v_k_1 - v_k_2)*0.5; // p_b
  double k_rel = v_k_rel.Mag();
  TVector3 v_k_cm = -v_k_3; // p_a
  double k_cm = v_k_cm.Mag();
  double sin_theta_prime = sin(v_k_cm.Angle(v_k_rel)); // theta_ab
  //Units of Degrees
  double theta_prime = asin(sin_theta_prime)*180/M_PI;
  
  double p_rel = k_rel * GeVfm;
  double p_cm = k_cm * GeVfm;
  double p_tot = p_1 + p_2 + p_3;
  
  // Check for potential division by zero
  if (p_tot <= 0) {
    return 0;
  }
  
  //double J = (8/sqrt(3)) * (p_cm * p_cm * p_rel * p_rel * sin(theta_prime*M_PI/180)) / (p_1 * p_2 * p_3 * p_tot * p_tot);
  double J = (8/sqrt(3)) / (p_1 * p_2 * p_3 * p_tot * p_tot);
  
  if((k_cm>10) || (k_rel>10)){
    return 0;
  }
  else{
    double result = interpolate(density_matrix_pp,theta_prime,k_cm,k_rel)*J;
    return result;
  }
}


void QEGenerator_3N::fill_array(){
  if(u==1){uType="array/AV8_1D.dnst";}
  else if(u==2){uType="array/N2LO_1D.dnst";}
  else if(u==3){uType="array/G3_1D.dnst";}
  else if(u==4){uType="array/AV4_1D.dnst";}
  std::ifstream densityFile(uType);
  if(!densityFile.is_open()){
    cout<<"3 nucleon distribution file failed to load.\n"
	<<"Aborting...\n\n\n";
    exit(-2);
  }

  //densityFile.ignore(1000,'\n');  
  for(int i = 0; i < theta_bins; i++) {
    for(int j = 0; j < k_cm_bins; j++) {
      for(int k = 0; k < k_rel_bins; k++){
        std::string C_pp;

        densityFile >> C_pp;
        density_matrix_pp[i][j][k] = std::stod(C_pp);
        densityFile.ignore(1000,'\n');
      }
    }
  }
}

double QEGenerator_3N::interpolate(double matrix[k_rel_bins], double k_rel){
  //Get the "bin" for the value
  double float_bin = (k_rel - k_rel_min)/k_rel_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return matrix[0];
  }
  if(upper_bin>k_rel_bins-1){
    return matrix[k_rel_bins-1];
  }
  //Otherwise return interpolation
  return delta*matrix[upper_bin] + (1-delta)*matrix[lower_bin];
}

double QEGenerator_3N::interpolate(double matrix[k_cm_bins][k_rel_bins], double k_cm, double k_rel){
  //Get the "bin" for the value
  double float_bin = (k_cm - k_cm_min)/k_cm_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return interpolate(matrix[0],k_rel);
  }
  if(upper_bin>k_cm_bins-1){
    return interpolate(matrix[k_cm_bins-1],k_rel);
  }
  //Otherwise return interpolation
  return delta*interpolate(matrix[upper_bin],k_rel) + (1-delta)*interpolate(matrix[lower_bin],k_rel);
}

double QEGenerator_3N::interpolate(double matrix[theta_bins][k_cm_bins][k_rel_bins], double theta, double k_cm, double k_rel){
  //Get the "bin" for the value
  double float_bin = (theta - theta_min)/theta_width;
  int lower_bin = float_bin;
  int upper_bin = lower_bin+1;
  double delta = float_bin - (double)lower_bin;

  //Check to see if bins are out of bounds.
  //If so, return boundary
  if(lower_bin<0){
    return interpolate(matrix[0],k_cm,k_rel);
  }
  if(upper_bin>theta_bins-1){
    return interpolate(matrix[theta_bins-1],k_cm,k_rel);
  }  
  //Otherwise return interpolation
  return delta*interpolate(matrix[upper_bin],k_cm,k_rel) + (1-delta)*interpolate(matrix[lower_bin],k_cm,k_rel);
}
