#include "electroGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

electroGenerator::electroGenerator(double E, gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  numin = 1.;
  numax = 20.;
  QSqmin = 1.0;
  QSqmax = 5.0;
  phikmin = 0.;
  phikmax = 2.*M_PI;
  
}

electroGenerator::~electroGenerator()
{
}

void electroGenerator::generate_event(double &weight, int &meson_type, int &baryon_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &r)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  int lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Decide what branching ratio to generate (Currently only implemented gamma p -> J/psi p and equivalent neutron channel)
  // The J/psi is neutral, so the outgoing baryon is the same species as the struck nucleon.
  double mMeson;
  double mBaryon;
  if (lead_type == pCode)
    {
      meson_type = JpsiCode;
      baryon_type = pCode;
    }
  else
    {
      meson_type = JpsiCode;
      baryon_type = nCode;
    }
  mMeson = mJpsi;
  mBaryon = mN;
  weight *= 1;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  // Pick random electron scattering
  double nu = numin + (numax - numin)*myRand->Rndm();
  double Ek = Ebeam - nu;
  double QSqmax_kine = 4*Ebeam*Ek;
  
  if (QSqmax_kine < QSqmin)
    {
      weight=0.;
      return;
    }
  
  double QSq = QSqmin + (min(QSqmax,QSqmax_kine) - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (numax - numin) * (min(QSqmax,QSqmax_kine) - QSqmin) * (phikmax - phikmin);
  
  // Determine scattered electron kinematics
  double qSq = QSq + sq(nu);
  double costhetak = 1. - QSq/(2.*Ebeam*Ek);
  double thetak = acos(costhetak);
  double epsilon = 1./(1. + 2.*qSq/QSq*sq(tan(thetak/2.)));
  TVector3 vk;
  vk.SetMagThetaPhi(Ek,thetak,phik);
  vk_target.SetVect(vk);
  vk_target.SetT(vk.Mag());
  
  TLorentzVector vq_target = vbeam_target - vk_target;
  r = (epsilon*myCS->R_Jpsi_p(QSq))/(1. + epsilon*myCS->R_Jpsi_p(QSq));
  
  // Factor in photon flux
  weight *= 1/(2.*M_PI) * photon_flux(nu,QSq);
  
  TVector3 v1, vRec;
  decay_function(weight, lead_type, rec_type, v1, vRec);
  
  if (weight <= 0.)
    return;
  
  TVector3 vAm2 = - v1 - vRec;
  double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);
  
  double Erec = sqrt(sq(mN) + vRec.Mag2());  
  vRec_target.SetVect(vRec);
  vRec_target.SetT(Erec);
  
  double E1 = mA - EAm2 - Erec;
  TLorentzVector v1_target(v1,E1);
  
  t_scatter(weight, mMeson, mBaryon, vq_target, v1_target, vMeson_target, vBaryon_target);
  
  if (weight <= 0.)
    return;
  
  double s = sq(mMeson) + sq(mBaryon) + 2.*vMeson_target.Dot(vBaryon_target);
  double t = sq(mMeson) - QSq - 2.*vq_target.Dot(vMeson_target);
  
  // Calculate the flux factor on the cross section
  double vgamma1 = vq_target.Dot(v1_target)/(nu*E1);
  
  // Calculate the weight
  weight *= vgamma1*myCS->sigma_Jpsi_p(s,t,QSq); // Photoproduction cross section
  
}

//double electroGenerator::photon_flux(double nu, double QSq)
//{
//  double qSq = QSq + sq(nu);
//  double Ek = Ebeam - nu;
//  double costhetak = 1. - QSq/(2.*Ebeam*Ek);
//  double thetak = acos(costhetak);
//  double epsilon = 1./(1. + 2.*qSq/QSq*sq(tan(thetak/2.)));
//  double K_gamma = nu - QSq/(2*mA);
//  
//  return alpha/(2.*M_PI)*K_gamma/sq(Ebeam)*1./QSq*1./(1. - epsilon);
//}

double electroGenerator::photon_flux(double nu, double QSq)
{
  double y = nu/Ebeam;
  double QSq_min = sq(y*me)/(1-y);

  if (QSq < QSq_min)
    return 0;
  else
    return alpha/M_PI*1/(nu*QSq)*(1. - y + sq(y)/2. - (1.-y)*QSq_min/QSq);
}

