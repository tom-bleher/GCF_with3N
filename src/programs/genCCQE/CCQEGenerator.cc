#include "TVector3.h"
#include "CCQEGenerator.hh"
#include "spectra/uBSpectrum.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

CCQEGenerator::CCQEGenerator(gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;

  neutrinoSpectrum = new TH1D("photonEnergy","photonEnergy",200,0.,10.);
  for (int i = 0; i < 200; i++) {
    neutrinoSpectrum->SetBinContent(i+1,mu_spectrum[i]);
  }

  QSqmin = 1.e-4;
  QSqmax = 2.5;
  phikmin = 0.;
  phikmax = 2.*M_PI;
  
}

CCQEGenerator::~CCQEGenerator()
{
}


void CCQEGenerator::generate_event(double &weight, double &Eneutrino, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  
  Eneutrino = neutrinoSpectrum->GetRandom(myRand);
  
  TLorentzVector vneutrino_target(0.,0.,Eneutrino,Eneutrino);
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with (here fixed at np pairs)
  int lead_type = nCode;
  rec_type = pCode;

  double mk = mmu;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  // Sample decay function
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
  double p1_minus = E1 - v1.Z();

  // Pick random electron scattering
  double QSqmax_kine = 2*Eneutrino*p1_minus;
  if (QSqmax_kine < QSqmin)
    {
      weight=0.;
      return;
    }
  
  double QSq = QSqmin + (min(QSqmax,QSqmax_kine) - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (min(QSqmax,QSqmax_kine) - QSqmin) * (phikmax - phikmin);
  
  // Calculate outgoing electron kinematics
  double k_minus = (QSq + sq(mk))/(2*Eneutrino);
  double plead_minus = p1_minus - k_minus;
  if (plead_minus < 0.)
    {
      weight=0.;
      return;
    }
  
  double p1_plus = E1 + v1.Z();
  double virt = v1_target.Mag2() - sq(mN);
  
  double A = k_minus/p1_minus;
  double c = A*(p1_plus*k_minus - 2*Eneutrino*plead_minus - virt - sq(mk)) + sq(mk);
  
  double delta_phi = phik - v1.Phi();
  double p1_perp = v1.Perp();
  double y = p1_perp*cos(delta_phi);
  double b = -2*A*y;
  double D = sq(b) - 4*c;
  if (D < 0.)
    {
      weight=0.;
      return;
    }
  
  double k_perp = (- b + sqrt(D))/2.;

  // May be two possible solutions
  if (sqrt(D) < -b)
    {
      weight *= 2.;
      if (myRand->Rndm() > 0.5)
	k_perp = (- b - sqrt(D))/2.;
    }

  // Both roots negative: that is the phik+pi two-root case, already generated there
  if (k_perp < 0.)
    {
      weight=0.;
      return;
    }

  double k_plus = (sq(mk) + sq(k_perp))/k_minus;
  double kz = (k_plus - k_minus)/2.;
  TVector3 vk(k_perp*cos(phik),k_perp*sin(phik),kz);
  double Ek = (k_plus + k_minus)/2.;
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);

  vLead_target = v1_target + vneutrino_target - vk_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();
  
  // Jacobian for delta function
  double J = 2.*Eneutrino*vk.Mag()*fabs(p1_minus - k_minus*y/k_perp)/Elead;
  weight *= 1./J;

  // Calculate the weight
  weight *= myCS->sigma_CC(Eneutrino, vk, vLead, (lead_type==pCode)); // Charged-current cross section

  TLorentzVector vA = vk_target+vLead_target+vRec_target+vAm2_target-vneutrino_target;
  
}
