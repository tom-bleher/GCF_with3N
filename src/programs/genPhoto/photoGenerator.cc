#include <iostream>
#include "TVector3.h"
#include "photoGenerator.hh"
#include "spectra/diamondSpectrum.hh"
#include "spectra/amorphousSpectrum.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

photoGenerator::photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{
  
  myCS = thisCS;
  myReaction = proton_piMinus;
  
  usingfixedE=false;
  fixedE=0;
  doLC = false;
  fillDiamond();

}

photoGenerator::photoGenerator(gcfNucleus * thisInfo, photoCrossSection * thisCS, TRandom3 * thisRand, spectrum thisSpectrum, reaction thisReaction) : gcfGenerator(thisInfo, thisRand)
{
  
  myCS = thisCS;
  myReaction = thisReaction;

  usingfixedE=false;
  fixedE=0;
  doLC = false;
  if (thisSpectrum == diamond)
    fillDiamond();
  else if (thisSpectrum == amorphous)
    fillAmorphous();

}

photoGenerator::~photoGenerator()
{
  delete photonSpectrum;
}

void photoGenerator::generate_event(double &weight, double &Ephoton, int &meson_type, double &mMeson, int &baryon_type, double &mBaryon, int &rec_type, TLorentzVector &vMeson_target, TLorentzVector &vBaryon_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  if (usingfixedE) {

    Ephoton=fixedE;
  }
  else{Ephoton = photonSpectrum->GetRandom(myRand);

  }	
  TLorentzVector vphoton_target(0.,0.,Ephoton,Ephoton);
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what branching ratio to generate.
  int lead_type;
  double mMesonMean;
  double gammaMeson;
  double mBaryonMean;
  double gammaBaryon;
  double mM_min;
  double mM_max;
  double mB_min=0;
  double mB_max=5;
  if (myReaction==proton_piMinus)
    {
      lead_type = nCode;
      meson_type = pimCode;
      baryon_type = pCode;
      mMesonMean = mpip;
      gammaMeson = 0.;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = 0;
      mM_max = 1.5;
    }
  else if (myReaction==proton_rho0)
    {
      lead_type = pCode; 
      meson_type = rho0Code;
      baryon_type = pCode;
      mMesonMean = mrho0;
      gammaMeson = gammarho0;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = 2*mpip;
      mM_max = 1.5;
    }
  else if (myReaction==proton_rhoMinus)
    {
      lead_type = nCode; 
      meson_type = rhomCode;
      baryon_type = pCode;
      mMesonMean = mrhop;
      gammaMeson = gammarhop;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = mpip + mpi0;
      mM_max = 1.5;
    }
  else if (myReaction==proton_omega)
    {
      lead_type = pCode; 
      meson_type = omegaCode;
      baryon_type = pCode;
      mMesonMean = momega;
      gammaMeson = gammaomega;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = 2*mpip + mpi0;
      mM_max = 1.5;
    }
  else if (myReaction==proton_phi)
    {
      lead_type = pCode; 
      meson_type = phiCode;
      baryon_type = pCode;
      mMesonMean = mphi;
      gammaMeson = gammaphi;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = 2*mKp;
      mM_max = 1.5;
    }
  else if (myReaction==neutron_phi)
    {
      lead_type = nCode; 
      meson_type = phiCode;
      baryon_type = nCode;
      mMesonMean = mphi;
      gammaMeson = gammaphi;
      mBaryonMean = mN;
      gammaBaryon = 0.;
      mM_min = 2*mKp;
      mM_max = 1.5;
    }
  else if (myReaction==DeltaPlusPlus_piMinus)
    {
      lead_type = pCode; 
      meson_type = pimCode;
      baryon_type = DeltappCode;
      mMesonMean = mpip;
      gammaMeson = 0;
      mBaryonMean = mDelta;
      gammaBaryon = gammaDelta;
      mM_min = 0;
      mM_max = 1.5;
      mB_min=mpip+mN;
    }
  else if (myReaction==DeltaPlus_piMinus)
    {
      lead_type = nCode; 
      meson_type = pimCode;
      baryon_type = DeltapCode;
      mMesonMean = mpip;
      gammaMeson = 0;
      mBaryonMean = mDelta;
      gammaBaryon = gammaDelta;
      mM_min = 0;
      mM_max = 1.5;
      mB_min=mpip+mN;
    }

  do
    mMeson = myRand->BreitWigner(mMesonMean,gammaMeson);
  while
    ((mMeson < mM_min) or (mMeson > mM_max));
  do
    mBaryon = myRand->BreitWigner(mBaryonMean,gammaBaryon);
  while
    ((mBaryon < mB_min) or (mBaryon > mB_max));
  
  // Decide recoil nucleon type
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 2.;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  TLorentzVector v1_target;
  if (!doLC)
    {
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
      v1_target.SetVect(v1);
      v1_target.SetT(E1);
    }
  else
    {
      TVector2 v1_perp, vRec_perp;
      double alpha1, alphaRec;
      decay_function_lc(weight, lead_type, rec_type, alpha1, v1_perp, alphaRec, vRec_perp);
      double alphaAm2 = Anum - alpha1 - alphaRec;
      TVector2 vAm2_perp = -1.*v1_perp - vRec_perp;

      double p1_minus = mbar*alpha1;
      double pRec_minus = mbar*alphaRec;
      double pAm2_minus = mbar*alphaAm2;
      double pRec_plus = (sq(mN) + vRec_perp.Mod2())/pRec_minus;
      double pAm2_plus = (sq(mAm2) + vAm2_perp.Mod2())/pAm2_minus;
      if (mAm2 == 0)
	pAm2_plus = 0.;
      double p1_plus = mA - pRec_plus - pAm2_plus;
      
      double E1 = 0.5*(p1_plus + p1_minus);
      double p1_z = 0.5*(p1_plus - p1_minus);
      double ERec = 0.5*(pRec_plus + pRec_minus);
      double pRec_z = 0.5*(pRec_plus - pRec_minus);
      double EAm2 = 0.5*(pAm2_plus + pAm2_minus);
      double pAm2_z = 0.5*(pAm2_plus - pAm2_minus);

      v1_target.SetXYZT(v1_perp.X(),v1_perp.Y(),p1_z,E1);
      vRec_target.SetXYZT(vRec_perp.X(),vRec_perp.Y(),pRec_z,ERec);
      vAm2_target.SetXYZT(vAm2_perp.X(),vAm2_perp.Y(),pAm2_z,EAm2);

      (v1_target+vRec_target+vAm2_target).Print();
      
    }

  double cosThetaCM;
  t_scatter(weight, mMeson, mBaryon, vphoton_target, v1_target, vMeson_target, vBaryon_target, cosThetaCM);
  
  if (weight <= 0.)
    return;
  
  double s = sq(mMeson) + sq(mBaryon) + 2.*vMeson_target.Dot(vBaryon_target);
  double t = sq(mMeson) - 2.*vMeson_target.Dot(vphoton_target);
  
  // Calculate the flux factor on the cross section
  double vgamma1 = vphoton_target.Dot(v1_target)/(Ephoton*v1_target.T());
  
  // Calculate the weight
  double thisCS=0;
  if (myReaction==proton_piMinus)
    thisCS=myCS->sigma_pim_p(s,cosThetaCM);
  else if (myReaction==proton_rho0)
    thisCS=myCS->sigma_rho0_p(s,t,cosThetaCM);
  else if (myReaction==proton_rhoMinus)
    thisCS=myCS->sigma_rhom_p(s,t,cosThetaCM);
  else if (myReaction==proton_omega)
    thisCS=myCS->sigma_omega_p(s,t,cosThetaCM);
  else if (myReaction==proton_phi)
    thisCS=myCS->sigma_phi_p(s,t,cosThetaCM);
  else if (myReaction==neutron_phi)
    thisCS=myCS->sigma_phi_n(s,t,cosThetaCM);
  else if (myReaction==DeltaPlusPlus_piMinus)
    thisCS=myCS->sigma_deltapp_pim(s,cosThetaCM);
  else if (myReaction==DeltaPlus_piMinus)
    thisCS=myCS->sigma_deltap_pim(s,cosThetaCM);

  weight *= vgamma1*thisCS; // Photoproduction cross section
  
}

void photoGenerator::setfixedE(double newfixedE){
  usingfixedE=true;
  fixedE=newfixedE;

}

void photoGenerator::print_beam_info()
{
  if (usingfixedE)
    {
      std::cerr << "photoGenerator: using a fixed beam energy of " << fixedE << " GeV.\n";
    }
  else
    {
      std::cerr << "photoGenerator: using the Hall-D photon source spectrum\n"
		<< "     " << 100.*photonSpectrum->Integral(120,160) / photonSpectrum->Integral() << " % falls within the coherent peak (8--9 GeV)\n";
    }
}

void photoGenerator::setDoLC(bool newDoLC)
{

  doLC = newDoLC;
  
}

void photoGenerator::fillDiamond()
{
  photonSpectrum = new TH1D("photonEnergy","photonEnergy",280,5.,12.);
  for (int i = 0; i < 280; i++)
    {
      photonSpectrum->SetBinContent(i+1,diamondSpectrum[i]);
    }
}

void photoGenerator::fillAmorphous()
{
  photonSpectrum = new TH1D("photonEnergy","photonEnergy",1000,2.,12.);
  for (int i = 0; i < 1000; i++)
    {
      photonSpectrum->SetBinContent(i+1,amorphousSpectrum[i]);
    }
}
