#include <iostream>
#include <fstream>
#include "TVector2.h"
#include "gcfGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

gcfGenerator::gcfGenerator(gcfNucleus * thisInfo, TRandom3 * thisRand)
{

  myInfo = thisInfo;
  myRand = thisRand;
  
  Anum = myInfo->get_A();
  mA = myInfo->get_mA();
  mbar = myInfo->get_mbar();
  mAmpp = myInfo->get_mAmpp(); // this includes the effect of Estar
  mAmpn = myInfo->get_mAmpn();
  mAmnn = myInfo->get_mAmnn();
  sigCM = myInfo->get_sigmaCM();
  random_Estar = myInfo->get_Estar_randomization();
  doRad = false;
  doCoul = false;

  pRelmin = 0.2;
  pRelmax = 1.05;
  phiRelmin = 0.;
  phiRelmax = 2.*M_PI;
  thetaRelmin = 0.;
  thetaRelmax = M_PI;  
  cosThetaRelmin = -1.;
  cosThetaRelmax = 1.;
  alphaRelmin = 0.;
  alphaRelmax = 2.;

  numin = 0.;
  numax = 0.;
  xBmin = 0.;
  xBmax = 0.;
  QSqmin = 0.;
  QSqmax = 0.;
  phikmin = 0.;
  phikmax = 0.;

}

gcfGenerator::~gcfGenerator()
{
}

void gcfGenerator::set_doRad(bool rad)
{
  doRad = rad;
}

void gcfGenerator::set_doCoul(bool coul)
{
  doCoul = coul;
}

void gcfGenerator::set_phiRel_range(double low, double high)
{
  phiRelmin = low;
  phiRelmax = high;
}

void gcfGenerator::set_phiRel_range_deg(double low, double high)
{
  phiRelmin = low*M_PI/180.;
  phiRelmax = high*M_PI/180.;
}

void gcfGenerator::set_thetaRel_range(double low, double high)
{
  thetaRelmin = low;
  thetaRelmax = high;
  cosThetaRelmin = cos(thetaRelmax);
  cosThetaRelmax = cos(thetaRelmin);
}

void gcfGenerator::set_thetaRel_range_deg(double low, double high)
{
  thetaRelmin = low*M_PI/180.;
  thetaRelmax = high*M_PI/180.;
  cosThetaRelmin = cos(thetaRelmax);
  cosThetaRelmax = cos(thetaRelmin);
}

void gcfGenerator::set_pRel_range(double low, double high)
{
  pRelmin = low;
  pRelmax = high;
}

void gcfGenerator::set_pRel_cut(double new_cutoff)
{
  pRel_cut = new_cutoff;
}

void gcfGenerator::set_nu_range(double low, double high)
{
  numin = low;
  numax = high;
}

void gcfGenerator::set_xB_range(double low, double high)
{
  xBmin = low;
  xBmax = high;
}

void gcfGenerator::set_QSq_range(double low, double high)
{
  QSqmin = low;
  QSqmax = high;
}

void gcfGenerator::set_phik_range(double low, double high)
{
  phikmin = low;
  phikmax = high;
}

void gcfGenerator::set_phik_range_deg(double low, double high)
{
  phikmin = low*M_PI/180.;
  phikmax = high*M_PI/180.;
}

void gcfGenerator::set_deuteron()
{
  pRelmin = 0.;
  pRel_cut = 0.;
}

void gcfGenerator::randomize_cutoff()
{
  pRel_cut = pRel_cut + (myRand->Uniform() - 0.5)*pRel_cut_range;
}

bool gcfGenerator::parse_phase_space_file(char* phase_space)
{

  ifstream ps_file(phase_space);
  string param;
  double low, high;
  while (ps_file >> param >> low >> high)
    {
      if (param == "phiRel" || param == "phirel")
	{
	  set_phiRel_range(low,high);
	}
      else if (param == "phiRel_deg" || param == "phirel_deg")
	{
	  set_phiRel_range_deg(low,high);
	}
      else if (param == "thetaRel" || param == "thetarel")
	{
	  set_thetaRel_range(low,high);
	}
      else if (param == "thetaRel_deg" || param == "thetarel_deg")
	{
	  set_thetaRel_range_deg(low,high);
	}
      else if (param == "pRel" || param == "prel")
	{
	  set_pRel_range(low,high);
	}
      else if (param == "nu" || param == "omega")
	{
	  set_nu_range(low,high);
	}
      else if (param == "x" || param == "xB")
	{
	  set_xB_range(low,high);
	}
      else if (param == "QSq" || param == "Qsq")
	{
	  set_QSq_range(low,high);
	}
      else if (param == "phik" || param == "phie")
	{
	  set_phik_range(low,high);
	}
      else if (param == "phik_deg" || param == "phie_deg")
	{
	  set_phik_range_deg(low,high);
	}
      else
	{
	  cerr << "Invalid phase space parameter provided. Aborting...\n";
	  return false;
	}
    }
  
  return true;
  
}


// ---------------------------------------------------------------------------
// |p_rel| importance sampling
//
// Uniform sampling of |p_rel| leaves the whole k^2 S(k) factor in the weight,
// and S(k) spans several decades, so a few low-k events dominate every sum.
// Sampling k from its own integrand instead makes that dimension contribute a
// constant (the integral over the allowed range) to the weight.
// ---------------------------------------------------------------------------

// 0 = pp, 1 = nn, 2 = pn — mirrors the branching inside gcfSRC::get_S.
int gcfGenerator::is_class(int lead_type, int rec_type) const
{
  if (lead_type == rec_type) return (lead_type == pCode) ? 0 : 1;
  return 2;
}

void gcfGenerator::build_pRel_IS_tables()
{
  fISdk = fISkmax / kISBins;
  const int types[3][2] = {{pCode, pCode}, {nCode, nCode}, {pCode, nCode}};
  for (int c = 0; c < 3; c++) {
    fISF[c].assign(kISBins + 1, 0.);
    double acc = 0., f_prev = 0.;
    for (int i = 1; i <= kISBins; i++) {
      const double k = i * fISdk;
      const double f = k * k * myInfo->get_S(k, types[c][0], types[c][1]);
      acc += 0.5 * (f + f_prev) * fISdk;   // trapezoid
      f_prev = f;
      fISF[c][i] = acc;
    }
  }
  fPRelIS_built = true;
}

double gcfGenerator::IS_F(int cls, double k) const
{
  if (k <= 0.) return 0.;
  if (k >= fISkmax) return fISF[cls][kISBins];
  const double x = k / fISdk;
  const int i = (int)x;
  const double t = x - i;
  return fISF[cls][i] + t * (fISF[cls][i + 1] - fISF[cls][i]);
}

double gcfGenerator::pRel_integral(int lead_type, int rec_type,
                                   double lo, double hi)
{
  if (!fPRelIS_built) build_pRel_IS_tables();
  const int cls = is_class(lead_type, rec_type);
  const double I = IS_F(cls, hi) - IS_F(cls, lo);
  return (I > 0.) ? I : 0.;
}

double gcfGenerator::IS_invF(int cls, double u, double &f_hat) const
{
  const std::vector<double> &F = fISF[cls];
  // binary search for the bin with F[i] <= u < F[i+1]
  int lo = 0, hi = kISBins;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (F[mid] <= u) lo = mid; else hi = mid;
  }
  const double dF = F[lo + 1] - F[lo];
  f_hat = dF / fISdk;                       // density implied by the interpolant
  if (dF <= 0.) return lo * fISdk;
  return (lo + (u - F[lo]) / dF) * fISdk;
}

void gcfGenerator::decay_function(double &weight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec)
{
  
  // Pick random CM motion
  TVector3 vCM(myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM));
  
  // Pick random relative motion
  TVector3 vRel;
  double phiRel = phiRelmin + (phiRelmax-phiRelmin)*myRand->Rndm();
  double cosThetaRel = cosThetaRelmin + (cosThetaRelmax - cosThetaRelmin)*myRand->Rndm();
  double thetaRel = acos(cosThetaRel);

  if ((phiRelmax < phiRelmin) or (cosThetaRelmax < cosThetaRelmin) or (pRelmax < pRelmin))
    {
      weight = 0.;
      return;
    }

  const double angular = (phiRelmax-phiRelmin) * (cosThetaRelmax - cosThetaRelmin);

  double pRel_Mag;
  if (fPRelIS)
    {
      // Draw |p_rel| from k^2 S(k) restricted to [pRelmin, pRelmax]. The
      // weight then picks up the integral over that range -- constant per
      // event -- times a residual f/f_hat that makes the estimator exactly
      // unbiased against the piecewise-linear cumulative.
      if (!fPRelIS_built) build_pRel_IS_tables();
      const int cls = is_class(lead_type, rec_type);
      const double Flo = IS_F(cls, pRelmin);
      const double Fhi = IS_F(cls, pRelmax);
      const double dF  = Fhi - Flo;
      if (dF <= 0.) { weight = 0.; return; }

      double f_hat = 0.;
      pRel_Mag = IS_invF(cls, Flo + dF*myRand->Rndm(), f_hat);
      if (f_hat <= 0.) { weight = 0.; return; }

      const double f = sq(pRel_Mag)*myInfo->get_S(pRel_Mag,lead_type,rec_type);
      weight *= angular * dF * (f / f_hat) / pow(2.*M_PI,3);
    }
  else
    {
      pRel_Mag = pRelmin + (pRelmax - pRelmin)*myRand->Rndm();
      // Factor universal functions into the weights
      weight *= angular * (pRelmax - pRelmin) // Phase Space
        * sq(pRel_Mag)*myInfo->get_S(pRel_Mag,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts
    }

  vRel.SetMagThetaPhi(pRel_Mag,thetaRel,phiRel);

  // Do a safeguard cut
  if (pRel_Mag < pRel_cut)
    {
      weight=0.;
      return;
    }

  // Determine initial nucleon momenta
  vi = 0.5*vCM + vRel;
  vRec = 0.5*vCM - vRel;

}

void gcfGenerator::decay_function_lc(double &weight, int lead_type, int rec_type, double &alphai, TVector2 &vi_perp, double &alphaRec, TVector2 &vRec_perp)
{
  
  // Pick random CM motion
  double alphaCM = myRand->Gaus(2.,sigCM/mbar);
  TVector2 vCM_perp(myRand->Gaus(0.,sigCM),myRand->Gaus(0.,sigCM));
  if (alphaCM < 0.)
    {
      weight=0.;
      return;
    }
  
  // Pick random relative motion
  double alphaRel = alphaRelmin + (alphaRelmax-alphaRelmin)*myRand->Rndm();
  double kmin = abs(1.-alphaRel)/sqrt(alphaRel*(2.-alphaRel))*mN;
  double k = max(pRelmin,kmin) + (pRelmax - max(pRelmin,kmin))*myRand->Rndm();
  double phiRel = phiRelmin + (phiRelmax-phiRelmin)*myRand->Rndm();
  double kperpSq = alphaRel*(2.-alphaRel)*(sq(k)+sq(mN)) - sq(mN);
  if (kperpSq < 0.)
    {
      weight=0.;
      return;
    }
  TVector2 vRel_perp;
  vRel_perp.SetMagPhi(sqrt(kperpSq),phiRel);

  if ((alphaRelmax < alphaRelmin) or (pRelmax < max(pRelmin,kmin)) or (phiRelmax < phiRelmin))
    {
      weight = 0.;
      return;
    }

  // Factor universal functions into the weights
  weight *= (alphaRelmax - alphaRelmin) * (pRelmax - max(pRelmin,kmin)) * (phiRelmax-phiRelmin) // Phase Space
    *k*sqrt(sq(k)+sq(mN))*myInfo->get_S(k,lead_type,rec_type)/pow(2.*M_PI,3); // Contacts
  
  // Do a safeguard cut
  if (k < pRel_cut)
    {
      weight=0.;
      return;
    }

  // Determine initial nucleon momenta
  alphai = alphaRel*alphaCM/2.;
  vi_perp = alphai/alphaCM*vCM_perp + vRel_perp;
  alphaRec = alphaCM - alphai;
  vRec_perp = alphaRec/alphaCM*vCM_perp - vRel_perp;
  
}

double gcfGenerator::get_mAm2(int lead_type, int rec_type)
{
  double Estar;
  return get_mAm2(lead_type,rec_type,Estar);
}

double gcfGenerator::get_mAm1(int lead_type)
{
  // Single-nucleon knockout: removing the struck nucleon leaves the A-1
  // residual. Proton struck -> A-1 with one fewer proton, and vice versa.
  return (lead_type == pCode) ? myInfo->get_mAm1p() : myInfo->get_mAm1n();
}

double gcfGenerator::get_mAm2(int lead_type, int rec_type, double &Estar)
{
  if (random_Estar)
    {
      if (lead_type == pCode and rec_type == pCode)
	return myInfo->get_mAmpp_random(Estar,myRand);
      else if (lead_type == nCode and rec_type == nCode)
	return myInfo->get_mAmnn_random(Estar,myRand);
      else
	return myInfo->get_mAmpn_random(Estar,myRand);
    }
  else
    {
      if (lead_type == pCode and rec_type == pCode)
	return mAmpp;
      else if (lead_type == nCode and rec_type == nCode)
	return mAmnn;
      else
	return mAmpn;
    }
}

void gcfGenerator::t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4)
{
  double cosThetaCM;
  t_scatter(weight, m3, m4, v1, v2, v3, v4, cosThetaCM);
}

void gcfGenerator::t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4, double &cosThetaCM)
{

  TLorentzVector Z_lab = v1 + v2;
  double s = Z_lab.M2();
  if (s < sq(m3 + m4))
    {
      weight=0.;
      return;
    }
  double W_cm = sqrt(s);
    
  // Boost vectors to scattering CM frame
  TVector3 m = Z_lab.BoostVector();
  TLorentzVector v1_cm = v1;
  TLorentzVector v2_cm = v2;
  v1_cm.Boost(-m);
  v2_cm.Boost(-m);

  // Rotate to scattering along z-axis
  double rot_phi = v1_cm.Vect().Phi();
  double rot_theta = v1_cm.Vect().Theta();

  // Determine scattered energy and momentum
  double E3_cm = (s + sq(m3) - sq(m4))/(2.*W_cm);
  double E4_cm = W_cm - E3_cm;
  double p_cm = sqrt(sq(E3_cm) - sq(m3));

  // Define Jacobian between cm scattering angle and t
  double J = 2.*p_cm*v1_cm.Vect().Mag()/(2.*M_PI);

  // Pick random CM scattering angle
  double phi_cm = 2.*M_PI*myRand->Rndm();
  cosThetaCM = -1. + 2.*myRand->Rndm();
  double theta_cm = acos(cosThetaCM);

  weight *= J * 4.*M_PI;

  // Set outgoing particle vectors
  TVector3 v_cm;
  v_cm.SetMagThetaPhi(p_cm,theta_cm,phi_cm);

  // Rotate back
  v_cm.RotateY(rot_theta);
  v_cm.RotateZ(rot_phi);
  TLorentzVector v3_cm(v_cm,E3_cm);
  TLorentzVector v4_cm(-v_cm,E4_cm);
    
  // Boost to lab system
  v3 = v3_cm;
  v4 = v4_cm;
  v3.Boost(m);
  v4.Boost(m);
    
}

TVector3 gcfGenerator::radiateElectron(TVector3 ve)
{
  double Ee = ve.Mag();
  double lambda_e = alpha/M_PI*(log(4*sq(Ee)/sq(me)) -1);
  double DeltaE = pow(myRand->Rndm(),1./lambda_e)*Ee;
  double Ee_rad = Ee - DeltaE;
  TVector3 ve_rad = ve;
  ve_rad.SetMag(Ee_rad);
  return ve_rad;
}

double gcfGenerator::deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/sq(me)) + 8./3.);
}

double gcfGenerator::radiationFactor(double Ebeam, double Ek, double QSq)
{
  double lambda_ei = alpha/M_PI*(log(4*sq(Ebeam)/sq(me)) -1);
  double lambda_ef = alpha/M_PI*(log(4*sq(Ek)/sq(me)) -1);

  return (1 - deltaHard(QSq)) * pow(Ebeam/sqrt(Ebeam*Ek),lambda_ei) * pow(Ek/sqrt(Ebeam*Ek),lambda_ef);
}

// True when the vertex kinematics (QSq, nu) lie inside the xB and nu windows
// of the phase-space file. A window left at its default (0, 0) is not applied.
// Like the QSq range, the windows act on the vertex (pre-radiation) electron.
bool gcfGenerator::in_xB_nu_window(double QSq, double nu) const
{
  if ((numin != 0. || numax != 0.) && (nu < numin || nu > numax))
    return false;
  if (xBmin != 0. || xBmax != 0.)
    {
      if (nu <= 0.)
        return false;
      double xB = QSq/(2.*mN*nu);
      if (xB < xBmin || xB > xBmax)
        return false;
    }
  return true;
}

double gcfGenerator::calcCoulombEnergy() {
  if (myInfo->get_Z() == 1 && myInfo->get_N() == 1) {
    return 0;
  }

  double radius = ((1.1 * cbrt(Anum)) + (0.86 / cbrt(Anum))); // fm
  double deltaE = 0.775 * ((3 * myInfo->get_Z() * alpha * GeVfm) / ( 2 * radius ));

  return deltaE;
}

void gcfGenerator::coulombCorrection(TLorentzVector &p, double deltaE) 
{
  // If p is 0, don't correct
  if(p.E() == 0) {
    return;
  }

  TVector3 p3 = p.Vect();
  double pMag = p3.Mag();

  if(deltaE < 0 && pMag < abs(deltaE)) {
    p.SetPxPyPzE(0, 0, 0, 0);
    return;
  }

  double deltaP = - pMag + sqrt(pow(pMag, 2) + 2 * p.E() * deltaE + pow(deltaE, 2));
  p3.SetMag(pMag + deltaP);

  p.SetPxPyPzE(p3.X(), p3.Y(), p3.Z(), p.E() + deltaE);
}
