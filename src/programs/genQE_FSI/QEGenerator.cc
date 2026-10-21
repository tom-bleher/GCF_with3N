#include "QEGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

QEGenerator::QEGenerator(double E, gcfNucleus * thisInfo, eNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;

  Ebeam = E;
  vbeam.SetXYZ(0.,0.,Ebeam);
  vbeam_target.SetXYZT(0.,0.,Ebeam,Ebeam);

  QSqmin = 1.0;
  QSqmax = 10.0;
  phikmin = 0.;
  phikmax = 2.*M_PI;
  
}

QEGenerator::~QEGenerator()
{
}

void QEGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  double Estar;
  generate_event(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,Estar);
}

void QEGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &Estar)
{
  TVector3 vRel_target;
  TLorentzVector q_target;
  generate_event(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,vRel_target,q_target,Estar);
}

void QEGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, TVector3 &vRel_target, TLorentzVector &q_target, double &Estar)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type, Estar);

  // Sample decay function
  TVector3 v1, vRec;
  decay_function(weight, lead_type, rec_type, v1, vRec);
  
  if (weight <= 0.)
    return;
  
  TVector3 vAm2 = - v1 - vRec;
  double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);

  vRel_target = 0.5*(v1 - vRec);
  
  double Erec = sqrt(sq(mN) + vRec.Mag2());  
  vRec_target.SetVect(vRec);
  vRec_target.SetT(Erec);
  
  double E1 = mA - EAm2 - Erec;
  TLorentzVector v1_target(v1,E1);
  double p1_minus = E1 - v1.Z();

  TVector3 vbeam_int = vbeam;

  // Initial Coulomb effect
  if(doCoul) {
    double deltaECoul = calcCoulombEnergy();
    vbeam_int.SetMag(vbeam_int.Mag() + deltaECoul);
  }

  // Initial radiation
  vbeam_int = (doRad ? radiateElectron(vbeam_int) : vbeam_int);
  double Ebeam_int = vbeam_int.Mag();
  TLorentzVector vbeam_int_target(vbeam_int,Ebeam_int);

  // Pick random electron scattering
  double QSqmax_kine = 2*Ebeam_int*p1_minus;
  if (QSqmax_kine < QSqmin)
    {
      weight=0.;
      return;
    }
  
  double QSq = QSqmin + (min(QSqmax,QSqmax_kine) - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (min(QSqmax,QSqmax_kine) - QSqmin) * (phikmax - phikmin);
  
  // Calculate outgoing electron kinematics
  double k_minus = QSq/(2*Ebeam_int);
  double plead_minus = p1_minus - k_minus;
  if (plead_minus < 0.)
    {
      weight=0.;
      return;
    }
  
  double p1_plus = E1 + v1.Z();
  double virt = v1_target.Mag2() - sq(mN);
  
  double A = k_minus/p1_minus;
  double c = A*(p1_plus*k_minus - 2*Ebeam_int*plead_minus - virt);
  
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

  double k_plus = sq(k_perp)/k_minus;
  double kz = (k_plus - k_minus)/2.;
  TVector3 vk_int(k_perp*cos(phik),k_perp*sin(phik),kz);
  double Ek_int = vk_int.Mag();
  TLorentzVector vk_int_target;
  vk_int_target.SetVect(vk_int);
  vk_int_target.SetT(Ek_int);

  q_target = vbeam_int_target - vk_int_target;

  vLead_target = v1_target + vbeam_int_target - vk_int_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();

  // Final Radiation
  TVector3 vk = (doRad ? radiateElectron(vk_int) : vk_int);
  double Ek = vk.Mag();
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);
  
  // Jacobian for delta function
  double J = 2.*Ebeam_int*Ek_int*fabs(1. - (v1.Z() + y * tan(vk.Theta()/2.) + Ebeam_int - Ek_int)/Elead);
  weight *= 1./J;

  // Calculate the weight
  weight *= myCS->sigma_eN(Ebeam_int, vk_int, vLead, (lead_type==pCode)); // eN cross section

  // Radiation factors
  if (doRad)
    weight *= radiationFactor(Ebeam, Ek_int, QSq);

  // Final Coulomb, shouldn't effect cross section
  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);

    if (lead_type == pCode) {
      coulombCorrection(vLead_target, deltaECoul); 
    }
    if (rec_type == pCode) {
      coulombCorrection(vRec_target, deltaECoul);
    }
  }
  
}

void QEGenerator::generate_event_lightcone(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  double Estar;
  generate_event_lightcone(weight,lead_type,rec_type,vk_target,vLead_target,vRec_target,vAm2_target,Estar);
}

void QEGenerator::generate_event_lightcone(double &weight, int &lead_type, int &rec_type, TLorentzVector& vk_target, TLorentzVector &vLead_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target, double &Estar)
{
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 4.;

  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type, Estar);

  // Sample decay function
  double alpha1, alphaRec;
  TVector2 v1_perp, vRec_perp;
  decay_function_lc(weight, lead_type, rec_type, alpha1, v1_perp, alphaRec, vRec_perp);
  
  if (weight <= 0.)
    return;
    
  double alphaCM = alpha1 + alphaRec;
  double alphaAm2 = Anum - alphaCM;
  TVector2 vAm2_perp = -1.*v1_perp - vRec_perp;

  TVector3 vbeam_int = vbeam;

  // Initial Coulomb effect
  if(doCoul) {
    double deltaE = calcCoulombEnergy();
    vbeam_int.SetMag(vbeam_int.Mag() + deltaE);
  }
  
  // Initial radiation
  vbeam_int = (doRad ? radiateElectron(vbeam_int) : vbeam_int);
  double Ebeam_int = vbeam_int.Mag();
  TLorentzVector vbeam_int_target(vbeam_int,Ebeam_int);

  // Pick random electron scattering
  double QSq = QSqmin + (QSqmax - QSqmin)*myRand->Rndm();
  double phik = phikmin + (phikmax - phikmin)*myRand->Rndm();
  weight *= (QSqmax - QSqmin) * (phikmax - phikmin);
  
  // Calculate kinematics
  double p1_minus = mbar*alpha1;
  double pRec_minus = mbar*alphaRec;
  double pAm2_minus = mbar*alphaAm2;
  double pRec_plus = (sq(mN) + vRec_perp.Mod2())/pRec_minus;
  double pAm2_plus = (sq(mAm2) + vAm2_perp.Mod2())/pAm2_minus;
  if (mAm2 == 0)
    pAm2_plus = 0.;
  double p1_plus = mA - pRec_plus - pAm2_plus;
  double virt = p1_plus*p1_minus- v1_perp.Mod2() - sq(mN);
  
  // Sovling for energy transfer
  double E1 = 0.5*(p1_plus + p1_minus);
  double p1_z = 0.5*(p1_plus - p1_minus);
  double A = 0.5*(QSq - virt);
  double a = p1_plus*p1_minus;
  double b = -2*E1*A;
  double c = sq(A) - QSq*sq(p1_z);
  double D = sq(b) - 4*a*c;
  
  if (D<0.)
    {
      weight = 0.;
      return;
    }

  // May be up to two valid solutions
  double omega1 = (-b - sqrt(D))/(2.*a);
  double omega2 = (-b + sqrt(D))/(2.*a);
  double omega;

  bool omega1Valid = (omega1>=0.) && (omega1<=Ebeam_int) && ((omega1*E1-A)*p1_z > 0);
  bool omega2Valid = (omega2>=0.) && (omega2<=Ebeam_int) && ((omega2*E1-A)*p1_z > 0);

  if ((!omega1Valid) && (!omega2Valid))
    {
      weight=0.;
      return;
    }
  if (!omega1Valid)
    omega = omega2;
  else if (!omega2Valid)
    omega = omega1;
  else
    {
      omega = (gRandom->Rndm()>0.5)? omega1 : omega2;
      weight *= 2.;
    }

  double Ek_int = Ebeam_int - omega;
  double cosThetak = 1. - QSq/(2.*Ebeam_int*Ek_int);

  if (fabs(cosThetak) > 1.)
    {
      weight=0.;
      return;
    }
  
  double sinThetak = sqrt(1. - sq(cosThetak));
  double kz = Ek_int*cosThetak;
  double kperp = Ek_int*sinThetak;  
  TVector3 vk_int(kperp*cos(phik),kperp*sin(phik),kz);
  
  TLorentzVector vk_int_target;
  vk_int_target.SetVect(vk_int);
  vk_int_target.SetT(Ek_int);

  // Orienting momenta along the q-vector
  TVector3 vq = vbeam_int - vk_int;
  double rot_phi = vq.Phi();
  double rot_theta = vq.Theta();

  TVector3 v1(v1_perp.X(),v1_perp.Y(),p1_z);
  v1.RotateY(rot_theta);
  v1.RotateZ(rot_phi);
  TLorentzVector v1_target(v1,E1);
  
  double ERec = 0.5*(pRec_plus + pRec_minus);
  double pRec_z = 0.5*(pRec_plus - pRec_minus);
  TVector3 vRec(vRec_perp.X(),vRec_perp.Y(),pRec_z);
  vRec.RotateY(rot_theta);
  vRec.RotateZ(rot_phi);
  vRec_target.SetVect(vRec);
  vRec_target.SetT(ERec);
  
  double EAm2 = 0.5*(pAm2_plus + pAm2_minus);
  double pAm2_z = 0.5*(pAm2_plus - pAm2_minus);
  TVector3 vAm2(vAm2_perp.X(),vAm2_perp.Y(),pAm2_z);
  vAm2.RotateY(rot_theta);
  vAm2.RotateZ(rot_phi);
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);
  
  vLead_target = v1_target + vbeam_int_target - vk_int_target;
  TVector3 vLead = vLead_target.Vect();
  double Elead = vLead_target.T();
    
  // Final Radiation
  TVector3 vk = (doRad ? radiateElectron(vk_int) : vk_int);
  double Ek = vk.Mag();
  vk_target.SetVect(vk);
  vk_target.SetT(Ek);
  
  // Jacobian for delta function
  double J = 2.*Ebeam_int*Ek_int*fabs(1. - (vLead.Dot(vq)*omega)/(vq.Mag2()*Elead));
  weight *= 1./J;

  // Calculate the weight
  weight *= myCS->sigma_eN(Ebeam_int, vk_int, vLead, (lead_type==pCode)); // eN cross section

  // Radiation factors
  if (doRad)
    weight *= radiationFactor(Ebeam, Ek_int, QSq);

  // Final Coulomb, shouldn't effect cross section
  if (doCoul) {
    double deltaECoul = calcCoulombEnergy();
    coulombCorrection(vk_target, -deltaECoul);

    if (lead_type == pCode) {
      coulombCorrection(vLead_target, deltaECoul); 
    }
    if (rec_type == pCode) {
      coulombCorrection(vRec_target, deltaECoul);
    }
  }
  
}
