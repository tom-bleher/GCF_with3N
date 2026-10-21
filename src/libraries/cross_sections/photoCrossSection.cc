#include <iostream>
#include <cmath>
#include "photoCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

using namespace std;

photoCrossSection::photoCrossSection()
{
}

photoCrossSection::~photoCrossSection(){}

double photoCrossSection::sigma_pi0_p(double s, double t)
{
  double a = 14.79;
  double b = 7.37;
  double c = 7.49;
  double d = 5.36;

  double m_m = mpi0;
  double m_B = mN;

  double A = s;
  double B = sq(m_B) - sq(m_m) - s;
  double C = sq(m_m);
  double D = sq(B) - 4*A*C;

  double k = (s - sq(mN))/(2.*mN);
  double tmax = sq(m_m) - k*mN/s*(- B - sqrt(D));
  double tmin = sq(m_m) - k*mN/s*(- B + sqrt(D));
  double x = (tmax - t)/(tmax - tmin);

  return pow(a/s,b)*pow(x,c)*exp(d*sq(log(x)));
}

double photoCrossSection::sigma_pip_n(double s, double cosThetaCM)
{

  double A = 9.490;
  double b = 5.329;
  double c = 4.638;

  return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
}

double photoCrossSection::sigma_pim_p(double s, double cosThetaCM)
{

  double A = 10.240;
  double b = 5.329;
  double c = 4.638;

  return pow(A/s,7)*pow(1-cosThetaCM,-b)*pow(1.+cosThetaCM,-c);
}


double photoCrossSection::sigma_rho0_p_old(double s, double cosThetaCM)
{
  const double b=-3.7;
  const double c=-2.2;
  const double a=5.82005e7;

  return 0.75*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
}

double photoCrossSection::sigma_rho0_p(double s, double t, double cosThetaCM)
{

  const double A = 67621.88263175558;
  const double B = 6.052575507613084;
  const double C = 69880481.1275816;
  const double D = 5.133425892196726;
  const double E = 1.885280028435792;

  return A*exp(B*t)+C*pow(s,-7.)*pow(1.2-cosThetaCM,-D)*pow(1.05+cosThetaCM,-E);
}

double photoCrossSection::sigma_rhom_p(double s, double t, double cosThetaCM)
{

  return (1.1/15.9)*sigma_rho0_p(s,t,cosThetaCM);
}

double photoCrossSection::sigma_omega_p_old(double s, double cosThetaCM)
{
  const double b=-3.7;
  const double c=-2.2;
  const double a=5.82005e7;

  return 0.25*(pow(s,-7)*a*pow(1-cosThetaCM,b)*pow(1+cosThetaCM,c));
}

double photoCrossSection::sigma_omega_p(double s, double t, double cosThetaCM)
{

  return sigma_rho0_p(s,t,cosThetaCM)/3.;
}

double photoCrossSection::sigma_phi_p(double s, double t, double cosThetaCM)
{

  const double A = 977.2505868903087;
  const double B = 3.069753904605845;
  const double C = 1.561280750635443;
  const double D = 0.16493020208550144;
  const double E = 1320143.8506911423;

  return A*exp(B*t)+E*pow(s,-7.)*exp(C*sq(cosThetaCM-D));
}

double photoCrossSection::sigma_phi_n(double s, double t, double cosThetaCM)
{

  return sigma_phi_p(s, t, cosThetaCM);
}

double photoCrossSection::sigma_Jpsi_p(double s, double t)
{
  double sig0 = 11.3; //nb
  double beta = 1.3;

  double x = (sq(mJpsi) + 2*mN*mJpsi)/(s - sq(mN));
  if ((x > 1.) or (x < 0.))
    return 0.;

  double sig = sig0 * pow(1. - x,beta);

  double pi = (sq(mN) - s)/(2.*sqrt(s));
  double pf = 0.5*sqrt(sq(sq(mN) - sq(mJpsi))/s - 2.*(sq(mN) + sq(mJpsi)) + s);
  double tmin = sq(sq(mJpsi))/(4.*s) - sq(pi + pf);
  double tmax = sq(sq(mJpsi))/(4.*s) - sq(pi - pf);

  return sig * dipole_F(t, tmin, tmax);
}

double photoCrossSection::sigma_Jpsi_p(double s, double t, double QSq)
{
  double n = 2.44;

  return pow(sq(mJpsi)/(QSq + sq(mJpsi)),n) * sigma_Jpsi_p(s, t);
}

double photoCrossSection::R_Jpsi_p(double QSq)
{
  double a = 2.164;
  double n = 2.131;
  return pow((a*sq(mJpsi) + QSq)/(a*sq(mJpsi)),n) - 1.;
}

double photoCrossSection::dipole_F(double t, double tmin, double tmax)
{
  // |F(t)|^2 with F(t) = 1/(1 - t/mSq)^2, normalized to 1 over [tmax, tmin].
  // mSq is the dipole mass squared in GeV^2: 0.71 is the proton dipole of
  // Lomnitz & Klein, arXiv:1803.06420 eq. (15). 1.14 is the code's original
  // value (d744e3e, no source given); it is close to the Frankfurt-Strikman
  // two-gluon scale m^2 ~ 1.1 GeV^2, which may be where it came from.
  double mSq = 0.71;
  //double mSq = 1.14;

  return pow(mSq - t,-4.) * 3./(pow(mSq - tmin,-3.) - pow(mSq - tmax,-3.));
}

double photoCrossSection::sigma_deltapp_pim(double s, double cosThetaCM)
{
  double a =50955200;
  double b =2.93657;
  double c =1.74578;

  return pow(s,-7)*a*pow(1-cosThetaCM,-b)*pow(1+cosThetaCM,-c);
}

double photoCrossSection::sigma_deltap_pim(double s, double cosThetaCM)
{

  return sigma_deltapp_pim(s,cosThetaCM);
}
