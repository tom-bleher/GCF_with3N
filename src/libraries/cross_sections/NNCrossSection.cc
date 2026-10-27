#include <iostream>
#include <cmath>
#include <cstdlib>
#include "NNCrossSection.hh"
#include "helpers.hh"
#include "constants.hh"

// pp elastic dsigma/dt, in mb/GeV^2, from two sources:
//
// 1. T_lab = 50-1000 MeV (s = 3.62-5.40 GeV^2): the pp elastic dsigma/dOmega_cm
//    tables (mb/sr) that GENIE R-3.06.02 INTRANUKE uses for hN scattering,
//    Generator-R-3_06_02/data/evgen/intranuke/diff_ang/pp/pp<T>.txt, loaded by
//    INukeHadroData2018.cxx. Per INukeHadroData2018.h they come "mostly from
//    the SAID (Arndt et al.) PWA fit" (R.A. Arndt et al., PRC 76, 025209 (2007),
//    doi:10.1103/PhysRevC.76.025209). genIon does not link GENIE, so the 0-90 deg
//    half of each table is copied below (the tables are symmetric about 90 deg
//    to 1e-3). As in GENIE, dsigma/dOmega is interpolated bilinearly in
//    (T, theta_cm). dsigma/dt = pi/p_cm^2 dsigma/dOmega. Integrated over
//    0-90 deg they give PDG sigma_el (sigma_tot below pion threshold) to 5% for
//    T = 150-1000 MeV. At T = 50 and 100 MeV they are 12% and 6-9% above it:
//    the 0 and 9 deg entries include Coulomb scattering (5.5x and 2.3x the
//    18 deg entry), and between 50 and 100 MeV the linear interpolation in T
//    overshoots (+21% at T = 63 MeV).
//
// 2. p_lab = 1.8-6.0 GeV/c (s = 5.5-13.2 GeV^2): a fit to 1667 published pp
//    dsigma/dt points from HEPData,
//      Ambats et al., PRD 9, 1179 (1974), doi:10.17182/hepdata.3409.v1
//      Williams et al., Nuovo Cim. A 8, 447 (1972), doi:10.17182/hepdata.37468.v1
//      Brabson et al., INSPIRE 1085234 (unpublished), doi:10.17182/hepdata.47019.v1
//      Coletti et al., Nuovo Cim. A 49, 479 (1967), doi:10.17182/hepdata.981.v1
//      Jenkins et al., PRL 40, 425 (1978), doi:10.17182/hepdata.3359.v1
//    with |t| >= 0.04 GeV^2 (above the Coulomb-nuclear interference region),
//    errors stat (+) 7% normalization, fitted in log space: chi2/ndf =
//    2567/1656. The form is symmetric in t <-> u:
//      dsigma/dt = F(s,t) + F(s,u),  u = 4m^2 - s - t,
//      F = A exp(B t) + C exp(D t) + G exp(H t),
//      A = 96.435 (s/9)^-0.47018    B = 7.3600 + 0.40354 (s-9)
//      C = 3.3269 (s/9)^2.5525      D = 2.1224 + 0.40754 (s-9)
//      G = 0.18926 (s/9)^-4.2211    H = 0.63573
//    (s in GeV^2, t in GeV^2). Model/data medians per data set are 0.81-1.20;
//    Brabson sits 10-15% below the others. Integrated over 0-90 deg it gives
//    sigma_el = 17.2, 14.3, 12.5, 11.3 mb at 3, 4, 5, 6 GeV/c (PDG: 17.2,
//    13.5, 12.7, 11.5). At 90 deg it falls as s^-9.1 near s = 9.4 GeV^2
//    (fixed-angle counting rule: s^-10).
//
// Between T = 850 and 1000 MeV the two are blended linearly in s; there they
// agree on sigma_el to 2% and on dsigma/dt to 15% at any angle.
//
// Outside s = 3.62-13.17 GeV^2 (T_lab = 50 MeV to p_lab = 6 GeV/c) there are
// no inputs here; see ppRangePolicy in NNCrossSection.hh. Above 13.17 GeV^2
// the extrapolated fit stays within a factor 1.1-2.0 (median per momentum) of
// the large-angle data of Jenkins et al. up to 8.8 GeV/c (s = 18.4 GeV^2).
// Below 3.62 GeV^2 holding dsigma/dt (rather than dsigma/dOmega) fixed keeps
// the weight finite as p_cm -> 0.
//
// This replaces the parameterization ported from V. Panin's R3BRoot QFS
// generator (qfs.cc, 2012). Its slope B(s) matched the data, but the
// normalization A(s) added in the GCF port was off by factors 0.16 to 1.8
// against PDG sigma_el, the polynomials diverge above s = 8.3 GeV^2, and it
// returned 0 outside 3.8 < s < 8.3 GeV^2 (issue #138).

namespace
{
  // GENIE pp tables: T = 50, 100, ..., 1000 MeV; theta_cm = 0, 9, ..., 90 deg.
  const int nTab_T = 20;
  const int nTab_th = 11;
  const double tab_T0 = 0.050; // GeV
  const double tab_dT = 0.050; // GeV
  const double tab_dth = 9.;   // deg
  const double dsdO_tab[nTab_T][nTab_th] = {
    {40.77, 40.77, 7.355, 7.981, 8.464, 8.628, 8.653, 8.624, 8.583, 8.55, 8.538}, // T = 50 MeV
    {9.422, 9.422, 4.099, 4.655, 4.814, 4.822, 4.774, 4.707, 4.645, 4.601, 4.585}, // T = 100 MeV
    {5.941, 5.261, 3.64, 3.925, 3.99, 3.99, 3.965, 3.931, 3.894, 3.865, 3.854}, // T = 150 MeV
    {5.539, 4.338, 3.556, 3.702, 3.751, 3.778, 3.795, 3.795, 3.779, 3.761, 3.753}, // T = 200 MeV
    {5.377, 4.15, 3.587, 3.667, 3.717, 3.763, 3.795, 3.8, 3.781, 3.757, 3.747}, // T = 250 MeV
    {5.392, 4.22, 3.714, 3.751, 3.79, 3.828, 3.844, 3.824, 3.78, 3.74, 3.723}, // T = 300 MeV
    {5.675, 4.463, 3.969, 3.931, 3.92, 3.92, 3.899, 3.845, 3.774, 3.718, 3.698}, // T = 350 MeV
    {6.089, 4.816, 4.295, 4.163, 4.08, 4.025, 3.953, 3.854, 3.754, 3.683, 3.658}, // T = 400 MeV
    {6.659, 5.292, 4.724, 4.467, 4.282, 4.145, 4.003, 3.85, 3.716, 3.628, 3.598}, // T = 450 MeV
    {7.382, 5.94, 5.308, 4.891, 4.561, 4.296, 4.045, 3.808, 3.622, 3.507, 3.469}, // T = 500 MeV
    {8.196, 6.802, 6.084, 5.484, 4.956, 4.496, 4.063, 3.682, 3.4, 3.231, 3.175}, // T = 550 MeV
    {9.155, 7.947, 7.053, 6.218, 5.431, 4.701, 4.013, 3.432, 3.012, 2.764, 2.681}, // T = 600 MeV
    {10.51, 9.479, 8.196, 7.002, 5.872, 4.823, 3.864, 3.088, 2.546, 2.227, 2.12}, // T = 650 MeV
    {12.17, 11.25, 9.422, 7.774, 6.246, 4.863, 3.662, 2.74, 2.121, 1.763, 1.643}, // T = 700 MeV
    {13.8, 12.97, 10.56, 8.46, 6.54, 4.849, 3.451, 2.436, 1.783, 1.414, 1.29}, // T = 750 MeV
    {15.24, 14.45, 11.52, 9.003, 6.733, 4.782, 3.241, 2.178, 1.524, 1.164, 1.045}, // T = 800 MeV
    {16.52, 15.72, 12.3, 9.399, 6.824, 4.664, 3.03, 1.958, 1.326, 0.9874, 0.8769}, // T = 850 MeV
    {17.71, 16.82, 12.94, 9.677, 6.832, 4.509, 2.821, 1.766, 1.17, 0.8596, 0.7599}, // T = 900 MeV
    {18.85, 17.84, 13.49, 9.865, 6.778, 4.326, 2.616, 1.595, 1.042, 0.7637, 0.6755}, // T = 950 MeV
    {19.99, 18.8, 13.97, 9.987, 6.679, 4.128, 2.417, 1.442, 0.9359, 0.6888, 0.6119}, // T = 1000 MeV
  };

  // Blend region between the tables and the fit, in T_lab (GeV)
  const double blend_T_lo = 0.850;
  const double blend_T_hi = 1.000;
  // Upper end of the fitted data, in p_lab (GeV/c)
  const double fit_p_max = 6.0;

  double s_of_T(double T)
  {
    return 4.*sq(mN) + 2.*mN*T;
  }

  // dsigma/dt (mb/GeV^2) from the GENIE/SAID tables at s and |cos theta_cm| = x
  double dsdt_table(double s, double x)
  {
    double T = (s - 4.*sq(mN))/(2.*mN);
    double pcm2 = (s - 4.*sq(mN))/4.;

    double fT = (T - tab_T0)/tab_dT;
    if (fT < 0.) fT = 0.;
    if (fT > nTab_T - 1) fT = nTab_T - 1;
    int i = (int)fT;
    if (i > nTab_T - 2) i = nTab_T - 2;
    fT -= i;

    double fth = acos(x)*180./M_PI/tab_dth;
    if (fth > nTab_th - 1) fth = nTab_th - 1;
    int j = (int)fth;
    if (j > nTab_th - 2) j = nTab_th - 2;
    fth -= j;

    double dsdO = (1.-fT)*((1.-fth)*dsdO_tab[i][j]   + fth*dsdO_tab[i][j+1])
                +     fT *((1.-fth)*dsdO_tab[i+1][j] + fth*dsdO_tab[i+1][j+1]);

    return dsdO*M_PI/pcm2;
  }

  double F_fit(double s, double t)
  {
    double A = 96.435*pow(s/9., -0.47018);
    double B = 7.3600 + 0.40354*(s - 9.);
    double C = 3.3269*pow(s/9., 2.5525);
    double D = 2.1224 + 0.40754*(s - 9.);
    double G = 0.18926*pow(s/9., -4.2211);
    double H = 0.63573;
    return A*exp(B*t) + C*exp(D*t) + G*exp(H*t);
  }

  // dsigma/dt (mb/GeV^2) from the HEPData fit at s and |cos theta_cm| = x
  double dsdt_fit(double s, double x)
  {
    double t = -0.5*(s - 4.*sq(mN))*(1. - x);
    double u = -0.5*(s - 4.*sq(mN))*(1. + x);
    return F_fit(s,t) + F_fit(s,u);
  }
}

NNCrossSection::NNCrossSection()
{
  // Set defaults
  myParam=ppData;
  myPolicy=ppExtend;
  warned_low=false;
  warned_high=false;
}

NNCrossSection::NNCrossSection(csParam thisParam)
{
  std::cerr << "NNCrossSection: you have selected parameterization: " << thisParam <<"\n";
  myParam=thisParam;
  myPolicy=ppExtend;
  warned_low=false;
  warned_high=false;
}

NNCrossSection::~NNCrossSection(){}

void NNCrossSection::set_pp_range_policy(ppRangePolicy thisPolicy)
{
  myPolicy=thisPolicy;
}

double NNCrossSection::pp_s_min()
{
  return s_of_T(tab_T0);
}

double NNCrossSection::pp_s_max()
{
  return 2.*sq(mN) + 2.*mN*sqrt(sq(fit_p_max) + sq(mN));
}

double NNCrossSection::sigma_pp(double s, double t)
{
  switch (myParam)
    {
    case ppData:
      return sigma_pp_data(s,t);
    default:
      {
	std::cerr << "Invalid cross section parameterization! Double check and fix!\n";
        exit(-1);
      }
    }
  return 0;
}

double NNCrossSection::sigma_pp_data(double s, double t)
{
  double q = s - 4.*sq(mN); // 4 p_cm^2
  if (q <= 0.)
    return 0.;

  // c.m. angle from t (on-shell relation), folded about 90 deg because the
  // protons are identical. For an off-shell struck nucleon t can lie outside
  // the on-shell range [4m^2 - s, 0]; it is clamped to that range.
  double x = 1. + 2.*t/q;
  if (x > 1.) x = 1.;
  if (x < -1.) x = -1.;
  x = fabs(x);

  double s_eval = s;
  if (s < pp_s_min())
    {
      if (not warned_low)
	{
	  std::cerr << "NNCrossSection: pp s = " << s << " GeV^2 is below the range of the pp data ("
		    << pp_s_min() << " GeV^2, T_lab = 50 MeV); "
		    << ((myPolicy==ppZero) ? "setting dsigma/dt to 0" : "holding dsigma/dt at its value there at the same c.m. angle")
		    << " for such s. Warning printed once.\n";
	  warned_low=true;
	}
      if (myPolicy==ppZero)
	return 0.;
      s_eval = pp_s_min();
    }
  else if (s > pp_s_max())
    {
      if (not warned_high)
	{
	  std::cerr << "NNCrossSection: pp s = " << s << " GeV^2 is above the range of the pp data ("
		    << pp_s_max() << " GeV^2, p_lab = " << fit_p_max << " GeV/c); "
		    << ((myPolicy==ppZero) ? "setting dsigma/dt to 0" : "extrapolating the fit")
		    << " for such s. Warning printed once.\n";
	  warned_high=true;
	}
      if (myPolicy==ppZero)
	return 0.;
    }

  double s_lo = s_of_T(blend_T_lo);
  double s_hi = s_of_T(blend_T_hi);
  if (s_eval <= s_lo)
    return dsdt_table(s_eval,x);
  if (s_eval >= s_hi)
    return dsdt_fit(s_eval,x);
  double w = (s_eval - s_lo)/(s_hi - s_lo);
  return (1.-w)*dsdt_table(s_eval,x) + w*dsdt_fit(s_eval,x);
}
