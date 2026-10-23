#ifndef __NN_CROSS_SECTIONS_H__
#define __NN_CROSS_SECTIONS_H__

// pp elastic parameterization. ppData: SAID-based tables (GENIE INTRANUKE)
// for T_lab <= 1 GeV joined to a fit to HEPData dsigma/dt for 1.8-6 GeV/c.
// See NNCrossSection.cc for sources and validity.
enum csParam {ppData};

// What sigma_pp does for s outside [pp_s_min(), pp_s_max()], the range
// covered by the pp data. ppExtend: below, dsigma/dt is held at its
// T_lab = 50 MeV value at the same c.m. angle; above, the fit is extrapolated.
// ppZero: dsigma/dt = 0. Either way a warning is printed once per side.
enum ppRangePolicy {ppExtend, ppZero};

class NNCrossSection
{
 public:
  NNCrossSection();
  NNCrossSection(csParam thisParam);
  ~NNCrossSection();
  // pp elastic dsigma/dt (mb/GeV^2) at Mandelstam s and t (GeV^2).
  // Symmetric under t <-> u (identical protons): integrating it over the
  // full sphere gives 2 sigma_el.
  double sigma_pp(double s, double t);
  double sigma_pp_data(double s, double t);
  void set_pp_range_policy(ppRangePolicy thisPolicy);
  ppRangePolicy get_pp_range_policy() const {return myPolicy;};
  static double pp_s_min();
  static double pp_s_max();

 private:
  csParam myParam;
  ppRangePolicy myPolicy;
  bool warned_low;
  bool warned_high;

};

#endif
