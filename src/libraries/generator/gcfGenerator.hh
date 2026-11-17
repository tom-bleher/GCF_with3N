#ifndef __GCF_GENERATOR_H__
#define __GCF_GENERATOR_H__

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "nucleus/gcfNucleus.hh"
#include <algorithm>
#include <vector>

class gcfGenerator
{
 public:
  gcfGenerator(gcfNucleus * thisInfo, TRandom3 * thisRand);
  ~gcfGenerator();
  void set_doRad(bool rad);
  void set_doCoul(bool coul);
  void set_phiRel_range(double low, double high);
  void set_phiRel_range_deg(double low, double high);
  void set_thetaRel_range(double low, double high);
  void set_thetaRel_range_deg(double low, double high);
  void set_pRel_range(double low, double high);
  void set_pRel_cut(double new_cutoff);
  // Lower end of the sampled |p_rel| range: the larger of pRelmin and the cut.
  double pRel_floor() const { return std::max(pRelmin, pRel_cut); }

  // Importance-sample |p_rel| from its own integrand k^2 * S(k) instead of
  // uniformly over [pRel_floor(), pRelmax].
  //
  // With uniform sampling the weight carries the full k^2 S(k) factor, and
  // S(k) spans several decades, so a handful of low-k events dominate every
  // sum. Measured on a 20k PWIA test: N_eff / N = 2.5% for a fixed k>0.25
  // floor and 0.11% once the floor follows the local Fermi momentum (which
  // reaches down to ~0.02 GeV/c at large radius) -- i.e. 22 effective events
  // out of 20000. Drawing k from k^2 S(k) instead makes that dimension's
  // weight contribution the integral over [pRel_floor(), pRelmax], which is
  // constant per event, so the spread vanishes.
  //
  // Expectation values are unchanged (a residual f/f_hat factor keeps the
  // estimator exactly unbiased against the piecewise-linear cumulative), only
  // the variance. Off by default so existing users are untouched.
  void set_pRel_importance_sampling(bool on = true) { fPRelIS = on; }
  bool pRel_importance_sampling() const { return fPRelIS; }
  void set_nu_range(double low, double high);
  void set_xB_range(double low, double high);
  void set_QSq_range(double low, double high);
  void set_phik_range(double low, double high);
  void set_phik_range_deg(double low, double high);
  void set_deuteron();

  double get_QSqmin() const { return QSqmin; }
  double get_QSqmax() const { return QSqmax; }
  double get_xBmin() const { return xBmin; }
  double get_xBmax() const { return xBmax; }

  void randomize_cutoff();

  bool parse_phase_space_file(char* phase_space);
  
  void decay_function(double &weight, int lead_type, int rec_type, TVector3 &vi, TVector3 &vRec);
  void decay_function_lc(double &weight, int lead_type, int rec_type, double &alphai, TVector2 &vi_perp, double &alphaRec, TVector2 &vRec_perp);
  
  protected:
  double get_mAm2(int lead_type, int rec_type);
  double get_mAm2(int lead_type, int rec_type, double &Estar);
  // A-1 residual mass for single-nucleon (mean-field) knockout of `lead_type`.
  double get_mAm1(int lead_type);

  void t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4);
  void t_scatter(double &weight, double m3, double m4, TLorentzVector v1, TLorentzVector v2, TLorentzVector &v3, TLorentzVector &v4, double &cosThetaCM);

  TVector3 radiateElectron(TVector3 ve);
  double deltaHard(double QSq);
  double radiationFactor(double Ebeam, double Ek, double QSq);

  double calcCoulombEnergy();
  void coulombCorrection(TLorentzVector &p, double deltaE);
  
  gcfNucleus * myInfo;
  TRandom3 * myRand;  

  int Anum;
  double mA;
  double mbar;
  double mAmpp;
  double mAmpn;
  double mAmnn;
  double sigCM;
  bool random_Estar;
  bool doRad;
  bool doCoul;

  double pRel_cut = 0.25;
  double pRel_cut_range = 0.05;
  
  double pRelmin;
  double pRelmax;
  double phiRelmin;
  double phiRelmax;
  double thetaRelmin;
  double thetaRelmax; 
  double cosThetaRelmin;
  double cosThetaRelmax;
  double alphaRelmin;
  double alphaRelmax;
  
  double numin;
  double numax;
  double xBmin;
  double xBmax;
  double QSqmin;
  double QSqmax;
  double phikmin;
  double phikmax;

  // ---- |p_rel| importance sampling ----
  bool fPRelIS = false;
  bool fPRelIS_built = false;
  static const int kISBins = 8000;      // grid over [0, kISmax]
  double fISkmax = 2.0;                 // GeV/c, beyond the AV18 table range
  double fISdk = 0.;
  // Cumulative F(k) = int_0^k k'^2 S(k') dk' for the three S() outcomes:
  // 0 = pp, 1 = nn, 2 = pn.
  std::vector<double> fISF[3];
  void build_pRel_IS_tables();
  int is_class(int lead_type, int rec_type) const;
  // Cumulative at k (linear interpolation on the grid).
  double IS_F(int cls, double k) const;
  // Inverse cumulative: returns k with F(k) = u, and sets f_hat to the
  // piecewise-constant density implied by the interpolated cumulative.
  double IS_invF(int cls, double u, double &f_hat) const;
  // int_lo^hi k^2 S(k) dk for the pair (lead_type, rec_type). Builds the
  // tables on first use. Used by derived generators that need to know how
  // much SRC strength a given |p_rel| window holds.
  double pRel_integral(int lead_type, int rec_type, double lo, double hi);

};

#endif
