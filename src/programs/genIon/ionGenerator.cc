#include "ionGenerator.hh"
#include "constants.hh"
#include "helpers.hh"

using namespace std;

ionGenerator::ionGenerator(double E, gcfNucleus * thisInfo, NNCrossSection * thisCS, TRandom3 * thisRand) : gcfGenerator(thisInfo, thisRand)
{

  myCS = thisCS;
  
  Ebeam = E;
  E1 = Ebeam + mN;
  p1 = sqrt(sq(E1) - sq(mN));
  v1.SetXYZ(0.,0.,p1);
  v1_target.SetXYZT(0.,0.,p1,E1);

  nCalls_pp = nLow_pp = nHigh_pp = 0;
  sumW_pp = sumWLow_pp = sumWHigh_pp = 0.;

}

ionGenerator::~ionGenerator()
{
}

void ionGenerator::generate_event(double &weight, int &lead_type, int &rec_type, TLorentzVector &v3_target, TLorentzVector &v4_target, TLorentzVector &vRec_target, TLorentzVector &vAm2_target)
{
  
  // Start with weight 1. Only multiply terms to weight. If trouble, set weight=0.
  weight = 1.;

  // Decide what kind of proton or neutron pair we are dealing with
  lead_type = pCode;
  rec_type = (myRand->Rndm() > 0.5) ? pCode:nCode;
  weight *= 2.;
  
  // Determine mass of A-2 system
  double mAm2 = get_mAm2(lead_type, rec_type);

  TVector3 v2, vRec;

  decay_function(weight, lead_type, rec_type, v2, vRec);
  
  if (weight <= 0.)
    return;
      
  TVector3 vAm2 = - v2 - vRec;
  double EAm2 = sqrt(vAm2.Mag2() + sq(mAm2));
  vAm2_target.SetVect(vAm2);
  vAm2_target.SetT(EAm2);
  
  double Erec = sqrt(sq(mN) + vRec.Mag2());
  vRec_target.SetVect(vRec);
  vRec_target.SetT(Erec);
  
  double E2 = mA - EAm2 - Erec;
  
  TLorentzVector v2_target(v2,E2);
  
  t_scatter(weight, mN, mN, v1_target, v2_target, v3_target, v4_target);
  
  if (weight <= 0.)
    return;
  
  double s = 2.*sq(mN) + 2.*v3_target.Dot(v4_target);
  double t = 2.*sq(mN) - 2.*v1_target.Dot(v3_target);
  double u = 2.*sq(mN) - 2.*v1_target.Dot(v4_target);
  
  // Calculate the flux factor on the cross section
  double v12 = sqrt(sq(v1_target.Dot(v2_target)) - sq(mN)*v2_target.Mag2())/(E1*E2);
  double v1A = p1/E1;
  
  // Calculate the weight. The two outgoing protons are identical and
  // t_scatter samples the full sphere, so every final state is generated
  // twice (3 <-> 4) and the full-sphere integral of the t <-> u symmetric
  // dsigma/dt is 2 sigma_el. The 1/2 makes the sum of weights count each
  // pp final state once.
  weight *= 0.5 * v12/v1A*myCS->sigma_pp(s,max(t,u)); // Scattering cross section

  // Bookkeeping of the weight that comes from s outside the pp data
  nCalls_pp++;
  sumW_pp += weight;
  if (s < NNCrossSection::pp_s_min())
    {
      nLow_pp++;
      sumWLow_pp += weight;
    }
  else if (s > NNCrossSection::pp_s_max())
    {
      nHigh_pp++;
      sumWHigh_pp += weight;
    }
  
}

double ionGenerator::get_s_nominal()
{
  // pp s for the beam nucleon on a free nucleon at rest
  return 2.*sq(mN) + 2.*mN*E1;
}

void ionGenerator::print_pp_summary()
{
  std::cerr << "genIon: pp cross section evaluated for " << nCalls_pp << " events; "
	    << "pp data cover s = " << NNCrossSection::pp_s_min() << "-" << NNCrossSection::pp_s_max() << " GeV^2.\n"
	    << "  s below: " << nLow_pp << " events, " << ((sumW_pp > 0.) ? 100.*sumWLow_pp/sumW_pp : 0.) << "% of the total weight\n"
	    << "  s above: " << nHigh_pp << " events, " << ((sumW_pp > 0.) ? 100.*sumWHigh_pp/sumW_pp : 0.) << "% of the total weight\n"
	    << "  outside the range dsigma/dt was "
	    << ((myCS->get_pp_range_policy()==ppZero) ? "set to 0" : "held at its lower-edge value (below) or extrapolated (above)") << ".\n";
}
