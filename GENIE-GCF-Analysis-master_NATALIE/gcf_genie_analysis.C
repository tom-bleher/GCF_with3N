//#define genie_loop_cxx
#include "TH1.h"
#include "TH2.h"

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "GeneratorComparison/AccMap.cpp"


#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <cmath>
#include "TVector3.h"

#include <iostream>
#include <fstream>
double sq(double x)
{
  return x*x;
}

// ############## Custom proton gap parameters ##################
double gap_function(double p, const double * params) { return (params[0] - params[1]*exp(-p/params[2])); };

const int nGaps_per_sector[6]={2,2,4,1,3,2};


const double gap_params_0[2][2][3] = {{{35.,10.,0.4},{48.,18.,0.6}},   // order is: gap, lower-upper, a-b-c
				      {{106.,10.,0.4},{111.,6.,0.4}}};
const double gap_params_1[2][2][3] = {{{88.,15.,0.6},{96.,15.,0.6}},
				      {{110.,15.,0.6},{180.,0.,1.}}};
const double gap_params_2[4][2][3] = {{{19.,50.,0.7},{24.,50.,0.7}},
				      {{35.,10.,0.6},{44.,13.,0.6}},
				      {{74.,25.,0.4},{108.,20.,0.4}},
				      {{115.,15.,0.6},{180.,0.,1.}}};
const double gap_params_3[1][2][3] = {{{100.,18.,0.6},{112,18,0.6}}};
const double gap_params_4[3][2][3] = {{{34.,60.,0.4},{38.,50.,0.4}},
				      {{40.,50.,0.4},{44.,40.,0.4}},
				      {{92.,15.,0.6},{106.,8.,0.6}}};
const double gap_params_5[2][2][3] = {{{82.,20.,0.4},{91.,20.,0.6}},
				      {{114.,2.,0.5},{180.,0.,1.}}};

const void * gap_params[6] = {gap_params_0,gap_params_1,gap_params_2,gap_params_3,gap_params_4,gap_params_5};

// ############## PROTON FIDUCIAL CUT PARAMETERS ################
const Double_t kThetaPar0PiPlus[6] = {  7.00823   ,  5.5        ,  7.06596   ,  6.32763   ,  5.5       ,  5.5      };
const Double_t kThetaPar1PiPlus[6] = {  0.207249  ,  0.1        ,  0.127764  ,  0.1       ,  0.211012  ,  0.281549 };
const Double_t kThetaPar2PiPlus[6] = {  0.169287  ,  0.506354   , -0.0663754 ,  0.221727  ,  0.640963  ,  0.358452 };
const Double_t kThetaPar3PiPlus[6] = {  0.1       ,  0.1        ,  0.100003  ,  0.1       ,  0.1       ,  0.1      };
const Double_t kThetaPar4PiPlus[6] = {  0.1       ,  3.30779    ,  4.499     ,  5.30981   ,  3.20347   ,  0.776161 };
const Double_t kThetaPar5PiPlus[6] = { -0.1       , -0.651811   , -3.1793    , -3.3461    , -1.10808   , -0.462045 };

const Double_t kFidPar0Low0PiPlus[6] = {  25.      ,  25.      ,  25.     ,  25.       ,  25.       ,  25.       };
const Double_t kFidPar0Low1PiPlus[6] = { -12.      , -12.      , -12.     , -12        , -12        , -12        };
const Double_t kFidPar0Low2PiPlus[6] = {   1.64476 ,   1.51915 ,   1.1095 ,   0.977829 ,   0.955366 ,   0.969146 };
const Double_t kFidPar0Low3PiPlus[6] = {   4.4     ,   4.4     ,   4.4    ,   4.4      ,   4.4      ,   4.4      };

const Double_t kFidPar1Low0PiPlus[6] = {  4.        ,  4.   ,  2.78427 ,  3.58539 ,  3.32277   ,  4.      };
const Double_t kFidPar1Low1PiPlus[6] = {  2.        ,  2.   ,  2.      ,  1.38233 ,  0.0410601 ,  2.      };
const Double_t kFidPar1Low2PiPlus[6] = { -0.978469  , -2.   , -1.73543 , -2.      , -0.953828  , -2.      };
const Double_t kFidPar1Low3PiPlus[6] = {  0.5       ,  0.5  ,  0.5     ,  0.5     ,  0.5       ,  1.08576 };

const Double_t kFidPar0High0PiPlus[6] = {  25.      , 24.8096  , 24.8758  ,  25.       , 25.       , 25.      };
const Double_t kFidPar0High1PiPlus[6] = { -11.9735  , -8.      , -8.      , -12.       , -8.52574  , -8.      };
const Double_t kFidPar0High2PiPlus[6] = {  0.803484 ,  0.85143 ,  1.01249 ,   0.910994 ,  0.682825 ,  0.88846 };
const Double_t kFidPar0High3PiPlus[6] = {  4.40024  ,  4.8     ,  4.8     ,   4.4      ,  4.79866  ,  4.8     };

const Double_t kFidPar1High0PiPlus[6] = {  2.53606  ,  2.65468  ,  3.17084 ,  2.47156 ,  2.42349  ,  2.64394 };
const Double_t kFidPar1High1PiPlus[6] = {  0.442034 ,  0.201149 ,  1.27519 ,  1.76076 ,  1.25399  ,  0.15892 };
const Double_t kFidPar1High2PiPlus[6] = { -2.       , -0.179631 , -2.      , -1.89436 , -2.       , -2.      };
const Double_t kFidPar1High3PiPlus[6] = {  1.02806  ,  1.6      ,  0.5     ,  1.03961 ,  0.815707 ,  1.31013 };

//how to read it in. Goes parameter by parameter and puts it in the proper array

double a(double mom, double p0, double p1, double p2, double p3)
{
  return p0 + p1*exp(p2*(mom-p3));
}

double b(double mom, double p0, double p1, double p2, double p3)
{
  return p0 + p1*mom*exp(p2*sq(mom-p3));
}

double theta_min(double mom, double p0, double p1, double p2, double p3, double p4, double p5)
{
  return p0 + p1/sq(mom) + p2*mom + p3/mom + p4*exp(mom*p5);
}

double deltaPhi(double theta, double a, double b, double thetaMin)
{
  return a*(1. - 1./((theta-thetaMin)/b + 1.));
}

bool gap_accept(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;

  int sector = (phi+30.)/60.;

  // Do a quick test for sector gaps
  for (int gap=0 ; gap < nGaps_per_sector[sector] ; gap++)
    {
      //std::cout << "Testing sector " << sector << " gap " << gap << "\n";

      const double * min_params = &(((const double*)gap_params[sector])[gap*6 + 0]); // Multi-dim arrays are 'made flat'
      const double * max_params = &(((const double*)gap_params[sector])[gap*6 + 3]);

      double min_theta = gap_function(mom,min_params);
      double max_theta = gap_function(mom,max_params);

      //double min_theta = gap_function(mom,(double*)((double**)((double ***)gap_params[sector])[gap])[0]);
      //double max_theta = gap_function(mom,(double*)((double**)((double ***)gap_params[sector])[gap])[1]);

      if ((theta > min_theta) && (theta < max_theta))
	return false;
    }

  return true;

}

bool accept_proton_simple(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;


  int sector = (phi+30.)/60.;


  // Now that we have calculated the gaps, figure out if we are in the conventional
  // fiducial wedges.

  double minTheta = theta_min(mom,kThetaPar0PiPlus[sector],
			      kThetaPar1PiPlus[sector],
			      kThetaPar2PiPlus[sector],
			      kThetaPar3PiPlus[sector],
			      kThetaPar4PiPlus[sector],
			      kThetaPar5PiPlus[sector]);

  // We must be on the correct side of the wedge's tip.
  if (theta < minTheta)
    return false;

  double phiCentral = 60.*sector;

  double aLow = a(mom,kFidPar0Low0PiPlus[sector],
		  kFidPar0Low1PiPlus[sector],
		  kFidPar0Low2PiPlus[sector],
		  kFidPar0Low3PiPlus[sector]);

  double aHigh = a(mom,kFidPar0High0PiPlus[sector],
		  kFidPar0High1PiPlus[sector],
		  kFidPar0High2PiPlus[sector],
		  kFidPar0High3PiPlus[sector]);

  double bLow = b(mom,kFidPar1Low0PiPlus[sector],
		  kFidPar1Low1PiPlus[sector],
		  kFidPar1Low2PiPlus[sector],
		  kFidPar1Low3PiPlus[sector]);

  double bHigh = b(mom,kFidPar1High0PiPlus[sector],
		  kFidPar1High1PiPlus[sector],
		  kFidPar1High2PiPlus[sector],
		  kFidPar1High3PiPlus[sector]);

  double deltaPhiLow = deltaPhi(theta, aLow, bLow, minTheta);
  double deltaPhiHigh = deltaPhi(theta, aHigh, bHigh, minTheta);

  //std::cout << mom << " " << theta << " " << phi << " " << sector << " " << aLow << " " << aHigh << "\n";

  if ((phi < phiCentral + deltaPhiHigh) && (phi > phiCentral - deltaPhiLow))
    return true;

  return false;
}

bool accept_proton(TVector3 p)
{
  if (gap_accept(p))
    return accept_proton_simple(p);
  return false;
}

// ############## ELECTRON FIDUCIAL CUT PARAMETERS (starting with e) #################
// sector, side
const double eFidAp0[6][2]={{25.9392,25.142},{26.3123,26.546},{25.2953,26.4484},{28.,23.},{28.,23.2826},{28.,23.8565}};
const double eFidAp1[6][2]={{-12.,-12.},{-12.,-12.},{-12.,-12.},{-12.,-12.},{-10.802,-12.},{-11.9703,-12.}};
const double eFidAp2[6][2]={{0.453685,0.357039},{0.40864,0.425442},{0.486011,0.519846},{0.366217,0.562825},{0.203741,1.15862},{0.283665,1.4323}};
const double eFidAp3[6][2]={{4.4,4.4},{4.4,4.4},{4.4,4.4},{4.4,4.4},{4.40456,4.4},{4.4497,4.4}};

const double eFidBp0[6][2]={{2.53279,2.25185},{2.9569,3.27384},{1.98429,1.5},{2.29784,1.65084},{2.79198,1.92076},{1.84305,3.18777}};
const double eFidBp1[6][2]={{1.19464,0.527297},{2.,2.},{1.914,1.81489},{0.351108,1.91044},{0.534785,0.90012},{0.883945,0.539951}};
const double eFidBp2[6][2]={{-2.,-0.1},{-2.,-2.},{-0.5088,-0.273122},{-0.1,-0.574873},{-2.,-0.112026},{-0.1,-2.}};
const double eFidBp3[6][2]={{1.25872,1.55534},{0.74865,0.519786},{0.721026,0.5},{1.6,0.5},{0.720147,0.692653},{1.6,1.02089}};

const double eFidThetaMinp0[6]={15.,13.,13.,13.,14.6669,14.0658};
const double eFidThetaMinp1[6]={4.,1.99082,3.82167,1.50112,1.79097,2.20492};
const double eFidThetaMinp2[6]={-0.828672,-0.159112,-0.0647212,-1.3,-0.355746,-1.3};
const double eFidThetaMinp3[6]={2.51825,0.100001,0.100008,0.1,0.1,0.1};
const double eFidThetaMinp4[6]={2.15567,13.6076,9.26858,12.3802,14.1265,8.90341};
const double eFidThetaMinp5[6]={-0.1,-0.554978,-0.448555,-0.200086,-0.777645,-0.1};
//##############################################################3

// For FidThetaMin calculation for electron
const Double_t kThetaPar0[6] = { 15        , 15        ,  15       , 15        ,  13       ,  13        };
const Double_t kThetaPar1[6] = { -0.425145 , -1.02217  , -0.7837   , -1.47798  ,   3.47361 ,   3.5714   };
const Double_t kThetaPar2[6] = { -0.666294 , -0.616567 , -0.673602 , -0.647113 ,  -0.34459 ,  -0.398458 };
const Double_t kThetaPar3[6] = {  5.73077  ,  5.51799  ,  8.05224  ,  7.74737  ,   8.45226 ,   9.54265  };
const Double_t kThetaPar4[6] = { 10.4976   , 14.0557   , 15.2178   , 16.7291   , -63.4556  , -22.649    };
const Double_t kThetaPar5[6] = { -1.13254  , -1.16189  , -2.08386  , -1.79939  ,  -3.3791  ,  -1.89746  };

// For parameter 0 of the FidPhiMin calculation for electron
const Double_t kFidPar0Low0[6] = {  25      ,  25        ,  25       ,  24.6345  ,  23.4731  ,  24.8599  };
const Double_t kFidPar0Low1[6] = { -12      , -12        , -12       , -12       , -12       , -12       };
const Double_t kFidPar0Low2[6] = {   0.5605 ,   0.714261 ,  0.616788 ,   0.62982 ,   1.84236 ,   1.00513 };
const Double_t kFidPar0Low3[6] = {  4.4     ,  4.4       ,  4.4      ,   4.4     ,   4.4     ,   4.4     };

// For parameter 1 of the FidPhiMin calculation for electron
const Double_t kFidPar1Low0[6] = {  2.1945   ,  4        ,  3.3352  ,  2.22769   ,  1.63143   ,  3.19807  };
const Double_t kFidPar1Low1[6] = {  1.51417  ,  1.56882  ,  2       ,  2         ,  1.90179   ,  0.173168 };
const Double_t kFidPar1Low2[6] = { -0.354081 , -2        , -2       , -0.760895  , -0.213751  , -0.1      };
const Double_t kFidPar1Low3[6] = {  0.5      ,  0.5      ,  1.01681 ,  1.31808   ,  0.786844  ,  1.6      };

// For parameter 0 of the FidPhiMax calculation for electron
const Double_t kFidPar0High0[6] = { 25       ,  25        ,  25        ,  25        ,  23.7067  ,  25       };
const Double_t kFidPar0High1[6] = { -8       , -10.3277   , -12        , -11.3361   , -12       , -11.4641  };
const Double_t kFidPar0High2[6] = {  0.479446 ,  0.380908 ,   0.675835 ,   0.636018 ,   2.92146 ,   0.55553 };
const Double_t kFidPar0High3[6] = {  4.8      ,  4.79964  ,   4.4      ,   4.4815   ,   4.4     ,   4.41327 };

// For parameter 1 of the FidPhiMax calculation for electron
const Double_t kFidPar1High0[6] = {  3.57349 ,  3.02279  ,  2.02102 ,  3.1948   ,  3.0934   ,  2.48828 };
const Double_t kFidPar1High1[6] = {  2       ,  0.966175 ,  2       ,  0.192701 ,  0.821726 ,  2       };
const Double_t kFidPar1High2[6] = { -2       , -2        , -1.70021 , -1.27578  , -0.233492 , -2       };
const Double_t kFidPar1High3[6] = {  0.5     ,  0.527823 ,  0.68655 ,  1.6      ,  1.6      ,  0.70261 };


bool accept_electron(TVector3 p)
{
  double mom = p.Mag();
  double theta = p.Theta() * 180./M_PI;
  double phi = p.Phi() * 180./M_PI;
  if (phi < -30.) phi+= 360.;

 // cout << "Phi in acceptance: " << phi << endl;

  int sector = (phi+30.)/60.;

  double minTheta = theta_min(mom,kThetaPar0[sector],
			      kThetaPar1[sector],
			      kThetaPar2[sector],
			      kThetaPar3[sector],
			      kThetaPar4[sector],
			      kThetaPar5[sector]);

  // Double check theta
  /* if (theta < minTheta)
    return false;
*/

  double phiCentral = 60.*sector;


  double aLow = a(mom,kFidPar0Low0[sector],
		  kFidPar0Low1[sector],
		  kFidPar0Low2[sector],
		  kFidPar0Low3[sector]);

  double aHigh = a(mom,kFidPar0High0[sector],
		  kFidPar0High1[sector],
		  kFidPar0High2[sector],
		  kFidPar0High3[sector]);

  double bLow = b(mom,kFidPar1Low0[sector],
		  kFidPar1Low1[sector],
		  kFidPar1Low2[sector],
		  kFidPar1Low3[sector]);

  double bHigh = b(mom,kFidPar1High0[sector],
		  kFidPar1High1[sector],
		  kFidPar1High2[sector],
		  kFidPar1High3[sector]);
  /*
  double aLow = a(mom,eFidAp0[sector][0],
		  eFidAp1[sector][0],
		  eFidAp2[sector][0],
		  eFidAp3[sector][0]);
  double aHigh = a(mom,eFidAp0[sector][1],
		  eFidAp1[sector][1],
		  eFidAp2[sector][1],
		  eFidAp3[sector][1]);
  double bLow = a(mom,eFidBp0[sector][0],
		  eFidBp1[sector][0],
		  eFidBp2[sector][0],
		  eFidBp3[sector][0]);
  double bHigh = a(mom,eFidBp0[sector][1],
		  eFidBp1[sector][1],
		  eFidBp2[sector][1],
		  eFidBp3[sector][1]);
  */

  double deltaPhiLow = deltaPhi(theta, aLow, bLow, minTheta);
  double deltaPhiHigh = deltaPhi(theta, aHigh, bHigh, minTheta);

  //std::cout << mom << " " << theta << " " << phi << " " << sector << " " << aLow << " " << aHigh << "\n";

  if ((phi < phiCentral + deltaPhiHigh) && (phi > phiCentral - deltaPhiLow))
    return true;

  return false;
}

// ############## NEUTRON FIDUCIAL CUT FOR ECAL #################

bool accept_neutron(TVector3 pm){
        double vtx_z = -22.25;
	double dist = 10.0;
        double r, pm_phi, x_rot, y_rot, angle, pm_theta, third_theta, d_intercept;

        const double ec_theta =  65.*M_PI/180.; // degrees
        const double z0       =563.1; // cm

        // Using law of sines to determine distance from vertex to EC
        pm_theta    = pm.Theta();
        third_theta = (M_PI - ec_theta - pm_theta);

        r = (z0-vtx_z)*sin(ec_theta)/sin(third_theta);

        double Xmiss = r*sin(pm.Theta())*cos(pm.Phi());
        double Ymiss = r*sin(pm.Theta())*sin(pm.Phi());

        pm_phi = 180./M_PI*pm.Phi();
        if (pm_phi < -30.) pm_phi += 360.;
        int sec = (int)(pm_phi+30)/60;
        if (sec>5) sec = 5;
        if (sec<0) sec = 0;

        angle = M_PI/180.*60.*(sec);

        x_rot = Xmiss*cos(angle) + Ymiss*sin(angle);
        y_rot = Xmiss*sin(angle) - Ymiss*cos(angle);

        d_intercept = dist/cos(atan(1.73));

        return((x_rot<390.-dist)&&(x_rot>1.73*fabs(y_rot)+55.+d_intercept));
}


// ############## NEUTRON FIDUCIAL CUT FOR TOF #################

bool accept_neutron_tof(TVector3 p)
{
  // This is our fiducial limit, 10 cm from the edge of each bar.
  const double distFromEdge=10.;

  // angle in degrees of the bar centers
  double thetaBars[57]={8.6,10.3,11.9,13.6,15.3,17,18.8,20.5,22.3,24,25.8,27.6,29.3,31.1,32.8,34.5,36.2,37.9,39.6,
			41.2,42.8,44.4,45.9,47.4,49.6,51.9,54.3,56.8,59.4,62,64.7,67.4,70.2,72.9,75.7,78.2,80.8,83.5,
			86.3,89.3,92.2,95.3,98.4,101.6,104.8,108.8,112,114.9,118.1,121.4,124.7,128,131.4,134.2,136.5,138.8,141};

  // width (length?) of the bars in cm
  double width[57]={32.3,48.1,64,79.8,95.7,106.6,122.4,138.3,154.1,170,185.8,201.7,217.6,233.4,249.3,265.1,281,296.8,312.7,
		    328.5,344.4,360.2,376.1,371.3,378.2,385,391.9,398.7,405.6,412.5,419.3,426.2,433,439.9,445.1,439.3,433.6,427.8,
		    422,416.3,410.5,404.8,399,393.3,387.5,380.1,363.5,347,330.4,313.9,297.3,280.8,264.2,246.8,235.4,224,212.7};

  // Distance to the target center in cm
  double distance[57]={513,509,503,500,498,496,495,494,493,493,493,494,495,496,498,501,504,507,511,
		       515,519,524,529,514,504,495,487,480,473,468,463,460,458,457,457,446,437,428,
		       421,414,409,405,402,400,399,402,395,389,384,381,378,377,378,379,380,383,385};


  double thetaDeg = p.Theta()*180./M_PI;
  if (thetaDeg<9 || thetaDeg>140)
    return false;

  // Sanitize phi to [-30,330]
  double phiDeg = p.Phi()*180./M_PI;
  if (phiDeg <-30.) phiDeg += 360.;

  // Get the value of phi within the relevant sector
  int sector = phiDeg/60;
  double PhiOneSector = phiDeg - 60.*sector;

  // Figure out the nearest bars
  int dsBar, usBar;
  for (dsBar=0, usBar=1; usBar < 57 ; dsBar++, usBar++)
    {
      if ((thetaBars[dsBar] < thetaDeg) && (thetaBars[usBar] > thetaDeg))
	break;
    }

  // Make a linear interpolation between the nearest bars
  double dsPhiMax = 180./M_PI*atan((width[dsBar] - distFromEdge)/(sin(thetaBars[dsBar]*M_PI/180.)*distance[dsBar]));
  double usPhiMax = 180./M_PI*atan((width[usBar] - distFromEdge)/(sin(thetaBars[usBar]*M_PI/180.)*distance[usBar]));
  double phiMax = dsPhiMax + (usPhiMax-dsPhiMax) * (thetaDeg - thetaBars[dsBar])/(thetaBars[usBar] - thetaBars[dsBar]);

  // Test if our neutron is within the phi bound
  if (fabs(PhiOneSector) < phiMax)
    return true;

  return false;
}



void do_SXC(int &lead_type, int &rec_type, double r)
{

  double pPP2NP = 0.041;
  double pPP2PN = 0.048;
  double pPP2NN = 0.0029;

  double pPN2NN = 0.035;
  double pPN2PP = 0.041;
  double pPN2NP = 0.0021;

  double pNP2PP = 0.035;
  double pNP2NN = 0.041;
  double pNP2PN = 0.0021;

  double pNN2PN = 0.041;
  double pNN2NP = 0.048;
  double pNN2PP = 0.0029;

  int pCode = 2212;
  int nCode = 2112;

  // Now we do a whole bunch of tests
  if ((lead_type==pCode) && (rec_type==pCode))
    {
      if (r < pPP2NN)
	{
	  lead_type=nCode;
	  rec_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP)
	{
	  lead_type=nCode;
	}
      else if (r < pPP2NN + pPP2NP + pPP2PN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==pCode))
    {
      if (r < pNP2PN)
	{
	  lead_type=pCode;
	  rec_type=nCode;
	}
      else if (r < pNP2PN + pNP2PP)
	{
	  lead_type=pCode;
	}
      else if (r < pNP2PN + pNP2PP + pNP2NN)
	{
	  rec_type=nCode;
	}
    }
  else if ((lead_type==pCode) && (rec_type==nCode))
    {
      if (r < pPN2NP)
	{
	  lead_type=nCode;
	  rec_type=pCode;
	}
      else if (r < pPN2NP + pPN2NN)
	{
	  lead_type=nCode;
	}
      else if (r < pPN2NP + pPN2NN + pPN2PP)
	{
	  rec_type=pCode;
	}
    }
  else if ((lead_type==nCode) && (rec_type==nCode))
    {
      if (r < pNN2PP)
	{
	  lead_type=pCode;
	  rec_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN)
	{
	  lead_type=pCode;
	}
      else if (r < pNN2PP + pNN2PN + pNN2NP)
	{
	  rec_type=pCode;
	}
    }
  else
    {
      std::cerr << "Invalid nucleon codes. Check and fix. Exiting\n\n\n";
    }
}


void Loop(TString filename1, TString filename2){
	using namespace genie;


      //  TString OutFileName = "/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_plead_500mil_f.root";
        TString OutFileName = "/pnfs/genie/scratch/users/nwright/devel_gaussian_q2/hists_emqe_gcf_G18_10a_eppstudy_plead_500mil_f.root";
        TFile *f_output = new TFile(OutFileName,"RECREATE");

        AccMap myAccMap("map_eg2_adin.root");

        TString filename[2] = {filename1,filename2};


        double mU = 0.9314941024;
        double reso_p = 0.01; // smearing for the proton
        double reso_e = 0.0125; // smearing for the electrons


        bool doFS = true;
        bool pn_observable = true;
        bool data_bins = true;
        bool physics_cuts_first = false;
        bool cx = false;
        bool smearing = true;

        int nevents = 0;

        int PNPairs = 0, PPPairs = 0, NNPairs = 0, NPPairs = 0;
        int SingleNucleon = 0, PairNucleons = 0;


        int bins = 50;


        int pLep_bins = bins;
        int thl_bins = bins;
        int pn_bins = bins;
        int pMiss_bins = bins;
        int En_bins = bins;
        int pLead_bins = bins;
        int pRec_bins = bins;
        int pcom_bins = bins;
        int xb_bins = 30;
        int q2_bins = 36;

        double pn_max = 1;
        double pn_min = 0;
        double pLead_min = 0;
        double pLead_max = 4;
        double pRec_min = 0;
        double pRec_max = 1;
        double xb_max = 2.;
        double xb_min = 1.2;
        double q2_min = 1.;
        double q2_max = 5.;
        double pMiss_max = 1;
        double pMiss_min = 0;

        if(data_bins){
                pLep_bins = 40;
                thl_bins = 60;
                pn_bins = 30;
                pMiss_bins = 9;
                pLead_bins = 40;
                pRec_bins = 17;
                xb_bins = 40;


                pn_min = .4;
                pLead_min = .4;
                pLead_max = 2.4;
                pRec_min = .35;
                pRec_max = 1.2;
                xb_max = 2.0;
        }

        vector<TString> cut;
        if(physics_cuts_first == false) {
                cut.push_back("");
                cut.push_back("_acceptance");
                cut.push_back("_bkx");
                cut.push_back("_thpq");
                cut.push_back("_pq");
                cut.push_back("_pmq");
                cut.push_back("_ep");
                cut.push_back("_epp");
        } else {
                cut.push_back("");
                cut.push_back("_bkx");
                cut.push_back("_thpq");
                cut.push_back("_pq");
                cut.push_back("_pmq");
                cut.push_back("_ep");
                cut.push_back("_epp");
                cut.push_back("_acceptance");

        }


       // cut.push_back("_thl");

        vector<TH1F *> h_list;

        TH1F* pLeadH[cut.size()];
        cout << 2 << endl;
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pLead"+cut[k]+"_hist";
		pLeadH[k] = new TH1F(PlotName,"",pLead_bins,pLead_min,pLead_max);
                h_list.push_back(pLeadH[k]);
        }

        TH1F* pLeadzH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pLeadz"+cut[k]+"_hist";
		pLeadzH[k] = new TH1F(PlotName,"",bins,-.5,4.);
        }

        TH1F* ELeadH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                cout << 6 << endl;
                TString PlotName = "ELead"+cut[k]+"_hist";
		ELeadH[k] = new TH1F(PlotName,"",bins,0,3.);
        }

        TH1F* p4LeadH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "p4Lead"+cut[k]+"_hist";
		p4LeadH[k] = new TH1F(PlotName,"",bins,0.5,1.);
        }

        TH1F* pRecH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pRec"+cut[k]+"_hist";
		pRecH[k] = new TH1F(PlotName,"",pRec_bins,pRec_min,pRec_max);
                h_list.push_back(pRecH[k]);

        }

        TH1F* p4RecH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "p4Rec"+cut[k]+"_hist";
		p4RecH[k] = new TH1F(PlotName,"",bins,0.,1.);
        }

        TH1F* ERecH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "ERec"+cut[k]+"_hist";
		ERecH[k] = new TH1F(PlotName,"",bins,.5,1.4);
        }

        TH1F* ELepH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "ELep"+cut[k]+"_hist";
		ELepH[k] = new TH1F(PlotName,"",bins*2,1.5,6.);
        }

        TH1F* pLepH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pLep"+cut[k]+"_hist";
		pLepH[k] = new TH1F(PlotName,"",pLep_bins,3,5.);
                h_list.push_back(pLepH[k]);

        }



        TH1F* pMissH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pMiss"+cut[k]+"_hist";
		pMissH[k] = new TH1F(PlotName,"",pMiss_bins,pMiss_min,pMiss_max);
                h_list.push_back(pMissH[k]);

        }

        TH1F* pMiss_eppH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pMiss_epp"+cut[k]+"_hist";
		pMiss_eppH[k] = new TH1F(PlotName,"",pn_bins,pn_min,pn_max);
                h_list.push_back(pMiss_eppH[k]);

        }


        TH1F* EMissH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "EMiss"+cut[k]+"_hist";
		EMissH[k] = new TH1F(PlotName,"",bins,-1.9,-1.8);
        }

        TH1F* theta_lrH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "theta_lr"+cut[k]+"_hist";
		theta_lrH[k] = new TH1F(PlotName,"",bins,0,180);
        }


    /*    TH1F* ptmiss_src_hist = new TH1F("ptmiss_src_hist", " ", bins, 0.,1.);

        TH1F* pmissm_src_hist = new TH1F("pmissm_src_hist", " ", bins, 1.,3.);

        TH1F* pmissz_src_hist = new TH1F("pmissz_src_hist", " ", bins, -.6,.6);
*/

        TH1F* pcomH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom"+cut[k]+"_hist";
		pcomH[k] = new TH1F(PlotName,"",bins,0,1);
                h_list.push_back(pcomH[k]);

        }

        TH1F* pcom_qH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom_q"+cut[k]+"_hist";
		pcom_qH[k] = new TH1F(PlotName,"",bins,-.85,.85);
        }

        TH1F* pcom_qtH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom_qt"+cut[k]+"_hist";
		pcom_qtH[k] = new TH1F(PlotName,"",bins,-.85,.85);
        }


        TH1F* pcom_xH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom_x"+cut[k]+"_hist";
		pcom_xH[k] = new TH1F(PlotName,"",bins,-.5,.5);
                h_list.push_back(pcom_xH[k]);

        }

        TH1F* pcom_yH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom_y"+cut[k]+"_hist";
		pcom_yH[k] = new TH1F(PlotName,"",bins,-.5,.5);
                h_list.push_back(pcom_yH[k]);

        }

        TH1F* pcom_zH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pcom_z"+cut[k]+"_hist";
		pcom_zH[k] = new TH1F(PlotName,"",bins,-.5,.5);
                h_list.push_back(pcom_zH[k]);

        }

        TH1F* prelH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "prel"+cut[k]+"_hist";
		prelH[k] = new TH1F(PlotName,"",bins,0,1.);
                h_list.push_back(prelH[k]);

        }

        TH1F* prel_xH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "prel_x"+cut[k]+"_hist";
		prel_xH[k] = new TH1F(PlotName,"",bins,-.9,.9);
                h_list.push_back(prel_xH[k]);

        }

        TH1F* prel_yH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "prel_y"+cut[k]+"_hist";
		prel_yH[k] = new TH1F(PlotName,"",bins,-.9,.9);
                h_list.push_back(prel_yH[k]);

        }

        TH1F* prel_zH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "prel_z"+cut[k]+"_hist";
		prel_zH[k] = new TH1F(PlotName,"",bins,-.9,.9);
                h_list.push_back(prel_zH[k]);

        }

        TH1F* proton_multH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "proton_mult"+cut[k]+"_hist";
		proton_multH[k] = new TH1F(PlotName,"",3,0,3.);
        }

        TH1F* q2H[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "q2"+cut[k]+"_hist";
		q2H[k] = new TH1F(PlotName,"",q2_bins,q2_min,q2_max);
                h_list.push_back(q2H[k]);

        }

        TH1F* bk_xH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                double xb_low = 0;
                if(k > 0){xb_low = xb_min;}
                TString PlotName = "bk_x"+cut[k]+"_hist";
		bk_xH[k] = new TH1F(PlotName,"",xb_bins,xb_low,xb_max);
                h_list.push_back(bk_xH[k]);

        }

        TH1F* theta_pqH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "theta_pq"+cut[k]+"_hist";
		theta_pqH[k] = new TH1F(PlotName,"",bins,0,80);
        }

        TH1F* pqH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pq"+cut[k]+"_hist";
		pqH[k] = new TH1F(PlotName,"",bins,0,1.1);
        }

        TH1F* qH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "q"+cut[k]+"_hist";
		qH[k] = new TH1F(PlotName,"",bins,0,5.1);
        }

        TH1F* qxH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "qx"+cut[k]+"_hist";
		qxH[k] = new TH1F(PlotName,"",bins,-2,2);
        }

        TH1F* qyH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "qy"+cut[k]+"_hist";
		qyH[k] = new TH1F(PlotName,"",bins,-2,2);
        }

        TH1F* qzH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "qz"+cut[k]+"_hist";
		qzH[k] = new TH1F(PlotName,"",bins,0,3);
        }



        TH1F* cthlH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cthl"+cut[k]+"_hist";
		cthlH[k] = new TH1F(PlotName,"",bins,.9,1.);
                h_list.push_back(cthlH[k]);

        }

        TH1F* cthrH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cthr"+cut[k]+"_hist";
		cthrH[k] = new TH1F(PlotName,"",bins,-1,1.);
        }

        TH1F* cthsH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cths"+cut[k]+"_hist";
		cthsH[k] = new TH1F(PlotName,"",bins,-1,1.);
        }

        TH1F* SRC_AngleH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "SRC_angle"+cut[k]+"_hist";
		SRC_AngleH[k] = new TH1F(PlotName,"",bins,0,181);
        }

        TH1F* thlH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "thl"+cut[k]+"_hist";
		thlH[k] = new TH1F(PlotName,"",thl_bins,10,40);
                h_list.push_back(thlH[k]);


        }

        TH1F* pmqH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pmq"+cut[k]+"_hist";
		if(cut[k] != "_epp") pmqH[k] = new TH1F(PlotName,"",40,100,180);
                else pmqH[k] = new TH1F(PlotName,"",20,100,180);
                h_list.push_back(pmqH[k]);

        }

        TH1F* pmiss_vqH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                double pmiss_vq_max = 1;
                if(cut[k] == "_bkx"){pmiss_vq_max = 0;}
                TString PlotName = "pmiss_vq"+cut[k]+"_hist";
		pmiss_vqH[k] = new TH1F(PlotName,"",bins,-.6,pmiss_vq_max);
        }

        TH1F* EmissMPmissH[cut.size()];
        for(int k = 0; k< cut.size(); k++){

                TString PlotName = "EmissMPmiss"+cut[k]+"_hist";
		EmissMPmissH[k] = new TH1F(PlotName,"",bins*4,-10,10);
        }

        TH1F* Emiss2MPmiss2H[cut.size()];
        for(int k = 0; k< cut.size(); k++){

                TString PlotName = "Emiss2MPmiss2"+cut[k]+"_hist";
		Emiss2MPmiss2H[k] = new TH1F(PlotName,"",bins,-1,1.5);
        }

        TH1F* wH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "w"+cut[k]+"_hist";
		wH[k] = new TH1F(PlotName,"",bins,0,1.8);
        }

        TH1F* omegaH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "omega"+cut[k]+"_hist";
		omegaH[k] = new TH1F(PlotName,"",bins,0,5.1);
        }

        TH1F* EcalH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "Ecal"+cut[k]+"_hist";
		EcalH[k] = new TH1F(PlotName,"",bins,4.5,5.5);
        }

        TH1F* Ecal_prH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "Ecal_pr"+cut[k]+"_hist";
		Ecal_prH[k] = new TH1F(PlotName,"",bins,4.5,5.5);
        }


        TH1F* cth_pcomH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cth_pcom"+cut[k]+"_hist";
		cth_pcomH[k] = new TH1F(PlotName,"",bins,-1,1.);
                h_list.push_back(cth_pcomH[k]);

        }

        TH1F* cth_prelH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cth_prel"+cut[k]+"_hist";
		cth_prelH[k] = new TH1F(PlotName,"",bins,-1,1.);
                h_list.push_back(cth_prelH[k]);

        }

        TH1F* phi_pcomH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "phi_pcom"+cut[k]+"_hist";
		phi_pcomH[k] = new TH1F(PlotName,"",30,-180,180);
        }

        TH1F* phi_prelH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "phi_prel"+cut[k]+"_hist";
		phi_prelH[k] = new TH1F(PlotName,"",30,-180,180);
        }

        TH1F* th_pmiqH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "th_pmiq"+cut[k]+"_hist";
		th_pmiqH[k] = new TH1F(PlotName,"",bins,0,181);
        }

        TH1F* pnH[cut.size()];
        cout << 2 << endl;
        for(int k = 0; k< cut.size(); k++){
                cout << 3 << endl;
                TString PlotName = "pn"+cut[k]+"_hist";
                cout << 4 << endl;
		pnH[k] = new TH1F(PlotName,"",pn_bins,pn_min,1.);
                cout << 5 << endl;
                h_list.push_back(pnH[k]);

        }

        TH1F* pn_genieH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pn_genie"+cut[k]+"_hist";
		pn_genieH[k] = new TH1F(PlotName,"",bins,0,1.);
        }

        TH1F* pnxH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pnx"+cut[k]+"_hist";
		pnxH[k] = new TH1F(PlotName,"",bins,-1.,1.);
        }

        TH1F* pnyH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pny"+cut[k]+"_hist";
		pnyH[k] = new TH1F(PlotName,"",bins,-1.,1.);
        }

        TH1F* pnzH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "pnz"+cut[k]+"_hist";
		pnzH[k] = new TH1F(PlotName,"",bins,-1.,1.);
        }

        TH1F* MmH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "Mm"+cut[k]+"_hist";
		MmH[k] = new TH1F(PlotName,"",bins,.5,1.5);
        }

        TH1F* pnVarH[cut.size()];
        cout << 2 << endl;
        for(int k = 0; k< cut.size(); k++){
                cout << 3 << endl;
                TString PlotName = "pnVar"+cut[k]+"_hist";
                cout << 4 << endl;
		pnVarH[k] = new TH1F(PlotName,"",bins,0,3.);
                cout << 5 << endl;
        }

        TH1F* p4nH[cut.size()];
        cout << 2 << endl;
        for(int k = 0; k< cut.size(); k++){
                cout << 3 << endl;
                TString PlotName = "p4n"+cut[k]+"_hist";
                cout << 4 << endl;
		p4nH[k] = new TH1F(PlotName,"",bins,0,2.);
                cout << 5 << endl;
        }

        TH1F* EnH[cut.size()];
        cout << 2 << endl;
        for(int k = 0; k< cut.size(); k++){
                cout << 3 << endl;
                TString PlotName = "En"+cut[k]+"_hist";
                cout << 4 << endl;
		EnH[k] = new TH1F(PlotName,"",bins,0,1.);
                cout << 5 << endl;
        }

        TH2F* EmissvPmissH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "EmissvPmiss"+cut[k]+"_hist";
		EmissvPmissH[k] = new TH2F(PlotName,"",bins,0,1.,bins,.5,.95);
        }

        TH2F* cthlvElepH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "cthlvElep"+cut[k]+"_hist";
		cthlvElepH[k] = new TH2F(PlotName,"",bins,1.5,6.,bins,.9,1);
        }

        TH2F* pmissvxBH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                 double xb_low = 0;
                if(k > 0){xb_low = 1.2;}
                TString PlotName = "pmissvxB"+cut[k]+"_hist";
		pmissvxBH[k] = new TH2F(PlotName,"",bins,0,2,bins,0,1.);
        }

        TH2F* Q2vpnH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "Q2vpn"+cut[k]+"_hist";
		Q2vpnH[k] = new TH2F(PlotName,"",bins,0,1.,bins,.9,3);
        }

        TH2F* thetavphipH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "thetavphip"+cut[k]+"_hist";
		thetavphipH[k] = new TH2F(PlotName,"",bins,0,360,bins,0,100);
        }

        TH2F* thetavphilepH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "thetavphilep"+cut[k]+"_hist";
		thetavphilepH[k] = new TH2F(PlotName,"",bins,0,360,bins,0,80);
        }

        TH2F* thetavppH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "thetavpp"+cut[k]+"_hist";
		thetavppH[k] = new TH2F(PlotName,"",bins,0,4,bins,0,180);
        }

        TH2F* thetavplepH[cut.size()];
        for(int k = 0; k< cut.size(); k++){
                TString PlotName = "thetavplep"+cut[k]+"_hist";
		thetavplepH[k] = new TH2F(PlotName,"",bins,1.5,5.2,bins,0,180);
        }


        TH1F* spectral_func = new TH1F("spectral_func","",bins*2,0,1.);
        //Event



        double m_p = .938;
        TRandom3 *gRandom = new TRandom3();

        cout << "before the file loop" << endl;

        int events_pass_epp = 0;
        int events_pass_ep = 0;



        for(int g = 0; g < 2; g++){
                // Open the GHEP/ROOT file
                TFile file(filename[g], "READ");
                // Get the tree header & print
                NtpMCTreeHeader * header =dynamic_cast<NtpMCTreeHeader*> (file.Get("header"));
        //  std::cout << *header << std::endl;
                // Get the GENIE GHEP tree and set its branch address
                TTree * tree = dynamic_cast<TTree*> (file.Get("gtree"));
                NtpMCEventRecord * mcrec = 0;
                tree->SetBranchAddress("gmcrec", &mcrec);

                double nentries = tree->GetEntries();
                cout << "entries: " << nentries << endl;

               // nentries = 30;
               // nentries = 5000000;

                if(g == 0 && cx == true) nentries = 10060000;
                if(g == 1 && cx == true) nentries -= 10000000;

                for(Long64_t i=0; i<nentries; i++){
                     /*   if(g == 0 && i == 10060000){
                                i += 7000000;
                        }
                       */
                        tree->GetEntry(i);
                        nevents++;

                        cout << "file: " << g+1 << endl;
                        if(g == 0) cout << i << endl;
                        else cout << i + nentries << endl;

                       // gRandom->SetSeed(i);
                        // print-out the event

                        if(mcrec == NULL){
                                mcrec->Clear();
                                continue;
                        }
                        EventRecord & event = *(mcrec->event);

                        const Interaction *summary = event.Summary();
                        const ProcessInfo &proc = summary->ProcInfo();
                        const Kinematics &kine = summary->Kine();

                        GHepParticle * nucleus = event.TargetNucleus();
                        TLorentzVector * nucleus_p4 = new TLorentzVector(0,0,0,0);
                        if(nucleus != NULL){nucleus_p4 = nucleus->P4();}

        //                cout << "Nucleus Mass: " << nucleus->Mass() << endl;


                        bool qelcc = proc.IsQuasiElastic() && proc.IsEM();


                        GHepParticle * lepton = event.FinalStatePrimaryLepton();
                        TLorentzVector * lepton_p4 = lepton->P4();
                        TVector3 lepton_p3 = lepton_p4->Vect();
                        double lepton_momentum = lepton_p3.Mag();


                        double cthl = lepton_p3.CosTheta();
                        double Ek = lepton_p4->E();

                        if(TMath::ACos(lepton_p3.CosTheta())*180/TMath::Pi() < 13){
                                mcrec->Clear();
                                continue;
                        }


                        double lep_mass;
                        if(proc.IsWeakCC()){lep_mass = 0.10566;}
                        if(proc.IsEM()){lep_mass = 0.000511;}



                        double theta_l = TMath::ACos(cthl)*180/TMath::Pi();
                        if (theta_l <= 0) { theta_l += 180;}
                        if (theta_l > 180) { theta_l -= 180;}


                        GHepParticle * probe = event.Probe();
                        TLorentzVector * probe_p4 = probe->P4();


                        vector<GHepParticle> initial_hadrons;
                        vector<GHepParticle> final_state_particles;

                        GHepParticle * remnant;

                        int p0_id = -3;
                        int p1_id = -3;


                        TVector3 pn;
                        GHepParticle * struck_nucleon;
                        GHepParticle * lead_nucleon;

                        double pn_genie = 0;
                        bool struck_isp = false;

                        for(int j = 0;j < event.GetEntries();j++){
                                GHepParticle *p = event.Particle(j);

                                int pdg = p->Pdg();
                                TVector3 p_mom = p->P4()->Vect();


                                if(p->Status() == 11){
                                        struck_nucleon = p;
                                        pn = struck_nucleon->P4()->Vect();
                                        pn_genie = pn.Mag();
                                        if(pn.Mag() > pn_max){pn_max = pn.Mag();}

                                        if(pdg == 2212){ struck_isp = true;}

                                }


                                if(p->Status() == kIStHadronInTheNucleus){

                                        if(pdg::IsNucleon(pdg)){
                                                if(p0_id == -3){p0_id = p->FirstDaughter();}
                                                else{p1_id = p->FirstDaughter();}
                                                initial_hadrons.push_back(*p);
                                                if(pdg::IsProton(pdg)){ pr_count++; }
                                        }
                                }

                                if(p->Status() == kIStStableFinalState){

                                        if(j == p0_id || j == p1_id){
                                                if(pdg::IsNucleon(pdg)){
                                                        final_state_particles.push_back(*p);
                                                }
                                        }

                                }

                                if(p->Status() == kIStFinalStateNuclearRemnant){
                                        remnant = p;
                                }

                        } // Looping over all particles

                        //Analysis



                        if(doFS){

                                if(final_state_particles.size() > 2){
                                        mcrec->Clear();
                                        continue;
                                }
                                initial_hadrons.clear();
                                for(int j = 0; j < final_state_particles.size(); j++){
                                        initial_hadrons.push_back(final_state_particles[j]);
                                }
                        }


                        int ni = initial_hadrons.size();

                        if(ni == 1){
                        spectral_func->Fill(pn.Mag());
                        }
                        if(ni != 2){
                                mcrec->Clear();
                                continue;
                        }



                        if(ni == 2){


                                PairNucleons++;


                                bool is_hit_p  = pdg::IsProton(initial_hadrons[1].Pdg());
                                bool is_rec_p  = pdg::IsProton(initial_hadrons[0].Pdg());
                                int Z = nucleus->Z();
                                if((is_hit_p && !is_rec_p) || (!is_hit_p && is_rec_p)){ Z = nucleus->Z()-1; }
                                if(is_hit_p && is_rec_p){ Z = nucleus->Z()-2; }
                                if(!is_hit_p && !is_rec_p){ Z = nucleus->Z(); }
                                int A = nucleus->A() - 2;

                                /*TParticlePDG * remn_nucleus_p = 0;
                                int ipdgc = pdg::IonPdgCode(A, Z);
                                remn_nucleus_p = PDGLibrary::Instance()->Find(ipdgc);

                                double remn_nucleus_mass = remn_nucleus_p->Mass();
                                */
                                double remn_nucleus_mass = 10.016853 * mU - 6*lep_mass; //agree with generator

                                int np = 0;
                                if(is_hit_p && is_rec_p){
                                        PPPairs++;
                                        np = 2;
                                }
                                if(is_hit_p && !is_rec_p){
                                        np = 1;
                                        PNPairs++;
                                }
                                if(!is_hit_p && is_rec_p){
                                        NPPairs++;
                                        np = 1;
                                }
                                if(!is_hit_p && !is_rec_p){
                                        NNPairs++;
                                }

                                TLorentzVector * recoil_p4 = initial_hadrons[0].P4();
                                TVector3 recoil_p3(recoil_p4->X(),recoil_p4->Y(),recoil_p4->Z());

                                TLorentzVector * leading_p4 = initial_hadrons[1].P4();
                                TVector3 leading_p3(leading_p4->X(),leading_p4->Y(),leading_p4->Z());

                                if(smearing){
                                        double SmearedPe = gRandom->Gaus(lepton_p3.Mag(),reso_e*lepton_p3.Mag());
                                        lepton_p3.SetMagThetaPhi(SmearedPe,lepton_p3.Theta(),lepton_p3.Phi());

                                        double SmearedPLead = gRandom->Gaus(leading_p3.Mag(),reso_p*leading_p3.Mag());
                                        leading_p3.SetMagThetaPhi(SmearedPLead,leading_p3.Theta(),leading_p3.Phi());

                                        if(is_rec_p){
                                                double SmearedPRec = gRandom->Gaus(recoil_p3.Mag(),reso_p*recoil_p3.Mag());
                                                recoil_p3.SetMagThetaPhi(SmearedPRec,recoil_p3.Theta(),recoil_p3.Phi());
                                        }
                                }



                                TVector3 vq = probe_p4->Vect() - lepton_p3;
                                double omega = probe_p4->E() - Ek; // Energy transfer
                                TLorentzVector Q2Vec(vq.X(),vq.Y(),vq.Z(),omega);
                                double QSq = Q2Vec.Mag2(); // Q-Squared
                                QSq = fabs(QSq);

                                double bk_x = QSq/(2*m_p*omega);
                                bk_x = fabs(bk_x);


                                if(pn_observable){
                                        pn = leading_p4->Vect() - vq;
                                        if(pn.Mag() < .25){
                                                mcrec->Clear();
                                                continue;
                                        }
                                }

                                double En = struck_nucleon->E();
                                TLorentzVector pn4(pn,En);

                                TVector3 pcom = recoil_p3 + pn;
                                TVector3 prel = .5*(recoil_p3 - pn);


                                spectral_func->Fill(prel.Mag());


                                double phi_pcom = pcom.Phi()*180/TMath::Pi();
                                double phi_prel = prel.Phi()*180/TMath::Pi();


                                double KE_1 = initial_hadrons[0].KinE();
                                double KE_2 = initial_hadrons[1].KinE();

                                double E_1 = KE_1 + .04;
                                double E_2 = KE_2 + .04;
                                double E_l = lepton->Energy();

                                double E_Cal = E_1 + E_2 + E_l;

                                if(initial_hadrons[0].Pdg() != 2212){
                                        E_Cal -= KE_1;
                                }
                                if(initial_hadrons[1].Pdg() != 2212){
                                        E_Cal -= KE_2;
                                }


                                double th_pq = vq.Angle(leading_p4->Vect())*180/TMath::Pi();
                                if (th_pq <= 0) { th_pq += 180;}
                                if (th_pq > 180) { th_pq -= 180;}

                                double pq_ratio = leading_p4->Vect().Mag()/vq.Mag();
                                TVector3 pn = leading_p4->Vect() - probe_p4->Vect() + lepton_p4->Vect();


                                TLorentzVector nuc_m(0,0,0,m_p);
                                TLorentzVector Mm2Vec = Q2Vec + 2*nuc_m - *leading_p4;
                                double Mm = Mm2Vec.M();


                                TLorentzVector pnM(0,0,0,Mm);
                                TLorentzVector wVec = Q2Vec + pnM;
                                double w = wVec.Mag();


                                if(prel.Mag() < .25){
                                        mcrec->Clear();
                                        continue;
                                }



                                double E_miss = leading_p4->E() - omega;
                                TVector3 P_miss = leading_p3 - vq;

                                /* Old vars
                                double P_miss_z = P_miss.Z();
                                double P_miss_m = E_miss - P_miss_z;
                                P_miss_m = fabs(P_miss_m);


                                double nuc_mass = pn4.M();
                                double nuc_mass_2 = nuc_mass*nuc_mass;

                                double kmiss2 = nuc_mass_2*((P_miss.Perp2() + nuc_mass_2)/(P_miss_m*(2*nuc_mass - P_miss_m))) - nuc_mass_2;
                                kmiss2 = fabs(kmiss2);
                                double kmiss = sqrt(kmiss2);
                                */

                                double src_angle = pn.Angle(recoil_p3)*180./TMath::Pi();
                                if (src_angle <= 0) { src_angle += 180;}
                                if (src_angle > 180) { src_angle -= 180;}

                                double theta_lr = leading_p4->Vect().Angle(recoil_p3)*180./TMath::Pi();
                                if (theta_lr <= 0) { theta_lr += 180;}
                                if (theta_lr > 180) { theta_lr -= 180;}


                                double lep_phi = lepton_p3.Phi()*180/TMath::Pi() + 210;
                                if(lep_phi < 0){lep_phi += 360;}
                                if(lep_phi > 360){lep_phi -= 360;}

                                double lead_phi = leading_p3.Phi()*180/TMath::Pi() + 210;
                                if(lead_phi < 0){lead_phi += 360;}
                                if(lead_phi > 360){lead_phi -= 360;}


                                ///FILLING HISTOGRAMS
                                ///We aren't looking at these right now and my disk space is precious
                             /*   MmH[0]->Fill(Mm);

                                pmqH[0]->Fill(pn.Angle(vq)*180/TMath::Pi());

                                proton_multH[0]->Fill(np);
                                SRC_AngleH[0]->Fill(src_angle);

                                qxH[0]->Fill(vq.X());
                                qyH[0]->Fill(vq.Y());
                                qzH[0]->Fill(vq.Z());
                                theta_lrH[0]->Fill(theta_lr);
                                EnH[0]->Fill(En);
                                pqH[0]->Fill(pq_ratio);
                                pnH[0]->Fill(pn.Mag());
                                p4nH[0]->Fill(pn4.Mag());
                                pnVarH[0]->Fill(pn.Perp()/(En - pn.Z()));
                                pnxH[0]->Fill(pn.X());
                                pnyH[0]->Fill(pn.Y());
                                pnzH[0]->Fill(pn.Z());


                                theta_pqH[0]->Fill(th_pq);
                                cthlH[0]->Fill(cthl);
                                EcalH[0]->Fill(E_Cal);


                                prelH[0]->Fill(prel.Mag());
                                prel_xH[0]->Fill(prel.X());
                                prel_yH[0]->Fill(prel.Y());
                                prel_zH[0]->Fill(prel.Z());

                                pcom_xH[0]->Fill(pcom.X());
                                pcom_yH[0]->Fill(pcom.Y());
                                pcom_zH[0]->Fill(pcom.Z());
                                pcomH[0]->Fill(pcom.Mag());
                                pcom_qH[0]->Fill(pcom.Dot(vq));
                                pcom_qtH[0]->Fill(pcom.Dot(vq.Orthogonal()));



                                cth_pcomH[0]->Fill(pcom.CosTheta());
                                cth_prelH[0]->Fill(prel.CosTheta());



                                phi_pcomH[0]->Fill(phi_pcom);
                                phi_prelH[0]->Fill(phi_prel);


                                wH[0]->Fill(w);
                                qH[0]->Fill(vq.Mag());
                                q2H[0]->Fill(QSq);

                                bk_xH[0]->Fill(bk_x);

                                ELepH[0]->Fill(Ek);
                                thlH[0]->Fill(theta_l);

                                p4LeadH[0]->Fill(leading_p4->Mag());
                                p4RecH[0]->Fill(recoil_p4->Mag());
                                pLeadH[0]->Fill(leading_p4->Vect().Mag());
                                pLeadzH[0]->Fill(leading_p4->Vect().Z());

                                pRecH[0]->Fill(recoil_p3.Mag());
                                ELeadH[0]->Fill(leading_p4->E());
                                ERecH[0]->Fill(recoil_p4->E());

                                cthrH[0]->Fill(recoil_p3.CosTheta());
                                cthsH[0]->Fill(pn.CosTheta());

                                pmiss_vqH[0]->Fill(P_miss.Dot(vq)/vq.Mag());
                                pMissH[0]->Fill(P_miss.Mag());
                                EMissH[0]->Fill(E_miss);

                                EmissMPmissH[0]->Fill(E_miss - P_miss.Dot(vq)/vq.Mag());
                                Emiss2MPmiss2H[0]->Fill(E_miss*E_miss - (P_miss.Dot(vq)/vq.Mag())*(P_miss.Dot(vq)/vq.Mag()));

                                th_pmiqH[0]->Fill(vq.Angle(P_miss)*180/TMath::Pi());

                                cthlvElepH[0]->Fill(Ek,cthl);
                                Q2vpnH[0]->Fill(pn.Mag(),QSq);





                                thetavphilepH[0]->Fill(lep_phi,lepton_p3.Theta()*180/TMath::Pi());
                                thetavplepH[0]->Fill(lepton_p3.Mag(),lepton_p3.Theta()*180/TMath::Pi());

                                pLepH[0]->Fill(lepton_p3.Mag());

                                if(pdg::IsProton(initial_hadrons[1].Pdg())){
                                        thetavphipH[0]->Fill(lead_phi,leading_p3.Theta()*180/TMath::Pi());
                                        thetavppH[0]->Fill(leading_p3.Mag(),leading_p3.Theta()*180/TMath::Pi());
                                }

                                EmissvPmissH[0]->Fill(P_miss.Mag(),E_miss);
                                pMiss_eppH[0]->Fill(P_miss.Mag());
                                pmissvxBH[0]->Fill(bk_x,P_miss.Mag());

                                pn_genieH[0]->Fill(pn_genie);

                                */
                                bool pass_epp = true;
                                bool pass_acceptance = true;



                                double lepton_phi = lepton_p3.Phi();
                                double leading_phi = leading_p3.Phi();
                                double recoil_phi = recoil_p3.Phi();
                                lepton_p3.SetPhi(lepton_phi + TMath::Pi());
                                leading_p3.SetPhi(leading_phi + TMath::Pi());
                                recoil_p3.SetPhi(recoil_phi + TMath::Pi());

                                double acc_weight = myAccMap.accept(leading_p3);



                                /*
                                if(doFS){
                                        vector<GHepParticle> final_hadrons;
                                        for(int j = 0; j < initial_hadrons.size(); j++){
                                                if(j == 1 || j == 2){
                                                        final_hadrons.push_back(initial_hadrons[j]);
                                                        continue;
                                                }
                                                else{
                                                        if(initial_hadrons[j].Pdg() == 2212)
                                                                final_hadrons.push_back(initial_hadrons[j]);
                                                        }
                                        }
                                                initial_hadrons.clear();
                                                for(int j = 0; j < final_hadrons.size(); j++){
                                                        initial_hadrons.push_back(final_hadrons[j]);

                                                        TVector3 p_vect = final_hadrons[j].P4()->Vect();
                                                        double p_phi = p_vect.Phi();
                                                        p_vect.SetPhi(p_phi + TMath::Pi());

                                                        if(!accept_proton(p_vect))
                                                                pass_acceptance = false;
                                                }

                                }
                                */

                                if(!accept_electron(lepton_p3)){
                                        pass_acceptance = false;
                                }

                                int lead_type = initial_hadrons[1].Pdg();
                                int rec_type = initial_hadrons[0].Pdg();

                                if(cx){
                                        double r = gRandom->Rndm();
                                        do_SXC(lead_type,rec_type,r);
                                }

                                if(lead_type != 2212 || !accept_proton(leading_p3)){
                                        pass_acceptance = false;

                                }

                                if(rec_type != 2212 || !accept_proton(recoil_p3)){
                                        pass_epp = false;
                                }




                                if(acc_weight == 0){
                                        mcrec->Clear();
                                        continue;
                                }

                                lepton_p3.SetPhi(lepton_phi);
                                leading_p3.SetPhi(leading_phi);
                                recoil_p3.SetPhi(recoil_phi);


                                //PCM Rotation
                                TVector3 v1 = pn;
                                TVector3 v2 = recoil_p3;
                                double Pmiss_phi     = v1.Phi();
                                double Pmiss_theta   = v1.Theta();
                                //--------------------------
                                vq.RotateZ(-Pmiss_phi);
                                vq.RotateY(-Pmiss_theta);
                                double q_phi = vq.Phi();
                                vq.RotateZ(-q_phi);
                                //--------------------------
                                v2.RotateZ(-Pmiss_phi);
                                v2.RotateY(-Pmiss_theta);
                                v2.RotateZ(-q_phi);
                                //--------------------------
                                v1.RotateZ(-Pmiss_phi);
                                v1.RotateY(-Pmiss_theta);
                                v1.RotateZ(-q_phi);

                                pcom = v1 + v2;


                                if(!physics_cuts_first){
                                        for(int index = 1; index < 8; index++){

                                                if(index == 1 && !pass_acceptance)
                                                        break;
                                                if(index == 2 && bk_x < 1.2)
                                                        break;
                                                if(index == 3 && th_pq  > 25)
                                                        break;
                                                if(index == 4){
                                                        if(pq_ratio < .62 || pq_ratio > .96)
                                                                break;
                                                }
                                                if(index == 5){
                                                        if(pn.Mag() > 1 || pn.Mag() < .4)
                                                                break;
                                                        }
                                                if(index == 6 && Mm > 1.1)
                                                        break;
                                                else {events_pass_ep++;}
                                                if(index == 7){
                                                        if(!pass_epp || recoil_p3.Mag() < .35)
                                                                break;
                                                        else{
                                                                events_pass_epp++;
                                                                acc_weight*= myAccMap.accept(recoil_p3);
                                                        }
                                                }

                                                pmqH[index]->Fill(pn.Angle(vq)*180/TMath::Pi(),acc_weight);

                                                MmH[index]->Fill(Mm,acc_weight);

                                                proton_multH[index]->Fill(np,acc_weight);
                                                SRC_AngleH[index]->Fill(src_angle,acc_weight);

                                                pLeadzH[index]->Fill(leading_p4->Vect().Z(),acc_weight);

                                                qxH[index]->Fill(vq.X(),acc_weight);
                                                qyH[index]->Fill(vq.Y(),acc_weight);
                                                qzH[index]->Fill(vq.Z(),acc_weight);
                                                theta_lrH[index]->Fill(theta_lr,acc_weight);
                                                EnH[index]->Fill(En,acc_weight);
                                                pnH[index]->Fill(pn.Mag(),acc_weight);
                                                p4nH[index]->Fill(pn4.Mag(),acc_weight);

                                                pqH[index]->Fill(pq_ratio,acc_weight);

                                                pnxH[index]->Fill(pn.X(),acc_weight);
                                                pnyH[index]->Fill(pn.Y(),acc_weight);
                                                pnzH[index]->Fill(pn.Z(),acc_weight);

                                                theta_pqH[index]->Fill(th_pq,acc_weight);
                                                cthlH[index]->Fill(cthl,acc_weight);
                                                EcalH[index]->Fill(E_Cal,acc_weight);



                                                if(pass_epp){
                                                        prelH[index]->Fill(prel.Mag(),acc_weight);
                                                        prel_xH[index]->Fill(prel.X(),acc_weight);
                                                        prel_yH[index]->Fill(prel.Y(),acc_weight);
                                                        prel_zH[index]->Fill(prel.Z(),acc_weight);

                                                        pcom_xH[index]->Fill(pcom.X(),acc_weight);
                                                        pcom_yH[index]->Fill(pcom.Y(),acc_weight);
                                                        pcom_zH[index]->Fill(pcom.Z(),acc_weight);
                                                        pcomH[index]->Fill(pcom.Mag(),acc_weight);

                                                        pcom_zH[index]->Fill(pcom.Z(),acc_weight);
                                                        pcomH[index]->Fill(pcom.Mag(),acc_weight);
                                                        pcom_qH[index]->Fill(pcom.Dot(vq),acc_weight);
                                                        pcom_qtH[index]->Fill(pcom.Dot(vq.Orthogonal()),acc_weight);
                                                        pRecH[index]->Fill(recoil_p3.Mag(),acc_weight);

                                                        cth_pcomH[index]->Fill(pcom.CosTheta(),acc_weight);
                                                        cth_prelH[index]->Fill(prel.CosTheta(),acc_weight);



                                                        phi_pcomH[index]->Fill(phi_pcom,acc_weight);
                                                        phi_prelH[index]->Fill(phi_prel,acc_weight);

                                                        ERecH[index]->Fill(recoil_p4->E(),acc_weight);
                                                        p4RecH[index]->Fill(recoil_p4->Mag(),acc_weight);

                                                        pMiss_eppH[index]->Fill(P_miss.Mag(),acc_weight);

                                                }




                                                wH[index]->Fill(w,acc_weight);
                                                qH[index]->Fill(vq.Mag(),acc_weight);
                                                q2H[index]->Fill(QSq,acc_weight);

                                                bk_xH[index]->Fill(bk_x,acc_weight);

                                                ELepH[index]->Fill(Ek,acc_weight);
                                                thlH[index]->Fill(theta_l,acc_weight);

                                                p4LeadH[index]->Fill(leading_p4->Mag(),acc_weight);
                                                pLeadH[index]->Fill(leading_p4->Vect().Mag(),acc_weight);
                                                ELeadH[index]->Fill(leading_p4->E(),acc_weight);

                                                cthrH[index]->Fill(recoil_p3.CosTheta(),acc_weight);
                                                cthsH[index]->Fill(pn.CosTheta(),acc_weight);

                                                pmiss_vqH[index]->Fill(P_miss.Dot(vq)/vq.Mag(),acc_weight);
                                                pMissH[index]->Fill(P_miss.Mag(),acc_weight);
                                                EMissH[index]->Fill(E_miss,acc_weight);

                                                EmissMPmissH[index]->Fill(E_miss - P_miss.Dot(vq)/vq.Mag(),acc_weight);
                                                Emiss2MPmiss2H[index]->Fill(E_miss*E_miss - (P_miss.Dot(vq)/vq.Mag())*(P_miss.Dot(vq)/vq.Mag()),acc_weight);


                                                th_pmiqH[index]->Fill(vq.Angle(P_miss)*180/TMath::Pi(),acc_weight);

                                                cthlvElepH[index]->Fill(Ek,cthl,acc_weight);

                                                Q2vpnH[index]->Fill(pn.Mag(),QSq,acc_weight);

                                                thetavphilepH[index]->Fill(lep_phi,lepton_p3.Theta()*180/TMath::Pi(),acc_weight);
                                                thetavplepH[index]->Fill(lepton_p3.Mag(),lepton_p3.Theta()*180/TMath::Pi(),acc_weight);

                                                pLepH[index]->Fill(lepton_p3.Mag(),acc_weight);


                                                if(pdg::IsProton(initial_hadrons[1].Pdg())){
                                                        thetavphipH[index]->Fill(lead_phi,leading_p3.Theta()*180/TMath::Pi(),acc_weight);
                                                        thetavppH[index]->Fill(leading_p3.Mag(),leading_p3.Theta()*180/TMath::Pi(),acc_weight);
                                                }

                                                pn_genieH[index]->Fill(pn_genie,acc_weight);


                                                pmissvxBH[index]->Fill(bk_x,P_miss.Mag(),acc_weight);
                                                EmissvPmissH[index]->Fill(P_miss.Mag(),E_miss,acc_weight);

                                        }
                                }



                                if(physics_cuts_first){
                                        for(int index = 1; index < 8; index++){
                                                if(index == 1 && bk_x < 1.2)
                                                        break;
                                                if(index == 2 && th_pq  > 25)
                                                        break;
                                                if(index == 3){
                                                        if(pq_ratio < .62 || pq_ratio > .96)
                                                                break;
                                                }
                                                if(index == 4){
                                                        if(pn.Mag() > 1 || pn.Mag() < .4)
                                                                break;
                                                }
                                                if(index == 5 && Mm > 1.1)
                                                        break;
                                                if(index == 6 && !pass_acceptance)
                                                        break;
                                                if(index == 7){
                                                        if(!pass_epp || recoil_p3.Mag() < .35)
                                                                break;
                                                }

                                                MmH[index]->Fill(Mm);
                                                pmqH[index]->Fill(pn.Angle(vq)*180/TMath::Pi());

                                                proton_multH[index]->Fill(np);
                                                SRC_AngleH[index]->Fill(src_angle);

                                                qxH[index]->Fill(vq.X());
                                                qyH[index]->Fill(vq.Y());
                                                qzH[index]->Fill(vq.Z());
                                                theta_lrH[index]->Fill(theta_lr);
                                                EnH[index]->Fill(En);
                                                pqH[index]->Fill(pq_ratio);
                                                pnH[index]->Fill(pn.Mag());
                                                p4nH[index]->Fill(pn4.Mag());
                                                pnVarH[index]->Fill(pn.Perp()/(En - pn.Z()));
                                                pnxH[index]->Fill(pn.X());
                                                pnyH[index]->Fill(pn.Y());
                                                pnzH[index]->Fill(pn.Z());


                                                theta_pqH[index]->Fill(th_pq);
                                                cthlH[index]->Fill(cthl);
                                                EcalH[index]->Fill(E_Cal);


                                                prelH[index]->Fill(prel.Mag());
                                                prel_xH[index]->Fill(prel.X());
                                                prel_yH[index]->Fill(prel.Y());
                                                prel_zH[index]->Fill(prel.Z());

                                                pcom_xH[index]->Fill(pcom.X());
                                                pcom_yH[index]->Fill(pcom.Y());
                                                pcom_zH[index]->Fill(pcom.Z());
                                                pcomH[index]->Fill(pcom.Mag());
                                                pcom_qH[index]->Fill(pcom.Dot(vq));
                                                pcom_qtH[index]->Fill(pcom.Dot(vq.Orthogonal()));



                                                cth_pcomH[index]->Fill(pcom.CosTheta());
                                                cth_prelH[index]->Fill(prel.CosTheta());



                                                phi_pcomH[index]->Fill(phi_pcom);
                                                phi_prelH[index]->Fill(phi_prel);


                                                wH[index]->Fill(w);
                                                qH[index]->Fill(vq.Mag());
                                                q2H[index]->Fill(QSq);

                                                bk_xH[index]->Fill(bk_x);

                                                ELepH[index]->Fill(Ek);
                                                thlH[index]->Fill(theta_l);

                                                p4LeadH[index]->Fill(leading_p4->Mag());
                                                p4RecH[index]->Fill(recoil_p4->Mag());
                                                pLeadH[index]->Fill(leading_p4->Vect().Mag());
                                                pLeadzH[index]->Fill(leading_p4->Vect().Z());

                                                pRecH[index]->Fill(recoil_p3.Mag());
                                                ELeadH[index]->Fill(leading_p4->E());
                                                ERecH[index]->Fill(recoil_p4->E());

                                                cthrH[index]->Fill(recoil_p3.CosTheta());
                                                cthsH[index]->Fill(pn.CosTheta());

                                                pmiss_vqH[index]->Fill(P_miss.Dot(vq)/vq.Mag());
                                                pMissH[index]->Fill(P_miss.Mag());
                                                EMissH[index]->Fill(E_miss);

                                                EmissMPmissH[index]->Fill(E_miss - P_miss.Dot(vq)/vq.Mag());
                                                Emiss2MPmiss2H[index]->Fill(E_miss*E_miss - (P_miss.Dot(vq)/vq.Mag())*(P_miss.Dot(vq)/vq.Mag()));

                                                th_pmiqH[index]->Fill(vq.Angle(P_miss)*180/TMath::Pi());

                                                cthlvElepH[index]->Fill(Ek,cthl);
                                                Q2vpnH[index]->Fill(pn.Mag(),QSq);


                                                thetavphilepH[index]->Fill(lep_phi,lepton_p3.Theta()*180/TMath::Pi());
                                                thetavplepH[index]->Fill(lepton_p3.Mag(),lepton_p3.Theta()*180/TMath::Pi());

                                                pLepH[index]->Fill(lepton_p3.Mag());

                                                if(pdg::IsProton(initial_hadrons[1].Pdg())){
                                                        thetavphipH[index]->Fill(lead_phi,leading_p3.Theta()*180/TMath::Pi());
                                                        thetavppH[index]->Fill(leading_p3.Mag(),leading_p3.Theta()*180/TMath::Pi());
                                                }

                                                EmissvPmissH[index]->Fill(P_miss.Mag(),E_miss);

                                                pMiss_eppH[index]->Fill(P_miss.Mag());

                                                pmissvxBH[index]->Fill(bk_x,P_miss.Mag());

                                                pn_genieH[index]->Fill(pn_genie);

                                        }
                                }




                        } else
                        if (ni == 1){

                                SingleNucleon++;


                        }


                        mcrec->Clear();
                }

        }

        cout << "N events: " << nevents << endl;
        cout << "epp events: " << events_pass_epp << endl;
        cout << "ep events: " << events_pass_ep << endl;

        std::cout << "Pairs: "  << PairNucleons << std::endl;
        std::cout << "Single: " << SingleNucleon << std::endl;
        std::cout << "PN Pairs: " << PNPairs << std::endl;
        std::cout << "PP Pairs: " << PPPairs << std::endl;
        std::cout << "NN Pairs: " << NNPairs << std::endl;
        std::cout << "NP Pairs: " << NPPairs << std::endl;


     /*   cout << "Pass Physics: " << pass_physics_cuts << std::endl;
        cout << "Pass Acceptance: " << pass_acceptance_cuts << std::endl;
        cout << "Percentage: " << pass_acceptance_cuts*100/pass_physics_cuts << std::endl;
    *//*    for (int i = 0; i < 50*kmiss_range; i++){
                double kmiss_entries = 0;
                double kmiss_src_entries = 0;
                double kmiss_pur = 0;


                for(int j = i; j < 50*kmiss_range; j++){
                        kmiss_entries += kmiss_hist->GetBinContent(j);
                        kmiss_src_entries += kmiss_src_hist->GetBinContent(j);

                }


                kmiss_pur = kmiss_src_entries/kmiss_entries;
                kmiss_eff_hist->SetBinContent(i+1,kmiss_src_entries*100/kmiss_src_hist->GetEntries());
                kmiss_pur_hist->SetBinContent(i+1,100*kmiss_pur);

        }
        */

        f_output->cd();

        for (int i = 0; i < h_list.size(); i++)
        h_list[i]->Write();

        f_output->Write();
        f_output->Close();

}
