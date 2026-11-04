#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include "constants/particleMasses.hh"
#include "constants/nuclearMasses.hh"
#include "constants/pdgCodes.hh"
#include "constants/geantCodes.hh"

// misc constants
const double GeVfm  =0.1973;
const double alpha = 0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;
const double nbGeVSq = cmSqGeVSq*1.E33;
const double Vud = 0.97417;
// Nucleon magnetic moments in nuclear magnetons, CODATA 2018
// (Tiesinga et al., Rev. Mod. Phys. 93, 025010 (2021), doi:10.1103/RevModPhys.93.025010)
const double mu_p=2.79284734463;
const double mu_n=-1.91304273;
const double carbonDeltaECoulombGeV=0.0029;

#endif
