#ifndef __NUCLEAR_MASSES_H__
#define __NUCLEAR_MASSES_H__

#include "particleMasses.hh"

const double mU = 0.9314941024;
const double m_1H = mN;
const double m_2H = 2.01410178 * mU - me;
const double m_3H = 3.01604928199 * mU - me;
const double m_3He = 3.0160293 * mU - 2*me;
const double m_4He = 4.00260325415 * mU - 2*me;
const double m_6Li = 6.015122795 * mU - 3*me;
const double m_8Be = 8.00530510 * mU - 4*me;
const double m_9Be = 9.0121831 * mU - 4*me;
const double m_10Be = 10.013534 * mU - 4*me;
const double m_10B = 10.0129370 * mU - 5*me;
const double m_11B = 11.0093054 * mU - 5*me;
const double m_11C = 11.0114336 * mU - 6*me;
const double m_10C = 10.016853 * mU - 6*me;
const double m_12C = 12. * mU - 6*me;
const double m_14C = 14.00324198862 * mU - 6*me; // AME2020
const double m_14N = 14.0030740048 * mU - 7*me;
const double m_14O = 14.008596706 * mU - 8*me;   // AME2020
const double m_15N = 15.00010889827 * mU - 7*me; // AME2020
const double m_15O = 15.003065636 * mU - 8*me;   // AME2020
const double m_16O = 15.99491461956 * mU - 8*me;

const double m_25Na = 24.9899540 * mU - 11*me;
const double m_25Mg = 24.98583696 * mU - 12*me;
const double m_25Al = 24.99042831 * mU - 13*me;
const double m_27Al = 26.98153841 * mU - 13*me;

const double m_38S = 37.971163 * mU - 16*me;
const double m_38Cl = 37.96801042 * mU - 17*me;
const double m_38Ar = 37.96273210 * mU - 18*me;
const double m_38K = 37.96908112 * mU - 19*me;
const double m_38Ca = 37.97631923 * mU - 20*me;
const double m_40Ar = 39.9623831238 * mU - 18*me;
const double m_40Ca = 39.962590866 * mU - 20*me;

const double m_54Cr = 53.9388804 * mU - 24*me;
const double m_54Mn = 53.9403589 * mU - 25*me;
const double m_54Fe = 53.9396090 * mU - 26*me;
const double m_56Fe = 55.9349363 * mU - 26*me;

const double m_206Hg = 205.977514 * mU - 80*me + 0.5682;
const double m_206Tl = 205.9761103 * mU - 81*me + 0.5682;
const double m_206Pb = 205.9744653 * mU - 82*me + 0.5682;
const double m_208Pb = 207.9766521 * mU - 82*me + 0.5682;

#endif
