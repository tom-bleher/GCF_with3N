#!/usr/bin/env python3
"""CLAS6 acceptance for simulated events.

Ports Acceptance/clas6acceptance.C (binary fiducial cuts for e, p) and
Acceptance/AccMap.cpp (weighted acceptance map lookup) to vectorised numpy,
plus Gaussian momentum smearing for finite detector resolution.

Parameters and logic are copied verbatim from those C++ files so the
filter is equivalent to what Natalie used to build the reference hists.

Usage:
    from clas_acceptance import ClasAcceptance
    acc = ClasAcceptance('Acceptance/map_eg2_adin.root', seed=42)
    e_ok = acc.accept_electron(p3_e)     # bool array
    p_ok = acc.accept_proton(p3_p)       # bool array
    w_e  = acc.map_weight(p3_e, 'e')     # float array in [0,1]
    w_p  = acc.map_weight(p3_p, 'p')
    p_smeared = acc.smear_mag(p3, frac_sigma=0.01)
"""
import numpy as np
import uproot

# =====================================================================
# Detector resolutions
# =====================================================================
RESO_ELECTRON = 0.0125   # fractional sigma_p / p
RESO_PROTON   = 0.01

# Minimum recoil-proton momentum for the weighted map (see AccMap.cpp,
# min_prec in constants.h — SRC-cut value of 0.35 GeV/c used elsewhere in
# the paper is consistent).
MIN_PREC = 0.30   # GeV/c — safe detection threshold below SRC cut

# =====================================================================
# PROTON fiducial parameters (from clas6acceptance.C, arrays k*PiPlus[6])
# Each row is a CLAS sector (0..5).
# =====================================================================
#  theta_min(mom) = p0 + p1/mom^2 + p2*mom + p3/mom + p4*exp(mom*p5)
_P_THETA = np.array([
    #  p0        p1          p2           p3         p4         p5
    [ 7.00823,  0.207249,   0.169287,    0.1,       0.1,      -0.1      ],
    [ 5.5,      0.1,        0.506354,    0.1,       3.30779,  -0.651811 ],
    [ 7.06596,  0.127764,  -0.0663754,   0.100003,  4.499,    -3.1793   ],
    [ 6.32763,  0.1,        0.221727,    0.1,       5.30981,  -3.3461   ],
    [ 5.5,      0.211012,   0.640963,    0.1,       3.20347,  -1.10808  ],
    [ 5.5,      0.281549,   0.358452,    0.1,       0.776161, -0.462045 ],
])

#  a(mom)   = p0 + p1*exp(p2*(mom-p3))
#  b(mom)   = p0 + p1*mom*exp(p2*(mom-p3)^2)
_P_A_LOW = np.array([
    [25., -12.,  1.64476,  4.4],
    [25., -12.,  1.51915,  4.4],
    [25., -12.,  1.1095,   4.4],
    [25., -12.,  0.977829, 4.4],
    [25., -12.,  0.955366, 4.4],
    [25., -12.,  0.969146, 4.4],
])
_P_B_LOW = np.array([
    [4.,        2.,        -0.978469, 0.5    ],
    [4.,        2.,        -2.,       0.5    ],
    [2.78427,   2.,        -1.73543,  0.5    ],
    [3.58539,   1.38233,   -2.,       0.5    ],
    [3.32277,   0.0410601, -0.953828, 0.5    ],
    [4.,        2.,        -2.,       1.08576],
])
_P_A_HIGH = np.array([
    [25.,      -11.9735,   0.803484, 4.40024],
    [24.8096,  -8.,        0.85143,  4.8    ],
    [24.8758,  -8.,        1.01249,  4.8    ],
    [25.,     -12.,        0.910994, 4.4    ],
    [25.,      -8.52574,   0.682825, 4.79866],
    [25.,      -8.,        0.88846,  4.8    ],
])
_P_B_HIGH = np.array([
    [2.53606,  0.442034,  -2.,        1.02806 ],
    [2.65468,  0.201149,  -0.179631,  1.6     ],
    [3.17084,  1.27519,   -2.,        0.5     ],
    [2.47156,  1.76076,   -1.89436,   1.03961 ],
    [2.42349,  1.25399,   -2.,        0.815707],
    [2.64394,  0.15892,   -2.,        1.31013 ],
])

# gap_function(mom, (a,b,c)) = a - b*exp(-mom/c)
# gap_params[sector] = [(lower_abc, upper_abc), ...]  for each gap in sector
_P_GAPS = [
    np.array([[[ 35., 10.,  0.4], [ 48., 18.,  0.6]],    # sector 0 (2 gaps)
              [[106., 10.,  0.4], [111.,  6.,  0.4]]]),
    np.array([[[ 88., 15.,  0.6], [ 96., 15.,  0.6]],    # sector 1 (2 gaps)
              [[110., 15.,  0.6], [180.,  0.,  1. ]]]),
    np.array([[[ 19., 50.,  0.7], [ 24., 50.,  0.7]],    # sector 2 (4 gaps)
              [[ 35., 10.,  0.6], [ 44., 13.,  0.6]],
              [[ 74., 25.,  0.4], [108., 20.,  0.4]],
              [[115., 15.,  0.6], [180.,  0.,  1. ]]]),
    np.array([[[100., 18.,  0.6], [112., 18.,  0.6]]]),  # sector 3 (1 gap)
    np.array([[[ 34., 60.,  0.4], [ 38., 50.,  0.4]],    # sector 4 (3 gaps)
              [[ 40., 50.,  0.4], [ 44., 40.,  0.4]],
              [[ 92., 15.,  0.6], [106.,  8.,  0.6]]]),
    np.array([[[ 82., 20.,  0.4], [ 91., 20.,  0.6]],    # sector 5 (2 gaps)
              [[114.,  2.,  0.5], [180.,  0.,  1. ]]]),
]

# =====================================================================
# ELECTRON fiducial parameters (from clas6acceptance.C, arrays k*[6])
# =====================================================================
_E_THETA = np.array([
    #  p0         p1          p2          p3         p4         p5
    [ 15.,     -0.425145,  -0.666294,   5.73077,   10.4976,  -1.13254],
    [ 15.,     -1.02217,   -0.616567,   5.51799,   14.0557,  -1.16189],
    [ 15.,     -0.7837,    -0.673602,   8.05224,   15.2178,  -2.08386],
    [ 15.,     -1.47798,   -0.647113,   7.74737,   16.7291,  -1.79939],
    [ 13.,      3.47361,   -0.34459,    8.45226,  -63.4556,  -3.3791 ],
    [ 13.,      3.5714,    -0.398458,   9.54265,  -22.649,   -1.89746],
])
_E_A_LOW = np.array([
    [25.,     -12.,  0.5605,   4.4],
    [25.,     -12.,  0.714261, 4.4],
    [25.,     -12.,  0.616788, 4.4],
    [24.6345, -12.,  0.62982,  4.4],
    [23.4731, -12.,  1.84236,  4.4],
    [24.8599, -12.,  1.00513,  4.4],
])
_E_B_LOW = np.array([
    [2.1945,  1.51417,  -0.354081, 0.5    ],
    [4.,      1.56882,  -2.,       0.5    ],
    [3.3352,  2.,       -2.,       1.01681],
    [2.22769, 2.,       -0.760895, 1.31808],
    [1.63143, 1.90179,  -0.213751, 0.786844],
    [3.19807, 0.173168, -0.1,      1.6    ],
])
_E_A_HIGH = np.array([
    [25.,      -8.,       0.479446, 4.8    ],
    [25.,     -10.3277,   0.380908, 4.79964],
    [25.,     -12.,       0.675835, 4.4    ],
    [25.,     -11.3361,   0.636018, 4.4815 ],
    [23.7067, -12.,       2.92146,  4.4    ],
    [25.,     -11.4641,   0.55553,  4.41327],
])
_E_B_HIGH = np.array([
    [3.57349, 2.,        -2.,       0.5    ],
    [3.02279, 0.966175,  -2.,       0.527823],
    [2.02102, 2.,        -1.70021,  0.68655],
    [3.1948,  0.192701,  -1.27578,  1.6    ],
    [3.0934,  0.821726,  -0.233492, 1.6    ],
    [2.48828, 2.,        -2.,       0.70261],
])


# =====================================================================
# Vectorised fiducial helpers
# =====================================================================
def _sector(phi_deg):
    """Map phi [-180,180] -> sector index [0,5].

    Follows clas6acceptance.C convention: wrap phi < -30 by +360,
    then sector = floor((phi+30)/60).
    """
    phi = np.where(phi_deg < -30., phi_deg + 360., phi_deg)
    sec = np.floor((phi + 30.) / 60.).astype(np.int64)
    return np.clip(sec, 0, 5), phi


def _theta_min(mom, par_row):
    """theta_min(mom) with one parameter row per event."""
    p0, p1, p2, p3, p4, p5 = [par_row[:, i] for i in range(6)]
    return p0 + p1 / mom**2 + p2 * mom + p3 / mom + p4 * np.exp(mom * p5)


def _a_fn(mom, par_row):
    p0, p1, p2, p3 = [par_row[:, i] for i in range(4)]
    return p0 + p1 * np.exp(p2 * (mom - p3))


def _b_fn(mom, par_row):
    p0, p1, p2, p3 = [par_row[:, i] for i in range(4)]
    return p0 + p1 * mom * np.exp(p2 * (mom - p3)**2)


def _delta_phi(theta, a, b, theta_min):
    return a * (1. - 1. / ((theta - theta_min) / b + 1.))


def _fiducial_pass(mom, theta_deg, phi_deg_signed,
                   theta_pars, aLow_pars, bLow_pars, aHigh_pars, bHigh_pars,
                   require_theta_min=True):
    """Return boolean mask of events inside the fiducial wedge.

    require_theta_min: accept_electron in the C++ code DOES NOT apply the
    theta>minTheta cut (it is commented out!). accept_proton_simple DOES.
    We preserve that difference via this flag.
    """
    sec, phi = _sector(phi_deg_signed)
    tp = theta_pars[sec]
    aL = _a_fn(mom, aLow_pars[sec])
    aH = _a_fn(mom, aHigh_pars[sec])
    bL = _b_fn(mom, bLow_pars[sec])
    bH = _b_fn(mom, bHigh_pars[sec])

    th_min = _theta_min(mom, tp)
    dphi_low  = _delta_phi(theta_deg, aL, bL, th_min)
    dphi_high = _delta_phi(theta_deg, aH, bH, th_min)

    phi_central = 60. * sec
    in_wedge = (phi < phi_central + dphi_high) & (phi > phi_central - dphi_low)
    if require_theta_min:
        in_wedge &= (theta_deg >= th_min)
    return in_wedge


def _gap_accept_proton(mom, theta_deg, phi_deg_signed):
    """True if NOT in any sector gap (proton-only)."""
    sec, _ = _sector(phi_deg_signed)
    keep = np.ones_like(mom, dtype=bool)
    for s in range(6):
        in_sec = (sec == s)
        if not in_sec.any():
            continue
        for gap in _P_GAPS[s]:
            lo_a, lo_b, lo_c = gap[0]
            hi_a, hi_b, hi_c = gap[1]
            min_t = lo_a - lo_b * np.exp(-mom[in_sec] / lo_c)
            max_t = hi_a - hi_b * np.exp(-mom[in_sec] / hi_c)
            in_gap = (theta_deg[in_sec] > min_t) & (theta_deg[in_sec] < max_t)
            idx = np.where(in_sec)[0]
            keep[idx[in_gap]] = False
    return keep


# =====================================================================
# 3-momentum -> (mom, theta, phi) helpers
# =====================================================================
def _vec_angles(px, py, pz):
    mom = np.sqrt(px**2 + py**2 + pz**2)
    safe = np.where(mom > 0, mom, 1e-30)
    theta_deg = np.degrees(np.arccos(np.clip(pz / safe, -1, 1)))
    phi_deg = np.degrees(np.arctan2(py, px))
    return mom, theta_deg, phi_deg


# =====================================================================
# Main class
# =====================================================================
class ClasAcceptance:
    """CLAS6 acceptance (binary fiducial cuts + weighted 3D map) with
    reproducible momentum smearing."""

    def __init__(self, mapfile='Acceptance/map_eg2_adin.root', seed=42,
                 electron_theta_min_deg=None, electron_theta_max_deg=None,
                 map_smoothing_radius=1):
        """map_smoothing_radius : minimum (2r+1)^3 box used in the *first*
        map_weight lookup pass.  The eg2 acceptance map has only ~10
        gen-throws/bin, giving severe per-bin acc/gen fluctuations
        (~50% of bins are exactly 1.0, ~33% exactly 0.0).  Smoothing
        from r=1 (3^3=27 bins) brings effective stats to ~270/bin and
        cures the resulting plot spikiness.  Set 0 to recover the
        bit-faithful AccMap.cpp behaviour (1^3 first pass)."""
        self.rng = np.random.default_rng(seed)
        self.e_th_min = electron_theta_min_deg
        self.e_th_max = electron_theta_max_deg
        self.map_smoothing_radius = int(map_smoothing_radius)
        f = uproot.open(mapfile)
        self._maps = {}
        for particle in ('p', 'e'):
            gen = f[f'solid_{particle}_gen']
            acc = f[f'solid_{particle}_acc']
            self._maps[particle] = dict(
                gen=gen.values(),
                acc=acc.values(),
                mom_edges=gen.axis(0).edges(),
                cos_edges=gen.axis(1).edges(),
                phi_edges=gen.axis(2).edges(),
            )

    # ---- fiducial cuts ---------------------------------------------------
    def accept_electron(self, px, py, pz):
        mom, th, ph = _vec_angles(px, py, pz)
        # Faithful to Natalie's C++ accept_electron: the sector-dependent
        # theta_min rejection is COMMENTED OUT there, so it is not applied
        # here either (require_theta_min=False). An optional hard theta
        # window can be requested via electron_theta_min/max_deg (default
        # off -- not part of her function).
        wedge = _fiducial_pass(mom, th, ph,
                               _E_THETA, _E_A_LOW, _E_B_LOW, _E_A_HIGH, _E_B_HIGH,
                               require_theta_min=False)
        if self.e_th_min is not None:
            wedge = wedge & (th >= self.e_th_min)
        if self.e_th_max is not None:
            wedge = wedge & (th <= self.e_th_max)
        return wedge

    def accept_proton(self, px, py, pz):
        mom, th, ph = _vec_angles(px, py, pz)
        wedge = _fiducial_pass(mom, th, ph,
                               _P_THETA, _P_A_LOW, _P_B_LOW, _P_A_HIGH, _P_B_HIGH,
                               require_theta_min=True)
        no_gap = _gap_accept_proton(mom, th, ph)
        return wedge & no_gap

    # ---- weighted acceptance map ---------------------------------------
    def map_weight(self, px, py, pz, particle='p'):
        """Return array of acc/gen weights for each event.

        Faithful port of AccMap::accept_map (Acceptance/AccMap.cpp):
          - Start with a 1x1x1 bin lookup; if gen==0, grow the summed
            region to 3^3, then 5^3, then 7^3 and re-sum.
          - Out-of-range mom/cos bins contribute 0 (like ROOT's
            GetBinContent for invalid indices); phi bins wrap.
          - If gen is still 0 after 7^3, return 1.0.
        """
        m = self._maps[particle]
        gen, acc = m['gen'], m['acc']
        mom_edges, cos_edges, phi_edges = m['mom_edges'], m['cos_edges'], m['phi_edges']
        nMom = gen.shape[0]
        nCos = gen.shape[1]
        nPhi = gen.shape[2]

        mom = np.sqrt(px**2 + py**2 + pz**2)
        cos_t = np.divide(pz, mom, out=np.zeros_like(mom), where=mom > 0)
        phi_deg = np.degrees(np.arctan2(py, px))
        phi_deg = np.where(phi_deg < phi_edges[0], phi_deg + 360., phi_deg)

        dmom = mom_edges[1] - mom_edges[0]
        dcos = cos_edges[1] - cos_edges[0]
        dphi = phi_edges[1] - phi_edges[0]

        # 0-indexed bin lookups (Python array == ROOT in-range bins 1..nBins).
        iMom = np.floor((mom   - mom_edges[0]) / dmom).astype(np.int64)
        iCos = np.floor((cos_t - cos_edges[0]) / dcos).astype(np.int64)
        iPhi = np.floor((phi_deg - phi_edges[0]) / dphi).astype(np.int64)

        # C++ AccMap.cpp: "if (momBin > 70) momBin=70;" sanitises
        # the upper edge for all particles (see constants.h).
        iMom = np.minimum(iMom, 69)  # 0-indexed -> ROOT bin 70

        iPhi_wrapped = iPhi % nPhi   # phi is periodic in CLAS

        def _sum_box(mom_centre, cos_centre, phi_centre, radius):
            """Sum gen/acc over a (2r+1)^3 box per event, faithfully:
            out-of-range mom/cos bins contribute 0; phi wraps modulo nPhi."""
            out_g = np.zeros(mom_centre.shape[0])
            out_a = np.zeros(mom_centre.shape[0])
            for dcm in range(-radius, radius + 1):
                mB_raw = mom_centre + dcm
                mb_ok  = (mB_raw >= 0) & (mB_raw < nMom)
                mB     = np.where(mb_ok, mB_raw, 0)
                for dcc in range(-radius, radius + 1):
                    cB_raw = cos_centre + dcc
                    cb_ok  = (cB_raw >= 0) & (cB_raw < nCos)
                    cB     = np.where(cb_ok, cB_raw, 0)
                    in_range = mb_ok & cb_ok
                    for dcp in range(-radius, radius + 1):
                        pB = (phi_centre + dcp) % nPhi
                        g_val = gen[mB, cB, pB]
                        a_val = acc[mB, cB, pB]
                        out_g += np.where(in_range, g_val, 0.0)
                        out_a += np.where(in_range, a_val, 0.0)
            return out_g, out_a

        # First-pass lookup: by default smooth over a 3^3 box (radius=1)
        # to overcome the eg2 map's ~10-trial-per-bin Poisson noise.
        # Setting map_smoothing_radius=0 reverts to bit-faithful 1^3.
        r0 = self.map_smoothing_radius
        if r0 <= 0:
            iMom_safe = np.clip(iMom, 0, nMom - 1)
            iCos_safe = np.clip(iCos, 0, nCos - 1)
            g = gen[iMom_safe, iCos_safe, iPhi_wrapped].astype(np.float64, copy=True)
            a = acc[iMom_safe, iCos_safe, iPhi_wrapped].astype(np.float64, copy=True)
            centre_in = ((iMom >= 0) & (iMom < nMom) &
                         (iCos >= 0) & (iCos < nCos))
            g[~centre_in] = 0.0
            a[~centre_in] = 0.0
        else:
            g, a = _sum_box(iMom, iCos, iPhi_wrapped, r0)

        def _expand(bad_idx, radius):
            return _sum_box(iMom[bad_idx], iCos[bad_idx],
                            iPhi_wrapped[bad_idx], radius)

        # Fallback expansion if any event still has gen=0 (rare with the
        # default r0=1 smoothing; expand from r0+1 upward, capped at 3
        # to mirror the AccMap.cpp behaviour of 3 expansion passes).
        for radius in range(max(r0 + 1, 1), 4):
            bad = np.where(g == 0)[0]
            if bad.size == 0:
                break
            g_new, a_new = _expand(bad, radius)
            g[bad] = g_new
            a[bad] = a_new

        # Map gen==0 -> weight = 1 (AccMap.cpp fallback).  Otherwise acc/gen.
        w = np.where(g > 0, a / np.maximum(g, 1e-30), 1.0)
        return w

    # ---- momentum smearing ---------------------------------------------
    def smear_mag(self, px, py, pz, frac_sigma):
        """Gaussian-smear the magnitude; direction unchanged (per Natalie).

        SmearedP = Gaus(|p|, frac_sigma * |p|),  p_vec *= SmearedP / |p|
        """
        mom = np.sqrt(px**2 + py**2 + pz**2)
        sm = mom + self.rng.normal(0.0, frac_sigma * mom)
        sm = np.clip(sm, 0.0, None)             # enforce non-negative momentum
        scale = np.divide(sm, mom, out=np.ones_like(mom), where=mom > 0)
        return px * scale, py * scale, pz * scale, sm


# =====================================================================
# High-level pipeline — smear + fiducial cut + multiply map weights
# =====================================================================
def apply_acceptance_pipeline(d, acc, channel='epp',
                              smear_electron=True, smear_proton=True,
                              e_reso=RESO_ELECTRON, p_reso=RESO_PROTON,
                              min_prec=MIN_PREC, ebeam=None,
                              include_electron_map_weight=True):
    """include_electron_map_weight=True (default): the electron acceptance
    map weight multiplies w_ep / w_epp, symmetric with the proton maps.
    Set False to drop it (binary electron fiducial still applies)."""
    """Given an event dict from compute_kinematics(...), smear and apply
    CLAS acceptance.  Adds:
       d['acc_mask_ep']  — events pass (e,e'p)  fiducial e + p-lead cuts
       d['acc_mask_epp'] — also pass p-recoil fiducial + p-recoil > min_prec
       d['acc_w_ep']     — per-event map weight for (e,e'p)   [map_e * map_pLead]
       d['acc_w_epp']    — per-event map weight for (e,e'pp)  [*= map_pRecoil]

    Smearing is ALSO applied to the kinematic arrays the downstream plot
    code reads (pN, pmiss, mmiss, prec, pcm_*, prel_mag, ...) so the
    histograms reflect the detector resolution."""
    n = len(d['weight'])
    elec_px = d['electron'][:, 0]
    elec_py = d['electron'][:, 1]
    elec_pz = d['electron'][:, 2]
    lead_px = d['lead_post'][:, 0]
    lead_py = d['lead_post'][:, 1]
    lead_pz = d['lead_post'][:, 2]
    rec_px = d['recoil_post'][:, 0]
    rec_py = d['recoil_post'][:, 1]
    rec_pz = d['recoil_post'][:, 2]

    if smear_electron:
        elec_px, elec_py, elec_pz, _ = acc.smear_mag(elec_px, elec_py, elec_pz, e_reso)
    if smear_proton:
        lead_px, lead_py, lead_pz, _ = acc.smear_mag(lead_px, lead_py, lead_pz, p_reso)
        rec_px,  rec_py,  rec_pz,  _ = acc.smear_mag(rec_px,  rec_py,  rec_pz,  p_reso)

    # Binary fiducial masks.
    ok_e    = acc.accept_electron(elec_px, elec_py, elec_pz)
    ok_lead = acc.accept_proton  (lead_px, lead_py, lead_pz)
    ok_rec  = acc.accept_proton  (rec_px,  rec_py,  rec_pz)

    # Weighted map per particle.
    w_e    = acc.map_weight(elec_px, elec_py, elec_pz, 'e')
    w_lead = acc.map_weight(lead_px, lead_py, lead_pz, 'p')
    w_rec  = acc.map_weight(rec_px,  rec_py,  rec_pz,  'p')

    # Combined event weights (no [0,1] clip — faithful to AccMap.cpp).
    # Electron map weight is gated by include_electron_map_weight; the
    # binary electron fiducial (ok_e) is always applied below.
    if include_electron_map_weight:
        w_ep  = w_e * w_lead
        w_epp = w_e * w_lead * w_rec
    else:
        w_ep  = w_lead
        w_epp = w_lead * w_rec

    # Recoil detection threshold (AccMap::recoil_accept uses min_prec)
    prec_mag = np.sqrt(rec_px**2 + rec_py**2 + rec_pz**2)

    # Masks (paired AND product is done lazily in plot code; both survive here)
    mask_ep  = ok_e & ok_lead
    mask_epp = mask_ep & ok_rec & (prec_mag > min_prec)

    d['acc_mask_ep']  = mask_ep
    d['acc_mask_epp'] = mask_epp
    d['acc_w_ep']  = w_ep
    d['acc_w_epp'] = w_epp

    # ---- recompute kinematics with smeared momenta (if smearing on) ----
    if smear_electron or smear_proton:
        # Rebuild 4-vectors: electrons are ultra-relativistic (E = |p|);
        # protons are massive so we redo E.
        M_P = 0.93827
        M_N = 0.93957
        M_AVG = 0.5 * (M_P + M_N)

        elec_E = np.sqrt(elec_px**2 + elec_py**2 + elec_pz**2)
        lead_E = np.sqrt(M_AVG**2 + lead_px**2 + lead_py**2 + lead_pz**2)
        rec_E  = np.sqrt(M_AVG**2 + rec_px**2  + rec_py**2  + rec_pz**2)

        # Recompute q = beam - electron_smeared.  Beam along +z; infer the
        # beam energy from the un-smeared branch (Ebeam = E_e + nu, constant
        # across the sample) if not supplied explicitly.
        if ebeam is None:
            ebeam_arr = d['q'][:, 3] + np.sqrt((d['electron'][:, :3]**2).sum(axis=1))
        else:
            ebeam_arr = np.full(n, float(ebeam))
        q_x = -elec_px
        q_y = -elec_py
        q_z = ebeam_arr - elec_pz
        nu  = ebeam_arr - elec_E
        q_vec = np.column_stack([q_x, q_y, q_z, nu])

        lead3 = np.column_stack([lead_px, lead_py, lead_pz])
        rec3  = np.column_stack([rec_px,  rec_py,  rec_pz])
        q3    = q_vec[:, :3]
        qmag  = np.sqrt((q3**2).sum(axis=1))

        pN_new     = np.sqrt(lead_px**2 + lead_py**2 + lead_pz**2)
        prec_new   = prec_mag
        pmiss_vec  = lead3 - q3
        pmiss_new  = np.sqrt((pmiss_vec**2).sum(axis=1))

        cos_pq = (lead3 * q3).sum(axis=1) / (pN_new * qmag + 1e-30)
        theta_pq_deg_new = np.degrees(np.arccos(np.clip(cos_pq, -1, 1)))

        Emiss  = 2.0 * M_AVG + nu - lead_E
        mmiss2 = Emiss**2 - pmiss_new**2
        mmiss_new = np.where(mmiss2 > 0, np.sqrt(np.clip(mmiss2, 0, None)),
                             -np.sqrt(np.clip(-mmiss2, 0, None)))

        pcm_vec   = lead3 + rec3 - q3
        pcm_mag   = np.sqrt((pcm_vec**2).sum(axis=1))
        prel_vec  = (pmiss_vec - rec3) / 2.0
        prel_mag  = np.sqrt((prel_vec**2).sum(axis=1))

        # CM-coordinate system: z = pmiss_hat, x in pmiss-q plane.
        pmh = pmiss_vec / (pmiss_new[:, None] + 1e-30)
        qz  = (q3 * pmh).sum(axis=1)
        qp  = q3 - qz[:, None] * pmh
        qpm = np.sqrt((qp**2).sum(axis=1))
        xh  = qp / (qpm[:, None] + 1e-30)
        yh  = np.cross(pmh, xh)

        pcm_x = (pcm_vec * xh).sum(axis=1)
        pcm_y = (pcm_vec * yh).sum(axis=1)
        pcm_z = (pcm_vec * pmh).sum(axis=1)

        Q2_new = qmag**2 - nu**2
        xB_new = np.divide(Q2_new, 2.0 * M_AVG * nu,
                           out=np.zeros_like(nu), where=nu > 0)

        d['pN']           = pN_new
        d['pmiss']        = pmiss_new
        d['pmiss_vec']    = pmiss_vec
        d['theta_pq_deg'] = theta_pq_deg_new
        d['mmiss']        = mmiss_new
        d['prec']         = prec_new
        d['pN_over_q']    = pN_new / (qmag + 1e-30)
        d['pcm_vec']      = pcm_vec
        d['pcm_mag']      = pcm_mag
        d['pcm_x']        = pcm_x
        d['pcm_y']        = pcm_y
        d['pcm_z']        = pcm_z
        d['prel_mag']     = prel_mag
        d['q3']           = qmag
        d['nu']           = nu
        d['Q2_calc']      = Q2_new
        d['xB']           = xB_new

    return d


def apply_acceptance_best_recoil(d, acc, smear_proton=True,
                                 p_reso=RESO_PROTON, min_prec=MIN_PREC):
    """(e,e'pp) CLAS acceptance evaluated on the BEST second proton — the
    gen-level recoil OR the highest-momentum FSI-secondary proton, as chosen
    by compute_best_recoil_kinematics (field d['recoil_best']).

    This is required for samples whose detected second proton is FSI-produced
    rather than a correlated recoil — most importantly the single-nucleon
    mean-field sample, which emits NO gen-level recoil (recoil_post = 0) and
    would otherwise be entirely rejected by the gen-recoil acc_mask_epp. For
    SRC samples the best proton is the SRC partner, so the result matches the
    gen-recoil acceptance.

    Requires apply_acceptance_pipeline to have run first (uses acc_mask_ep /
    acc_w_ep for the electron+lead part). Adds:
       d['acc_mask_epp_best'] — e + p-lead + best-p-recoil fiducial, prec>min
       d['acc_w_epp_best']    — acc_w_ep * map_weight(best recoil)
    """
    rb = d['recoil_best']
    rpx, rpy, rpz = rb[:, 0].copy(), rb[:, 1].copy(), rb[:, 2].copy()
    if smear_proton:
        rpx, rpy, rpz, _ = acc.smear_mag(rpx, rpy, rpz, p_reso)
    ok_rec = acc.accept_proton(rpx, rpy, rpz)
    w_rec  = acc.map_weight(rpx, rpy, rpz, 'p')
    prec_mag = np.sqrt(rpx**2 + rpy**2 + rpz**2)
    has = d.get('has_recoil_best', np.ones(len(d['weight']), dtype=bool))

    d['acc_mask_epp_best'] = d['acc_mask_ep'] & has & ok_rec & (prec_mag > min_prec)
    d['acc_w_epp_best']    = d['acc_w_ep'] * w_rec
    return d
