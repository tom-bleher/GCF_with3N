#!/usr/bin/env python3
"""
natalie_paper_plots.py — Reproduce plots from Wright et al., arXiv:2104.05090
(doi:10.48550/arXiv.2104.05090)
"Transport Estimations of Final State Interaction Effects on Short-range
 Correlation Studies Using the (e,e'p) and (e,e'pp) Reactions"

Reads PWIA and Full-FSI ROOT TTrees produced by SRC_analysis_2N,
applies kinematic cuts, computes T+SCX reweighting, and generates
all figures from the paper. CLAS experimental data is overlaid on
right (SRC-cut) panels.

Data normalization (paper Fig 2 caption, page 4):
    "The right panel y-axis scale correspond to the measured data counts
     and the calculations are individually area normalized to the data.
     The left panel y-axis scale is arbitrary."

In practice: the CLAS (e,e'pp) data file `SRC_e2p_C_GoodRuns_coulomb.root`
already has the SRC selection cuts applied. For SRC-cut panels (right
column), each PWIA / T+SCX / Full curve is scaled so its integral matches
the data integral over the same bin range. Data points are shown as black
filled circles with sqrt(N) Poisson errors. Left panels (no SRC cuts)
have no data and the curves are area-normalized to the first non-zero
calculation curve.

Sample definitions:
    * Figs 2-7 (PWIA / T+SCX / FULL) use the GCF SRC sample only: AV18
      spectral function with relative momentum above the SRC onset
      (generator wf_mode=2; the onset is either a constant kF or the local
      kF(r), see the generator's src_onset switch). These are the
      `--pwia` / `--fsi` files.
    * Figs 8-10 compare FULL GCF (the same SRC sample) to FULL MF, where the
      mean-field sample is a SINGLE struck nucleon drawn from a Fermi gas
      (generator wf_mode=4): a flat momentum distribution n(p)=3/(4 pi kF^3)
      up to kF (global; local Fermi gas selectable in the generator), with no
      correlated recoil and no high-momentum tail, transported through GENIE
      hA FSI in the A-1 residual. This is a faithful stand-in for the paper's
      GENIE Fermi-gas MF. Any (e,e'pp) second proton comes only from FSI
      secondaries. This is the `--fsi-mf` file.

Event building is detector-style and operates on the generator's flat
final-state particle list (part_* branches): lead = highest-momentum
accepted proton, recoil = highest-momentum other accepted proton, and no
event is ever discarded for having extra protons. The paper's theta_pq,
|p_N|/|q| and |p_recoil| requirements are SRC selection CUTS, not part of
that definition -- see the long note above flatlist_resolve.

Usage:
    python natalie_paper_plots.py [--pwia FILE] [--fsi FILE] [--outdir DIR]
"""

import argparse
import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

# clas_acceptance lives beside this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clas_acceptance import ClasAcceptance, RESO_ELECTRON, RESO_PROTON

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
M_P = 0.93827  # proton mass [GeV]
M_N = 0.93957  # neutron mass [GeV]
M_AVG = 0.5 * (M_P + M_N)  # average nucleon mass
PROTON = 2212
NEUTRON = 2112

# ---------------------------------------------------------------------------
# T+SCX parameters — Eq. 3 of the paper
#
# ALL VALUES ARE EXACT, FOR 12C, FROM THE PRIMARY REFERENCES:
#
# Transparencies: Colle, Cosyn & Ryckebusch, PRC 93, 034608 (2016),
#   doi:10.1103/PhysRevC.93.034608,
#   Table I (page 13 of the preprint), KinB kinematics
#   (CLAS-like: x_B>=1.2, theta_{p1,q}<=25deg, 0.62<|p1|/|q|<0.96)
#
# SCX probabilities: Duer et al., PRL 122, 172502 (2019),
#   doi:10.1103/PhysRevLett.122.172502,
#   Supplemental Material Table III (page 9 of the preprint)
#   computed via RMSGA framework of Colle et al.
# ---------------------------------------------------------------------------
T_SCX = {
    # Single-nucleon transparency T_N for 12C, from Natalie:
    "T_N":  0.53,   # +/- 0.05
    # Two-nucleon transparency T_NN for 12C, from Natalie:
    "T_NN": 0.44,   # +/- 0.04
    # Channel-specific SCX probabilities (Natalie's do_SXC code).
    # Each entry is P(final state given initial state) for a single SCX event;
    # the probability of "no flip" is 1 - sum of the three flipping probs for
    # that initial state.
    "PP2NP": 0.041, "PP2PN": 0.048, "PP2NN": 0.0029,
    "PN2NN": 0.035, "PN2PP": 0.041, "PN2NP": 0.0021,
    "NP2PP": 0.035, "NP2NN": 0.041, "NP2PN": 0.0021,
    "NN2PN": 0.041, "NN2NP": 0.048, "NN2PP": 0.0029,
}

# ---------------------------------------------------------------------------
# SRC selection cuts — page 4 of the paper
# ---------------------------------------------------------------------------
SRC_CUTS = dict(
    xB_min=1.2,
    theta_pq_max_deg=25.0,
    pN_over_q_min=0.62,
    pN_over_q_max=0.92,
    pmiss_min=0.4,   # GeV/c
    pmiss_max=1.0,   # GeV/c
    mmiss_max=1.1,   # GeV
    precoil_min=0.35, # GeV/c — only for (e,e'pp) channel
)

# The same cuts as N. Wright's code actually applies them
# (genie_gcf_loop.C:512-524, pass_selection_cuts). Note the |p_N|/|q| upper
# bound is 0.96 there, while the paper's text says 0.92; the code is what made
# the figures. Everything else is identical to SRC_CUTS.
SRC_CUTS_WRIGHT = dict(SRC_CUTS, pN_over_q_max=0.96)

# Rescatter fate codes, from the legend in genie_gcf_loop.C:1541-1552:
#   0 Undefined  1 No Interaction  2 SCX  3 Inelastic  4 Absorption
#   5 Knock out  6 Compound nucleus  7 Pion Production ...
# Verified independently: on a 2M hA sample, fate 1 was the only value for
# which the final-state particle matched its input nucleon's momentum.
FATE_NO_INTERACTION = 1
FATE_ABSORPTION = 4


# ===================================================================
#  Data loading
# ===================================================================

def load_tree(filepath, treename="events", max_events=None):
    """Load TTree into dict of numpy arrays.

    The per-event final state is the flat particle list written by
    SRC_analysis_2N: part_pdg / part_px / part_py / part_pz / part_E plus
    the provenance columns part_parentRole and part_rescatterCode. Those
    are variable length, so they are kept as awkward arrays; everything
    else becomes a numpy column.

    Note the list holds the COMPLETE final state -- GENIE's cascade output
    verbatim, or the generated nucleons when FSI is off -- not just the
    FSI secondaries. There are no privileged lead/recoil slots to add on
    top of it (the tree's own `lead`/`recoil` branches are built with the
    3N paper's rule and are deliberately ignored here)."""
    f = uproot.open(filepath)
    tree = f[treename]
    all_keys = list(tree.keys())
    if "part_pdg" not in all_keys:
        raise KeyError(
            f"{filepath}: no part_* branches. This file predates the "
            "flat-list generator output; regenerate it with the current "
            "SRC_analysis_2N.")
    # Read only what the analysis uses. The rest -- lead/recoil (built here
    # instead), part_E (rebuilt from |p| and the proton mass), and the tree's
    # xB/nu/pmiss/scattering_angle (all recomputed by compute_kinematics from
    # the SMEARED electron) -- would otherwise be loaded and thrown away. At
    # 150M events that is tens of GB of RAM.
    WANT_FLAT = ("weight", "electron", "q", "gen_lead_type", "gen_rec_type",
                 "lead_pre", "recoil_pre", "Q2")
    WANT_JAG = ("part_pdg", "part_px", "part_py", "part_pz",
                "part_parentRole", "part_rescatterCode", "part_generation")
    jag_keys = [k for k in all_keys if k in WANT_JAG]
    flat_keys = [k for k in all_keys if k in WANT_FLAT]
    missing = [k for k in WANT_FLAT if k not in all_keys]
    if missing:
        raise KeyError(f"{filepath}: missing required branches {missing}")
    stop = None if not max_events else min(int(max_events), tree.num_entries)
    if stop is not None and stop < tree.num_entries:
        print(f"    [cap] reading {stop} of {tree.num_entries} entries")
    out = tree.arrays(flat_keys, library="np", entry_stop=stop)
    jag = tree.arrays(jag_keys, library="ak", entry_stop=stop)
    for k in jag_keys:
        out[k] = jag[k]
    return out


def load_meta(root_path):
    """Load companion <root_path>.meta.txt written by SRC_analysis_2N.
    Returns a dict (strings); returns None if the file is absent."""
    meta_path = root_path + ".meta.txt"
    if not os.path.exists(meta_path):
        return None
    meta = {}
    with open(meta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                meta[k.strip()] = v.strip()
    return meta


def apply_weight_cap(d, percentile, label=""):
    """Cap per-event weights at the given percentile (variance reduction).

    Any event whose weight exceeds the threshold has its weight set to the
    threshold. This is a biased variance-reduction technique: the few
    highest-weight events are down-weighted to reduce histogram noise, at
    the cost of a small bias in shape and in Σw.

    Only positive weights are used to compute the percentile (weight=0
    events — typically FSI-killed ones — are left untouched).
    """
    w = d["weight"]
    pos = w > 0
    if not pos.any():
        print(f"  weight cap on {label}: no positive weights, skipped")
        return
    threshold = float(np.percentile(w[pos], percentile))
    n_capped = int((w > threshold).sum())
    max_before = float(w.max())
    sum_before = float(w.sum())
    d["weight"] = np.where(w > threshold, threshold, w)
    max_after = float(d["weight"].max())
    sum_after = float(d["weight"].sum())
    delta_pct = 100.0 * (sum_after - sum_before) / sum_before if sum_before > 0 else 0.0
    print(f"  cap @ P{percentile:.2f} on {label}: "
          f"w_cut={threshold:.3g}, capped {n_capped}/{len(w)} "
          f"({100.0*n_capped/len(w):.3f}%), "
          f"w_max {max_before:.3g}\u2192{max_after:.3g}, "
          f"\u03a3w {sum_before:.4g}\u2192{sum_after:.4g} ({delta_pct:+.2f}%)")


# ===================================================================
#  Kinematic variable computation
# ===================================================================

def compute_kinematics(d, Ebeam=6.0):
    """Add derived kinematic columns to the dict *d* (in-place).

    Expects branches: electron[4], lead_post[4], recoil_post[4], q[4],
    lead_type, rec_type, weight.
    """
    # Beam 4-vector (along z)
    beam = np.zeros((len(d["weight"]), 4))
    beam[:, 2] = Ebeam
    beam[:, 3] = Ebeam

    # Shortcuts — 3-momentum magnitudes
    def mag3(arr4):
        return np.sqrt(arr4[:, 0]**2 + arr4[:, 1]**2 + arr4[:, 2]**2)

    def dot3(a4, b4):
        return a4[:, 0]*b4[:, 0] + a4[:, 1]*b4[:, 1] + a4[:, 2]*b4[:, 2]

    # Build 4-vectors  (px, py, pz, E) stored column-major
    elec = np.column_stack([d["electron"][:, i] for i in range(4)])
    lead = np.column_stack([d["lead_post"][:, i] for i in range(4)])
    rec  = np.column_stack([d["recoil_post"][:, i] for i in range(4)])
    qvec = np.column_stack([d["q"][:, i] for i in range(4)])

    # Virtual photon
    nu = qvec[:, 3]                            # energy transfer
    q3 = mag3(qvec)                            # |q|
    Q2 = q3**2 - nu**2                         # Q^2 (should match branch)

    # x_B
    xB = Q2 / (2.0 * M_AVG * nu)

    # Leading proton 3-mom
    pN = mag3(lead)

    # Missing momentum: p_miss = p_lead - q
    pmiss_vec = lead[:, :3] - qvec[:, :3]
    pmiss = np.sqrt(pmiss_vec[:, 0]**2 + pmiss_vec[:, 1]**2 + pmiss_vec[:, 2]**2)

    # theta_pq — angle between leading nucleon and q
    cos_theta_pq = dot3(lead, qvec) / (pN * q3 + 1e-30)
    cos_theta_pq = np.clip(cos_theta_pq, -1, 1)
    theta_pq_deg = np.degrees(np.arccos(cos_theta_pq))

    # Missing mass: m_miss = sqrt( (2m + nu - E_lead)^2 - |p_lead - q|^2 )
    # "missing mass of the recoil nucleon assuming breakup of a stationary
    #  two-nucleon pair"
    Emiss = 2.0 * M_AVG + nu - lead[:, 3]
    mmiss_sq = Emiss**2 - pmiss**2
    mmiss = np.where(mmiss_sq > 0, np.sqrt(mmiss_sq), -np.sqrt(-mmiss_sq))

    # Recoil momentum
    prec = mag3(rec)

    # |p_N| / |q|
    pN_over_q = pN / (q3 + 1e-30)

    # CM momentum: p_cm = p_lead + p_recoil - q  (for e,e'pp)
    pcm_vec = lead[:, :3] + rec[:, :3] - qvec[:, :3]
    pcm_mag = np.sqrt(pcm_vec[:, 0]**2 + pcm_vec[:, 1]**2 + pcm_vec[:, 2]**2)

    # Relative momentum: p_rel = (p_miss - p_recoil) / 2
    prel_vec = (pmiss_vec - rec[:, :3]) / 2.0
    prel_mag = np.sqrt(prel_vec[:, 0]**2 + prel_vec[:, 1]**2 + prel_vec[:, 2]**2)

    # Angle between p_miss and q
    cos_pmiss_q = dot3(
        np.column_stack([pmiss_vec, np.zeros(len(pmiss_vec))]),
        qvec
    ) / (pmiss * q3 + 1e-30)
    cos_pmiss_q = np.clip(cos_pmiss_q, -1, 1)
    theta_pmiss_q_deg = np.degrees(np.arccos(cos_pmiss_q))

    # ---- Coordinate system for p_cm components (Fig 6) ----
    # z-axis = p_miss direction, q lies in x-z plane
    # Following Ref [21]: z along p_miss_hat, q in x-z plane
    pmiss_hat = pmiss_vec / (pmiss[:, None] + 1e-30)
    # z = pmiss_hat
    z_hat = pmiss_hat
    # q_perp = q - (q.z_hat) z_hat
    q_dot_z = np.sum(qvec[:, :3] * z_hat, axis=1)
    q_perp = qvec[:, :3] - q_dot_z[:, None] * z_hat
    q_perp_mag = np.sqrt(np.sum(q_perp**2, axis=1))
    x_hat = q_perp / (q_perp_mag[:, None] + 1e-30)
    # y = z x x
    y_hat = np.cross(z_hat, x_hat)

    pcm_x = np.sum(pcm_vec * x_hat, axis=1)
    pcm_y = np.sum(pcm_vec * y_hat, axis=1)
    pcm_z = np.sum(pcm_vec * z_hat, axis=1)

    # Store everything back
    d["xB"] = xB
    d["Q2_calc"] = Q2
    d["nu"] = nu
    d["q3"] = q3
    d["pN"] = pN
    d["pmiss"] = pmiss
    d["pmiss_vec"] = pmiss_vec
    d["theta_pq_deg"] = theta_pq_deg
    d["mmiss"] = mmiss
    d["prec"] = prec
    d["pN_over_q"] = pN_over_q
    d["pcm_vec"] = pcm_vec
    d["pcm_mag"] = pcm_mag
    d["pcm_x"] = pcm_x
    d["pcm_y"] = pcm_y
    d["pcm_z"] = pcm_z
    d["prel_mag"] = prel_mag
    d["theta_pmiss_q_deg"] = theta_pmiss_q_deg

    return d


# ===================================================================
#  Flat-list final-state selection (detector-style)
# ===================================================================
#
# The generator's lead/recoil slots and parentRole tags are cascade
# bookkeeping: after FSI (e.g. an elastic pp scatter) "which proton was
# the original" is not an observable, so the analysis must not privilege
# any particle in the list. Every final-state nucleon -- the surviving
# struck nucleon, its correlated partner, and anything the cascade
# knocked out -- enters ONE flat candidate list, and the channels are
# built the way a detector builds them:
#
#   lead   = the highest-momentum accepted proton;
#   recoil = the highest-momentum OTHER accepted proton;
#   (e,e'p)  = electron accepted and >= 1 accepted proton;
#   (e,e'pp) = electron accepted and >= 2 accepted protons.
#
# An event is never discarded for having "too many" protons.
#
# The paper gives no tie-break rule -- it only names "the leading
# (knocked-out) proton" and "the detected second proton", and says the
# two samples "differ only in the requirement of detecting a recoil
# proton in the final state". It does settle where the kinematic cuts
# belong, though: theta_pq < 25 deg and 0.62 < |p_N|/|q| < 0.92 are listed
# under *event selection cuts*, and Fig. 2's left panel ("all events that
# pass the CLAS acceptance cuts") carries substantial FULL yield at
# p_lead ~ 0.5 GeV/c. With |q| ~ 1.5-2.5 GeV/c that is |p_N|/|q| ~ 0.2-0.33,
# far below 0.62, so those events cannot exist if the window were part of
# the lead DEFINITION. They live in apply_src_cuts only, together with
# |p_recoil| > 0.35 GeV/c.
#
# The same convention applies to EVERY sample (PWIA, FSI SRC, FSI MF).
# In the MF sample the second proton can only be an FSI product -- which
# is the whole point of that comparison.
#
# For the PWIA sample the candidate list must include NEUTRONS
# (include_neutrons=True): the T+SCX reweighting assigns each event the
# probability that zero-angle SCX turns its (lead, recoil) isospins into
# the measured proton final state, so neutron legs carry (e,e'p[n]) and
# (e,e'[n]p) contributions with the same kinematics. The chosen
# candidates' true types are kept in lead_type / rec_type, and the
# measured-channel selections (select_eep / select_eepp) still demand
# protons.

_KMAX_CAND = 24   # cap on the per-event final-state list

# NOTE: no explicit low-momentum proton threshold is applied. The eg2
# acceptance map already handles it, and does so probabilistically rather
# than with a hard edge -- measured mean map weight for an isotropic proton
# throw: 0.000 at 0.05 GeV/c, 0.005 at 0.10, 0.045 at 0.15, 0.42 at 0.20,
# 0.63 at 0.25, 0.80 at 0.60. The gen==0 -> weight 1.0 fallback inside
# map_weight never fires down there (0.0% of accepted protons below
# 0.30 GeV/c), so there is no soft-proton leak to guard against. A detected
# proton is a detected proton, whatever produced it; the paper's
# |p_recoil| > 0.35 GeV/c requirement is an SRC selection cut and lives in
# apply_src_cuts.


def _free_raw_columns(d):
    """Drop the per-particle arrays once lead/recoil have been resolved.

    The part_* branches exist only to build the event: every quantity the
    figures need is already extracted into fixed-length columns (lead_post,
    recoil_post, lead_role/rec_role, lead_fate/rec_fate). For a 150M-event
    sample the jagged arrays and their awkward overhead are tens of GB, and
    the analysis holds three samples at once, so keeping them is what pushed
    a 60M/60M/150M run to 169 GB RSS and nearly OOMed the node.

    lead_pre / recoil_pre go too: they are used only for the truth k_rel cut,
    which has already been applied by the time this is called.
    """
    import gc
    for k in ("part_pdg", "part_px", "part_py", "part_pz",
              "part_parentRole", "part_rescatterCode", "part_generation",
              "lead_pre", "recoil_pre"):
        d.pop(k, None)
    gc.collect()
    return d


def wright_resolve(d, acc, ebeam, seed=None, chunk=2_000_000, phi_rot=True):
    """Event building following N. Wright's analysis macro exactly.

    This is a SECOND, selectable convention -- flatlist_resolve is untouched
    and remains the default. The point is to test whether her conventions
    reproduce her published numbers, not to assert that they are the right
    ones.

    Transcribed from genie_gcf_loop.C (the macro whose output filename matches
    the file gst_plotter.C reads for the FULL curve):

      * 1503-1521  the final-state list holds PROTONS ONLY. Each is smeared
                   (1% on the magnitude), then must pass the proton fiducial
                   AND have a NON-ZERO acceptance-map value -- the map is used
                   as a binary filter here, not a weight.
      * 1554-1557  "Goodbye Absorption": if EITHER SRC nucleon has rescatter
                   fate 4 (absorption), the event is discarded.
      * 1592-1648  lead = the highest-momentum accepted proton. Her first loop
                   intends to apply the SRC cuts to lead candidates, but 1618
                   compares against `pLead_max` -- a histogram axis limit of
                   2.4 GeV -- instead of `lead_max`, so it almost never fires
                   and the fallback loop (1628-1640) decides. We reproduce the
                   effective behaviour, which is simply highest-|p|.
      * 1652-1740  recoil candidates = all other protons above 0.35 GeV/c,
                   falling back to all other protons if none qualify.
                   Provenance is never enforced: a nucleon the cascade knocked
                   out of the medium is as valid a recoil as the SRC partner.
                   1 candidate -> take it. 2 or 3 -> treat each candidate's
                   map value as an independent detection probability and throw
                   ONE random number; if two or more come out detected,
                   DISCARD THE WHOLE EVENT (including its (e,e'p) entry).
                   4+ candidates -> her dispatch is commented out, so the
                   event silently becomes (e,e'p)-only.
      * 1691       when the lottery was used, the recoil map weight is NOT
                   applied afterwards -- it has already been consumed.
      * 1824-1828  an analysis-level truth cut, p_rel > 0.25 GeV/c.
      * 1425, 1514 the electron and every proton are rotated by +pi before
                   the fiducials and the map gate (phi_rot). The lead weight
                   (1853) and the lottery (1697-1709) read the map at the
                   UNSMEARED, UNROTATED momentum instead -- reproduced as is.

    Sets the same columns as flatlist_resolve, plus d["wright_drop"] marking
    events her code discards.
    """
    import awkward as ak
    rng = np.random.default_rng(seed)
    n = len(d["weight"])
    K = _KMAX_CAND

    d["gen_lead_type"] = np.asarray(d["gen_lead_type"]).copy()
    d["gen_rec_type"]  = np.asarray(d["gen_rec_type"]).copy()

    def _pad(field, fill):
        return ak.to_numpy(ak.fill_none(ak.pad_none(field, K, clip=True), fill))

    out = dict(
        e3=np.zeros((n, 3)), ok_e=np.zeros(n, bool),
        lead3=np.zeros((n, 3)), rec3=np.zeros((n, 3)),
        ltype=np.zeros(n, np.int64), rtype=np.zeros(n, np.int64),
        lrole=np.full(n, -1, np.int64), rrole=np.full(n, -1, np.int64),
        lfate=np.full(n, -1, np.int64), rfate=np.full(n, -1, np.int64),
        wlead=np.ones(n), wrec=np.ones(n), ncand=np.zeros(n, np.int64),
        drop=np.zeros(n, bool),
    )

    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n); m = hi - lo
        rows = np.arange(m)

        pdg   = _pad(d["part_pdg"][lo:hi], 0).astype(np.int64)
        cand  = np.stack([_pad(d["part_px"][lo:hi], 0.0),
                          _pad(d["part_py"][lo:hi], 0.0),
                          _pad(d["part_pz"][lo:hi], 0.0)], axis=2)
        crole = _pad(d["part_parentRole"][lo:hi], -1).astype(np.int64)
        cfate = _pad(d["part_rescatterCode"][lo:hi], -1).astype(np.int64)
        cand_true = cand.copy()          # unsmeared, as her final_state_p4

        # --- "Goodbye Absorption" (1554-1557) ---------------------------
        # The fate is a property of the INPUT nucleon, so it is constant
        # across every final-state particle sharing a parentRole.
        absorbed = np.zeros(m, bool)
        for role in (0, 1):
            has = (crole == role)
            fate_role = np.where(has.any(axis=1),
                                 np.max(np.where(has, cfate, -1), axis=1), -1)
            absorbed |= (fate_role == FATE_ABSORPTION)
        out["drop"][lo:hi] |= absorbed

        # --- proton list: smear, fiducial, non-zero map (1503-1521) ------
        e3 = d["electron"][lo:hi, :3].copy()
        is_p = (pdg == PROTON)
        if acc is not None:
            ex, ey, ez, _ = acc.smear_mag(e3[:, 0], e3[:, 1], e3[:, 2],
                                          RESO_ELECTRON)
            e3 = np.column_stack([ex, ey, ez])
            fx, fy, fz, _ = acc.smear_mag(cand[:, :, 0].ravel(),
                                          cand[:, :, 1].ravel(),
                                          cand[:, :, 2].ravel(), RESO_PROTON)
            cand = np.stack([fx.reshape(m, K), fy.reshape(m, K),
                             fz.reshape(m, K)], axis=2)
            with np.errstate(divide="ignore", invalid="ignore"):
                fid = acc.accept_proton(cand[:, :, 0].ravel(),
                                        cand[:, :, 1].ravel(),
                                        cand[:, :, 2].ravel(),
                                        rot=phi_rot).reshape(m, K)
                wgate = acc.map_weight(cand[:, :, 0].ravel(),
                                       cand[:, :, 1].ravel(),
                                       cand[:, :, 2].ravel(),
                                       'p', rot=phi_rot).reshape(m, K)
                # 1697-1709: the lottery reads the map at the unsmeared,
                # unrotated momentum (final_state_p4), not at the gate's.
                wmap = acc.map_weight(cand_true[:, :, 0].ravel(),
                                      cand_true[:, :, 1].ravel(),
                                      cand_true[:, :, 2].ravel(),
                                      'p').reshape(m, K)
            ok = is_p & fid & (wgate > 0)         # map used as a BINARY filter
            ok_e = acc.accept_electron(e3[:, 0], e3[:, 1], e3[:, 2],
                                       rot=phi_rot)
        else:
            wmap = np.ones((m, K)); ok = is_p
            ok_e = np.ones(m, bool)

        pmag = np.sqrt((cand ** 2).sum(axis=2))
        p_ok = np.where(ok, pmag, -1.0)

        # --- lead: highest-|p| accepted proton (effective 1628-1640) -----
        lead_idx = np.argmax(p_ok, axis=1)
        has_lead = p_ok[rows, lead_idx] > 0.0
        out["drop"][lo:hi] |= ~has_lead            # her no_lead++ / continue

        # --- recoil candidates (1652-1687) ------------------------------
        others = ok.copy(); others[rows, lead_idx] = False
        above = others & (pmag > SRC_CUTS["precoil_min"])
        n_above = above.sum(axis=1)
        pool = np.where(n_above[:, None] > 0, above, others)   # her fallback
        ncand = pool.sum(axis=1)
        out["ncand"][lo:hi] = ncand

        rec_idx = np.full(m, -1, np.int64)
        use_rec_acc = np.zeros(m, bool)

        one = ncand == 1
        if one.any():
            rec_idx[one] = np.argmax(pool[one], axis=1)
            use_rec_acc[one] = True                # 1691: weight only here

        # 2 or 3 candidates -> her assign_recoil lottery (540-623, 1695-1731)
        for k in (2, 3):
            sel = ncand == k
            if not sel.any():
                continue
            idx = np.argsort(~pool[sel], axis=1)[:, :k]        # candidate slots
            probs = np.take_along_axis(wmap[sel], idx, axis=1) # her map values
            hit = rng.random(probs.shape) < np.clip(probs, 0, 1)
            nhit = hit.sum(axis=1)
            rowsel = np.where(sel)[0]
            multi = nhit >= 2
            out["drop"][lo + rowsel[multi]] = True             # 1728-1731
            single = nhit == 1
            if single.any():
                pick = np.argmax(hit[single], axis=1)
                rec_idx[rowsel[single]] = idx[single, pick]
            # nhit == 0 -> no recoil, event stays (e,e'p)

        # 4+ candidates: her dispatch is commented out -> no recoil assigned.

        has_rec = rec_idx >= 0
        ri = np.where(has_rec, rec_idx, 0)

        out["e3"][lo:hi] = e3
        out["ok_e"][lo:hi] = ok_e
        out["lead3"][lo:hi] = np.where(has_lead[:, None], cand[rows, lead_idx], 0.0)
        out["rec3"][lo:hi] = np.where(has_rec[:, None], cand[rows, ri], 0.0)
        out["ltype"][lo:hi] = np.where(has_lead, pdg[rows, lead_idx], 0)
        out["rtype"][lo:hi] = np.where(has_rec, pdg[rows, ri], 0)
        out["lrole"][lo:hi] = np.where(has_lead, crole[rows, lead_idx], -1)
        out["rrole"][lo:hi] = np.where(has_rec, crole[rows, ri], -1)
        out["lfate"][lo:hi] = np.where(has_lead, cfate[rows, lead_idx], -1)
        out["rfate"][lo:hi] = np.where(has_rec, cfate[rows, ri], -1)
        # 1853: lead weight from the UNSMEARED momentum
        if acc is not None:
            lt = cand_true[rows, lead_idx]
            out["wlead"][lo:hi] = np.where(
                has_lead, acc.map_weight(lt[:, 0], lt[:, 1], lt[:, 2], 'p'), 1.0)
            rt = cand_true[rows, ri]
            wr = acc.map_weight(rt[:, 0], rt[:, 1], rt[:, 2], 'p')
            out["wrec"][lo:hi] = np.where(has_rec & use_rec_acc, wr, 1.0)

    has_lead = out["ltype"] != 0
    has_rec = out["rtype"] != 0
    d["lead_type"] = out["ltype"].astype(np.int32)
    d["rec_type"] = out["rtype"].astype(np.int32)
    d["lead_post"] = np.column_stack(
        [out["lead3"], np.sqrt((out["lead3"] ** 2).sum(axis=1) + M_P ** 2)])
    d["recoil_post"] = np.column_stack(
        [out["rec3"], np.sqrt((out["rec3"] ** 2).sum(axis=1) + M_P ** 2)])
    d["lead_role"], d["rec_role"] = out["lrole"], out["rrole"]
    d["lead_fate"], d["rec_fate"] = out["lfate"], out["rfate"]
    d["lead_gen"] = np.full(len(out["lrole"]), -1)
    d["rec_gen"] = np.full(len(out["rrole"]), -1)

    Ee = np.sqrt((out["e3"] ** 2).sum(axis=1))
    d["electron"] = np.column_stack([out["e3"], Ee])
    d["q"] = np.column_stack([-out["e3"][:, 0], -out["e3"][:, 1],
                              ebeam - out["e3"][:, 2], ebeam - Ee])
    compute_kinematics(d, Ebeam=ebeam)

    # 1824-1828: analysis-level TRUTH cut, prel > 0.25 GeV/c.
    #
    #     TVector3 prel = .5*(pn - true_recoil->P4()->Vect()); //truth IS
    #     if(prel.Mag() < .25){ ... continue; }
    #
    # At that point `pn` is NOT the reconstructed missing momentum -- it is the
    # struck nucleon's initial momentum, set at 1476-1479 from the status-11
    # particle, and only reassigned to the smeared p_miss on the line AFTER the
    # cut (1832). `true_recoil` is the first kIStHadronInTheNucleus nucleon,
    # which for an SRC pair is the partner (SRCNuclearRecoil adds it at
    # Module-2, before the struck nucleon's daughter at Module-5). So for the
    # SRC samples this is exactly the generated relative momentum,
    # 0.5|p1_gen - p2_gen| -- i.e. her analysis imposes k_rel > 0.25 on top of
    # whatever her generator produced.
    #
    # For the mean-field sample there is no partner: `true_recoil` is then the
    # OUTGOING nucleon (p1 + q), so her prel is 0.5|q| ~ 0.75-1.2 GeV/c and the
    # cut never fires. It must therefore not be applied to MF.
    p1_gen = np.asarray(d["lead_pre"])[:, :3] - np.asarray(d["q"])[:, :3]
    p2_gen = np.asarray(d["recoil_pre"])[:, :3]
    prel_truth = 0.5 * np.linalg.norm(p1_gen - p2_gen, axis=1)
    has_partner = np.isin(np.asarray(d["gen_rec_type"]), (PROTON, NEUTRON))
    d["wright_prel_truth"] = np.where(has_partner, prel_truth, np.inf)
    out["drop"] |= has_partner & (prel_truth < 0.25)

    keep = ~out["drop"]
    d["wright_drop"] = out["drop"]
    d["acc_mask_ep"] = out["ok_e"] & has_lead & keep
    d["acc_mask_epp"] = out["ok_e"] & has_rec & keep
    d["acc_w_ep"] = out["wlead"]
    d["acc_w_epp"] = out["wlead"] * out["wrec"]
    d["n_lead_cand"] = out["ncand"] + has_lead.astype(np.int64)
    d["n_recoil_cand"] = out["ncand"]
    _free_raw_columns(d)
    return d


def flatlist_resolve(d, acc, ebeam, include_neutrons=False, chunk=2_000_000,
                     phi_rot=True, absorption_veto=True, epp_recoil="unique"):
    """Detector-style event building on the flat final-state particle list.

    Smears the electron and every candidate (when `acc` is given), applies
    the CLAS fiducials, then picks lead = highest-|p| accepted proton and
    recoil = highest-|p| other accepted proton (see the note above).
    Overwrites d's electron / q / lead_post / lead_type / recoil_post /
    rec_type plus the acceptance masks (acc_mask_ep/epp) and map weights
    (acc_w_ep/epp), and recomputes every derived kinematic column.

    Also records, for the chosen lead and recoil, the provenance columns
    lead_role / rec_role (part_parentRole: 0 = descends from the generated
    struck nucleon, 1 = from the generated partner) and lead_gen / rec_gen
    (part_generation: 0 = IS that nucleon, untouched by the cascade; >0 = an
    FSI product). These drive the candidate-origin diagnostic.

    Diagnostics d["n_lead_cand"] / d["n_recoil_cand"] count accepted
    candidates. Processed in chunks of `chunk` events to bound peak memory.

    phi_rot          evaluate every fiducial and map lookup at phi + 180 deg,
                     as N. Wright's macro does for the electron and protons.
    absorption_veto  throw the event away -- weight set to 0 and removed from
                     every channel -- when either generated nucleon has hA
                     fate 4 (absorption). ON by default. The generator-level
                     weight is kept in weight_gen, used ONLY by the checks
                     that compare generator configurations (pre-cut Q^2
                     stage, Sigma-w vs attempts).
    epp_recoil       how the (e,e'pp) recoil is chosen. (e,e'p) never looks
                     at recoils under any setting.
        "unique"    recoil candidates are the accepted protons other than
                    the lead with |p| >= precoil_min. Exactly one -> it is
                    the recoil; none -> no (e,e'pp); two or more ->
                    ambiguous, the event is NOT (e,e'pp) but stays in
                    (e,e'p). Softer protons and any other particles are
                    allowed as extras.
        "exactly2"  (e,e'pp) only when exactly two protons are accepted in
                    total (lead + one other, any momentum). This is the
                    topology of the CLAS data file we compare to, which
                    contains only two-proton events.
        "highest"   the previous rule: recoil = second-highest accepted
                    proton, extras ignored.
        "all"       every candidate with |p| >= precoil_min is a recoil: the
                    highest takes the event's own row and each further one
                    is APPENDED as an extra (e,e'pp)-only row (epp_extra
                    True) with that proton as the recoil and its own map
                    weight. (e,e'p) still counts the event once.

    Also sets n_recoil_hard, the number of possible recoils (accepted protons
    other than the lead with |p| >= precoil_min), for every rule.
    """
    import awkward as ak
    n = len(d["weight"])
    K = _KMAX_CAND
    allowed = (PROTON, NEUTRON) if include_neutrons else (PROTON,)

    # Preserve the GENERATED struck/partner isospins before the
    # detector-style resolution overwrites lead_type / rec_type with the
    # DETECTED final-state types. The T+SCX reweighting is a model
    # correction applied to the generated PWIA pair, so it must use these
    # generated isospins, not the detected recoil (which is absent for the
    # events with no detected second nucleon).
    d["gen_lead_type"] = np.asarray(d["gen_lead_type"]).copy()
    d["gen_rec_type"]  = np.asarray(d["gen_rec_type"]).copy()

    n_over = int(ak.sum(ak.num(d["part_pdg"], axis=1) > K))
    if n_over:
        print(f"    WARNING: {n_over} event(s) have more than {K} final-state "
              f"particles; the list is clipped (raise _KMAX_CAND).")

    out_e3    = np.zeros((n, 3))
    out_ok_e  = np.zeros(n, dtype=bool)
    out_lead3 = np.zeros((n, 3))
    out_rec3  = np.zeros((n, 3))
    out_ltype = np.zeros(n, dtype=np.int64)
    out_rtype = np.zeros(n, dtype=np.int64)
    out_lrole = np.full(n, -1, dtype=np.int64)
    out_rrole = np.full(n, -1, dtype=np.int64)
    out_lgen  = np.full(n, -1, dtype=np.int64)
    out_rgen  = np.full(n, -1, dtype=np.int64)
    out_lfate = np.full(n, -1, dtype=np.int64)
    out_rfate = np.full(n, -1, dtype=np.int64)
    out_nlead = np.zeros(n, dtype=np.int64)
    out_nrec  = np.zeros(n, dtype=np.int64)
    out_abs   = np.zeros(n, dtype=bool)
    out_ambig = np.zeros(n, dtype=bool)
    out_nhard = np.zeros(n, dtype=np.int64)
    xr_l, x3_l, xt_l, xro_l, xg_l, xf_l = [], [], [], [], [], []

    def _pad(field, m, fill):
        return ak.to_numpy(
            ak.fill_none(ak.pad_none(field, K, clip=True), fill))

    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        m  = hi - lo

        # ---- candidate list: the flat final state, verbatim ----
        ctype = _pad(d["part_pdg"][lo:hi], m, 0).astype(np.int64)
        cand  = np.stack([_pad(d["part_px"][lo:hi], m, 0.0),
                          _pad(d["part_py"][lo:hi], m, 0.0),
                          _pad(d["part_pz"][lo:hi], m, 0.0)], axis=2)
        crole = _pad(d["part_parentRole"][lo:hi], m, -1).astype(np.int64)
        cgen  = _pad(d["part_generation"][lo:hi], m, -1).astype(np.int64)
        cfate = _pad(d["part_rescatterCode"][lo:hi], m, -1).astype(np.int64)
        valid = np.isin(ctype, allowed)

        # ---- optional absorption veto (genie_gcf_loop.C:1554-1557) ----
        # The fate is stamped on the INPUT nucleon, so it is constant across
        # every final-state particle sharing a parentRole.
        if absorption_veto:
            absorbed = np.zeros(m, bool)
            for role in (0, 1):
                has = (crole == role)
                fate_role = np.where(has.any(axis=1),
                                     np.max(np.where(has, cfate, -1), axis=1),
                                     -1)
                absorbed |= (fate_role == FATE_ABSORPTION)
            out_abs[lo:hi] = absorbed

        # ---- detector-resolution smearing ----
        e3 = d["electron"][lo:hi, :3].copy()
        if acc is not None:
            ex, ey, ez, _ = acc.smear_mag(e3[:, 0], e3[:, 1], e3[:, 2],
                                          RESO_ELECTRON)
            e3 = np.column_stack([ex, ey, ez])
            fx, fy, fz, _ = acc.smear_mag(cand[:, :, 0].ravel(),
                                          cand[:, :, 1].ravel(),
                                          cand[:, :, 2].ravel(), RESO_PROTON)
            cand = np.stack([fx.reshape(m, K), fy.reshape(m, K),
                             fz.reshape(m, K)], axis=2)

        # ---- "detected" = valid candidate inside the CLAS proton fiducials.
        # Neutron candidates (T+SCX bookkeeping) use the same proton
        # fiducials: zero-angle SCX preserves the flipped leg's kinematics.
        pmag = np.sqrt((cand**2).sum(axis=2))
        detected = valid
        if acc is not None:
            # Padded slots have zero momentum; the fiducial parametrisations
            # divide by |p|, so silence the resulting inf/nan. Those slots are
            # already excluded by `valid`.
            with np.errstate(divide="ignore", invalid="ignore"):
                fid = acc.accept_proton(cand[:, :, 0].ravel(),
                                        cand[:, :, 1].ravel(),
                                        cand[:, :, 2].ravel(),
                                        rot=phi_rot).reshape(m, K)
            detected = detected & fid
            ok_e = acc.accept_electron(e3[:, 0], e3[:, 1], e3[:, 2],
                                       rot=phi_rot)
        else:
            ok_e = np.ones(m, dtype=bool)

        # ---- lead / recoil: the two highest-momentum accepted nucleons ----
        p_ok = np.where(detected, pmag, -1.0)
        rows = np.arange(m)

        lead_idx = np.argmax(p_ok, axis=1)
        has_lead = p_ok[rows, lead_idx] > 0.0

        p_rest = p_ok.copy()
        p_rest[rows, lead_idx] = -1.0
        rec_idx = np.argmax(p_rest, axis=1)       # highest-|p| other nucleon
        n_other = (p_rest > 0.0).sum(axis=1)
        n_hard = (p_rest >= SRC_CUTS["precoil_min"]).sum(axis=1)
        if epp_recoil == "unique":
            has_rec = has_lead & (n_hard == 1)
            ambig = has_lead & (n_hard >= 2)
        elif epp_recoil == "exactly2":
            has_rec = has_lead & (n_other == 1)
            ambig = has_lead & (n_other >= 2)
        elif epp_recoil == "highest":
            has_rec = has_lead & (n_other > 0)
            ambig = np.zeros(m, bool)
        elif epp_recoil == "all":
            # The highest hard candidate (== rec_idx whenever n_hard >= 1)
            # takes this row; every further hard candidate becomes an extra row.
            has_rec = has_lead & (n_hard >= 1)
            ambig = np.zeros(m, bool)
            extra = (p_rest >= SRC_CUTS["precoil_min"]) & has_lead[:, None]
            extra[rows, rec_idx] = False
            er, ec = np.nonzero(extra)
            if er.size:
                xr_l.append(lo + er)
                x3_l.append(cand[er, ec])
                xt_l.append(ctype[er, ec])
                xro_l.append(crole[er, ec])
                xg_l.append(cgen[er, ec])
                xf_l.append(cfate[er, ec])
        else:
            raise ValueError(f"unknown epp_recoil={epp_recoil!r}")
        out_ambig[lo:hi] = ambig
        out_nhard[lo:hi] = np.where(has_lead, n_hard, 0)

        out_e3[lo:hi]    = e3
        out_ok_e[lo:hi]  = ok_e
        out_lead3[lo:hi] = np.where(has_lead[:, None], cand[rows, lead_idx], 0.0)
        out_rec3[lo:hi]  = np.where(has_rec[:, None],  cand[rows, rec_idx],  0.0)
        out_ltype[lo:hi] = np.where(has_lead, ctype[rows, lead_idx], 0)
        out_rtype[lo:hi] = np.where(has_rec,  ctype[rows, rec_idx],  0)
        out_lrole[lo:hi] = np.where(has_lead, crole[rows, lead_idx], -1)
        out_rrole[lo:hi] = np.where(has_rec,  crole[rows, rec_idx],  -1)
        out_lgen[lo:hi]  = np.where(has_lead, cgen[rows, lead_idx], -1)
        out_rgen[lo:hi]  = np.where(has_rec,  cgen[rows, rec_idx],  -1)
        out_lfate[lo:hi] = np.where(has_lead, cfate[rows, lead_idx], -1)
        out_rfate[lo:hi] = np.where(has_rec,  cfate[rows, rec_idx],  -1)
        out_nlead[lo:hi] = detected.sum(axis=1)
        out_nrec[lo:hi]  = np.maximum(detected.sum(axis=1) - 1, 0)

    has_lead = out_ltype != 0
    has_rec  = out_rtype != 0
    d["lead_type"] = out_ltype.astype(np.int32)
    d["rec_type"]  = out_rtype.astype(np.int32)
    d["lead_post"] = np.column_stack(
        [out_lead3, np.sqrt((out_lead3**2).sum(axis=1) + M_P**2)])
    d["recoil_post"] = np.column_stack(
        [out_rec3, np.sqrt((out_rec3**2).sum(axis=1) + M_P**2)])
    d["lead_role"]  = out_lrole
    d["rec_role"]   = out_rrole
    d["lead_gen"] = out_lgen
    d["rec_gen"]  = out_rgen
    d["lead_fate"] = out_lfate
    d["rec_fate"]  = out_rfate

    # ---- overwrite electron / q, recompute all derived kinematics ----
    Ee = np.sqrt((out_e3**2).sum(axis=1))
    d["electron"] = np.column_stack([out_e3, Ee])
    d["q"] = np.column_stack([-out_e3[:, 0], -out_e3[:, 1],
                              ebeam - out_e3[:, 2], ebeam - Ee])
    compute_kinematics(d, Ebeam=ebeam)

    # ---- acceptance masks & map weights on the CHOSEN particles ----
    keep = ~out_abs
    d["acc_mask_ep"]  = out_ok_e & has_lead & keep
    d["acc_mask_epp"] = out_ok_e & has_rec & keep
    d["abs_drop"] = out_abs
    d["epp_ambiguous"] = out_ambig
    d["weight_gen"] = np.asarray(d["weight"]).copy()
    d["weight"] = np.where(out_abs, 0.0, d["weight_gen"])
    w_e    = np.ones(n)
    w_lead = np.ones(n)
    w_rec  = np.ones(n)
    if acc is not None:
        w_e = acc.map_weight(out_e3[:, 0], out_e3[:, 1], out_e3[:, 2], 'e',
                             rot=phi_rot)
        if has_lead.any():
            w_lead[has_lead] = acc.map_weight(
                out_lead3[has_lead, 0], out_lead3[has_lead, 1],
                out_lead3[has_lead, 2], 'p', rot=phi_rot)
        if has_rec.any():
            w_rec[has_rec] = acc.map_weight(
                out_rec3[has_rec, 0], out_rec3[has_rec, 1],
                out_rec3[has_rec, 2], 'p', rot=phi_rot)
    d["acc_w_ep"]  = w_e * w_lead
    d["acc_w_epp"] = w_e * w_lead * w_rec

    d["n_lead_cand"]   = out_nlead
    d["n_recoil_cand"] = out_nrec
    d["n_recoil_hard"] = out_nhard
    d["epp_extra"] = np.zeros(n, dtype=bool)
    _free_raw_columns(d)

    # ---- --epp-recoil all: one extra (e,e'pp)-only row per additional
    # recoil candidate. Every per-event column is copied from the parent
    # event, then the recoil columns are replaced. Extra rows are kept out of
    # (e,e'p) twice over: acc_mask_ep is False and select_eep checks epp_extra.
    if epp_recoil == "all" and xr_l:
        xr = np.concatenate(xr_l)
        x3 = np.concatenate(x3_l)
        nx = len(xr)
        for k in list(d.keys()):
            v = d[k]
            if k.startswith("_"):
                continue
            if not isinstance(v, np.ndarray):
                if not hasattr(v, "__len__") or len(v) != n:
                    continue
                v = np.asarray(v)
            if v.ndim >= 1 and v.shape[0] == n:
                d[k] = np.concatenate([v, v[xr]], axis=0)
        sl = slice(n, n + nx)
        d["recoil_post"][sl] = np.column_stack(
            [x3, np.sqrt((x3 ** 2).sum(axis=1) + M_P ** 2)])
        d["rec_type"][sl] = np.concatenate(xt_l).astype(d["rec_type"].dtype)
        d["rec_role"][sl] = np.concatenate(xro_l)
        d["rec_gen"][sl]  = np.concatenate(xg_l)
        d["rec_fate"][sl] = np.concatenate(xf_l)
        d["acc_mask_ep"][sl] = False
        d["epp_extra"][sl] = True
        if acc is not None:
            d["acc_w_epp"][sl] = d["acc_w_ep"][sl] * acc.map_weight(
                x3[:, 0], x3[:, 1], x3[:, 2], 'p', rot=phi_rot)
        compute_kinematics(d, Ebeam=ebeam)
    return d


# ===================================================================
#  Event selection
# ===================================================================

def select_eep(d):
    """(e,e'p): at least one accepted proton (the lead). Rows appended by
    --epp-recoil all are extra recoil pairings of an event already counted,
    so they never enter (e,e'p)."""
    m = d["lead_type"] == PROTON
    if "epp_extra" in d:
        m = m & ~d["epp_extra"]
    return m


def _primary(d):
    """One row per generated event: False on rows added by --epp-recoil all."""
    if "epp_extra" in d:
        return ~d["epp_extra"]
    return np.ones(len(d["weight"]), dtype=bool)


def _print_recoil_multiplicity(d):
    """Distribution of the number of possible recoils per event: accepted
    protons other than the lead with |p| >= precoil_min. One entry per event,
    weighted by generator weight x electron and lead map weights."""
    base = select_eep(d) & d["acc_mask_ep"]
    w = d["weight"] * _acc_w_ep(d)
    k = d["n_recoil_hard"]
    for label, sel in (("accepted e + lead", base),
                       ("+ SRC cuts on e and lead", base & apply_src_cuts(d))):
        W = w[sel].sum()
        parts = []
        for lab, kk in (("0", k == 0), ("1", k == 1), ("2", k == 2),
                        ("3", k == 3), (">=4", k >= 4)):
            s_ = sel & kk
            frac = 100 * w[s_].sum() / W if W > 0 else 0.0
            parts.append(f"{lab}: {int(s_.sum())} ({frac:.2f}%)")
        print(f"    possible recoils per event [{label}]   " + "   ".join(parts))


def select_eepp(d):
    """(e,e'pp): at least two accepted protons (lead and recoil)."""
    return (d["lead_type"] == PROTON) & (d["rec_type"] == PROTON)


_ACTIVE_SRC_CUTS = SRC_CUTS   # swapped to SRC_CUTS_WRIGHT by --convention


def apply_src_cuts(d):
    """Return boolean mask for SRC selection cuts (paper page 4)."""
    c = _ACTIVE_SRC_CUTS
    mask = (
        (d["xB"] > c["xB_min"]) &
        (d["theta_pq_deg"] < c["theta_pq_max_deg"]) &
        (d["pN_over_q"] > c["pN_over_q_min"]) &
        (d["pN_over_q"] < c["pN_over_q_max"]) &
        (d["pmiss"] > c["pmiss_min"]) &
        (d["pmiss"] < c["pmiss_max"]) &
        (d["mmiss"] < c["mmiss_max"])
    )
    return mask

# ===================================================================
#  T+SCX reweighting  (Eq. 3 of paper)
# ===================================================================

def tscx_weight_eep(d):
    """Per-event T+SCX weight for (e,e'p) channel using Natalie's parameters.

    Deterministic equivalent of Natalie's do_SXC() event-by-event flip:
      w = P(final lead = proton) * T_N
    where P(final lead = p) depends on the initial (lead, recoil) types and
    sums over no-flip + recoil-only-flip (if initial lead = p) or
    lead-only-flip + both-flip (if initial lead = n).
    """
    w = np.zeros(len(d["weight"]))
    lt = d.get("gen_lead_type", d["lead_type"])
    rt = d.get("gen_rec_type",  d["rec_type"])
    pp = (lt == PROTON)  & (rt == PROTON)
    pn = (lt == PROTON)  & (rt == NEUTRON)
    np_= (lt == NEUTRON) & (rt == PROTON)
    nn = (lt == NEUTRON) & (rt == NEUTRON)
    # P(final lead = p):
    #   (p,p): no flip OR recoil-only flip  = 1 - PP2NP - PP2NN  = 0.9561
    #   (p,n): no flip OR recoil-only flip  = 1 - PN2NN - PN2NP  = 0.9629
    #   (n,p): lead-only flip OR both flip  = NP2PP + NP2PN      = 0.0371
    #   (n,n): lead-only flip OR both flip  = NN2PP + NN2PN      = 0.0439
    T_N = T_SCX["T_N"]
    w[pp]  = (1.0 - T_SCX["PP2NP"] - T_SCX["PP2NN"]) * T_N
    w[pn]  = (1.0 - T_SCX["PN2NN"] - T_SCX["PN2NP"]) * T_N
    w[np_] = (T_SCX["NP2PP"] + T_SCX["NP2PN"])       * T_N
    w[nn]  = (T_SCX["NN2PP"] + T_SCX["NN2PN"])       * T_N
    return w


def tscx_weight_eepp(d):
    """Per-event T+SCX weight for (e,e'pp) channel using Natalie's parameters.

      w = P(final = (p, p)) * T_NN
    """
    w = np.zeros(len(d["weight"]))
    lt = d.get("gen_lead_type", d["lead_type"])
    rt = d.get("gen_rec_type",  d["rec_type"])
    pp = (lt == PROTON)  & (rt == PROTON)
    pn = (lt == PROTON)  & (rt == NEUTRON)
    np_= (lt == NEUTRON) & (rt == PROTON)
    # P(final = (p,p)):
    #   (p,p): no flip = 1 - PP2NN - PP2NP - PP2PN  = 0.9081
    #   (p,n): recoil flips n->p  = PN2PP            = 0.041
    #   (n,p): lead   flips n->p  = NP2PP            = 0.035
    #   (n,n): both-flip nn->pp is OMITTED, matching Eq. 3 of the paper
    #          (three terms only; no sigma_(e,e'nn) contribution).
    T_NN = T_SCX["T_NN"]
    w[pp]  = (1.0 - T_SCX["PP2NN"] - T_SCX["PP2NP"] - T_SCX["PP2PN"]) * T_NN
    w[pn]  = T_SCX["PN2PP"] * T_NN
    w[np_] = T_SCX["NP2PP"] * T_NN
    return w


# ===================================================================
#  CLAS data loading (pre-made paper histograms)
# ===================================================================
#
# Per Natalie Mail 2: Acceptance/SRC_eg2_Hists_C.root "contains premade
# histograms that are exactly those used in the original paper".
# Each entry below maps a paper figure panel to its data histogram or
# TGraphAsymmErrors in that file.

_CLAS_PANEL_MAP = {
    # (fig, panel): (hist_name, object_type)
    # panel = 'src' only since the paper's left ("all events") panels
    # show no data points (paper Fig 2 caption: "left panel y-axis scale
    # is arbitrary").
    ("fig2",       "src"): ("ep_mom1",       "TH1"),
    ("fig3_lead",  "src"): ("epp_mom1",      "TH1"),
    ("fig3_rec",   "src"): ("epp_mom2",      "TH1"),
    ("fig4_ep",    "src"): ("ep_Pm_graph",   "TGraph"),
    ("fig4_epp",   "src"): ("epp_Pm_graph",  "TGraph"),
    # Paper Fig 5 uses coarse binning (9 pts, y<=0.14); the fine
    # 30-point `pp_to_p` has a noisy tail up to 0.27.
    ("fig5_ratio", "src"): ("pp_to_p_coarse", "TGraph"),
    # Fig 6 Pcm components in the pmiss frame (z || pmiss, q in x-z plane).
    # Each is stored only as 36 sub-binned sub-histograms `*_i_j_k`. Summing
    # over i=0..3 with j=k=0 reproduces the full (e,e'pp) SRC-cut sample
    # (sub-bins were verified: j=0 is the "all" slice, k=0 is the "all"
    # slice, so summing i for fixed j=0,k=0 gives the total 364-event
    # (e,e'pp) sample).  The 4th column of Fig 6 (|p_cm|) has no
    # matching 1D histogram in the file, so we leave it sim-only.
    ("fig6_x",     "src"): ("epp_Pcmn1m",    "TH1_sumi"),
    ("fig6_y",     "src"): ("epp_Pcmn2m",    "TH1_sumi"),
    ("fig6_z",     "src"): ("epp_Pcmzm",     "TH1_sumi"),
    # Fig 7 p_rel (e,e'pp) SRC-cut: paper uses `epp_km`, NOT `epp_k`.
    # Both have N≈362; `epp_k` truncates to 0 above 0.62 GeV/c
    # (different definition — possibly CM-frame), while `epp_km`
    # extends smoothly to 0.78 GeV/c, matching the paper Fig 7
    # right-panel tail shape.
    ("fig7_prel",  "src"): ("epp_km",        "TH1"),
    # Fig 9 theta(p_miss, q) for (e,e'p) and (e,e'pp) after SRC cuts.
    # `ep_Pmq` covers 100-180 deg (5604 events), `epp_Pmq` 100-180 deg (364).
    ("fig9_ep",    "src"): ("ep_Pmq",        "TH1"),
    ("fig9_epp",   "src"): ("epp_Pmq",       "TH1"),
}


def _edges_from_centers(x):
    """Reconstruct bin edges from an array of bin centers (TGraph).

    Inner edges are midpoints between consecutive centers; the outer
    edges are extrapolated so that the first/last bin has the same
    half-width as its immediate neighbour. Works for both uniform and
    mildly non-uniform spacing.
    """
    x = np.asarray(x, dtype=float)
    mid = 0.5 * (x[:-1] + x[1:])
    first = x[0]  - (mid[0]  - x[0])
    last  = x[-1] + (x[-1]  - mid[-1])
    return np.concatenate([[first], mid, [last]])


def _read_data_hist(rootfile, key, kind):
    """Return (x, y, yerr_lo, yerr_hi, edges) for a single data entry.

    For TH1 we use the native bin edges. For TGraphs we reconstruct
    plausible bin edges from the point centers so the sim can be
    histogrammed on identical bins (crucial for area normalization —
    without this, sim and data would sit at different bin heights
    whenever the number of sim bins differs from the number of data
    points).
    """
    # TH1_sumi: histogram is only stored as 36 sub-histograms `<key>_i_j_k`.
    # Summing over i=0..3 with j=k=0 recovers the aggregate (e,e'pp) SRC sample
    # (sub-bin structure: i is an exclusive 4-bin partition, j=0/k=0 are the
    # "all events" slices along the two inclusive axes).
    if kind == "TH1_sumi":
        axis_edges = None
        total = None
        for i in range(4):
            h = rootfile[f"{key}_{i}_0_0"]
            v = h.values().astype(float)
            if axis_edges is None:
                axis_edges = h.axis(0).edges()
            total = v if total is None else total + v
        edges = axis_edges
        centers = 0.5 * (edges[:-1] + edges[1:])
        err = np.sqrt(np.clip(total, 0, None))
        return centers, total, err, err, edges

    obj = rootfile[key]
    if kind == "TH1":
        edges = obj.axis(0).edges()
        centers = 0.5 * (edges[:-1] + edges[1:])
        counts = obj.values()
        # Use sqrt(N) for data counts (Poisson), matching paper figures
        err = np.sqrt(np.clip(counts, 0, None))
        return centers, counts.astype(float), err, err, edges
    elif kind == "TGraph":
        x = obj.values("x")
        y = obj.values("y")
        # TGraphAsymmErrors: .errors('low'/'high') returns (x_err, y_err)
        try:
            _, yerr_lo = obj.errors("low")
            _, yerr_hi = obj.errors("high")
        except Exception:
            # TGraphErrors (symmetric): one-shot call
            try:
                _, yerr = obj.errors()
                yerr_lo = yerr_hi = yerr
            except Exception:
                yerr_lo = yerr_hi = np.zeros_like(y)
        edges = _edges_from_centers(x) if len(x) >= 2 else None
        return x, y, np.asarray(yerr_lo), np.asarray(yerr_hi), edges
    else:
        raise ValueError(f"Unknown data kind: {kind}")


def load_paper_extracted(filepath="analysis/paper_extracted_data.json"):
    """Load data points extracted from the paper PDF's vector graphics
    (extract_paper_points.py).  These are the exact points Natalie plotted —
    verified against the Hists file on fig 9 (median |dy| < 0.3 counts).

    Used for figs 6/7/8, where the paper drew the Ref. [21] (Cohen et al.,
    PRL 121, 092501, doi:10.1103/PhysRevLett.121.092501) (e,e'pp) sample
    at finer binning (~0.04 GeV/c) than any histogram stored in the Hists
    file (30 bins over [-1,1]), and for the |p_cm| column which has no
    stored 1D histogram at all.

    Returns { key: data_entry } with the same structure load_clas_hists()
    produces, or None if the JSON is missing.  Bin edges are rebuilt on the
    uniform grid implied by the point spacing (the paper omits zero bins,
    so edges cannot come from _edges_from_centers directly).
    """
    if not os.path.isfile(filepath):
        return None
    import json
    with open(filepath) as fh:
        raw = json.load(fh)
    out = {}
    for key, pts in raw.items():
        if not pts:
            continue
        x = np.array([p["x"] for p in pts])
        order = np.argsort(x)
        x = x[order]
        y = np.array([p["y"] for p in pts])[order]
        elo = np.array([p["yerr_lo"] for p in pts])[order]
        ehi = np.array([p["yerr_hi"] for p in pts])[order]
        if len(x) >= 2:
            s = float(np.median(np.diff(x)))
            n = int(round((x[-1] - x[0]) / s)) + 1
            edges = x[0] - s / 2 + s * np.arange(n + 1)
        else:
            edges = None
        out[key] = dict(x=x, y=y, yerr_lo=elo, yerr_hi=ehi, edges=edges,
                        hist_name=key, kind="paper_pdf")
    return out


def load_clas_hists(filepath="Acceptance/SRC_eg2_Hists_C.root"):
    """Load all paper data histograms/graphs needed for the figures.

    Returns a dict: { (fig, panel): { "x","y","yerr_lo","yerr_hi","edges" } }.
    Returns None if the file is missing.
    """
    if not os.path.isfile(filepath):
        return None
    f = uproot.open(filepath)
    out = {}
    for key, (hname, kind) in _CLAS_PANEL_MAP.items():
        try:
            x, y, yerr_lo, yerr_hi, edges = _read_data_hist(f, hname, kind)
            out[key] = dict(x=x, y=y, yerr_lo=yerr_lo, yerr_hi=yerr_hi,
                            edges=edges, hist_name=hname, kind=kind)
        except Exception as e:
            print(f"  [WARN] could not load {hname}: {e}")
    return out


# ===================================================================
#  CLAS raw-event data  (SRC_e2p_C_GoodRuns_coulomb.root)
# ===================================================================
#
# The eg2 "Hists" file stores pre-made 1D hists/graphs at fixed, sometimes
# coarse binning (e.g. `epp_Pcm*_i_0_0` are 30-bin 0.067-wide hists over
# [-1,1] although the paper shows the same events at finer binning).  To
# reproduce the paper's appearance we re-bin the (e,e'pp) data from the
# raw event file at any resolution we like.  (The (e,e'p) sample has
# ~5600 events in a separate ntuple not included here, so we keep the
# pre-made Hists file for the (e,e'p) figures.)

def load_clas_events(filepath="data/clas/SRC_e2p_C_GoodRuns_coulomb.root"):
    """Load raw CLAS (e,e'pp) events, compute all observables, mark SRC cuts.

    The file contains 605 measured (e,e'pp) events, one TTree `T`.  For
    each event we compute (pLead, pRec, pmiss, pcm_{x,y,z,mag}, prel, mmiss,
    xB, theta_pq, pN/|q|, Q2) and two boolean masks — `src_mask_ep` (all
    SRC cuts on the leading proton) and `src_mask_epp` (the same plus the
    |precoil|>0.35 GeV/c cut).  p_cm components are in the Ref. [21]
    (Cohen et al., PRL 121, 092501, doi:10.1103/PhysRevLett.121.092501)
    pmiss-frame (z along pmiss_hat, q in the x-z plane) — same convention
    the simulation uses for pcm_x/y/z.
    """
    if not os.path.isfile(filepath):
        return None

    T = uproot.open(filepath)["T"]
    arrs = T.arrays()

    Pe   = np.asarray(arrs["Pe"])          # (N, 3)
    Pp   = np.asarray(arrs["Pp"])          # (N, 2, 3) two protons per event
    qvec = np.asarray(arrs["q"])           # (N, 3)
    Pmiss_all = np.asarray(arrs["Pmiss"])  # (N, 2, 3)  — [:,0] is lead pmiss
    Nu   = np.asarray(arrs["Nu"])
    Xb   = np.asarray(arrs["Xb"])
    Q2   = np.asarray(arrs["Q2"])
    pq_angle_lead = np.asarray(arrs["pq_angle"])[:, 0]   # degrees

    p_lead = Pp[:, 0]
    p_rec  = Pp[:, 1]
    pN_mag   = np.linalg.norm(p_lead, axis=1)
    prec_mag = np.linalg.norm(p_rec,  axis=1)
    q_mag    = np.linalg.norm(qvec,   axis=1)

    pmiss_vec = Pmiss_all[:, 0]
    pmiss_mag = np.linalg.norm(pmiss_vec, axis=1)

    # 2-body missing mass (struck SRC pair approximation):
    # m_miss^2 = (2 m_avg + nu - E_lead)^2 - |p_miss|^2
    Ep1 = np.sqrt(pN_mag**2 + M_P**2)
    mmiss2 = (2.0 * M_AVG + Nu - Ep1)**2 - pmiss_mag**2
    mmiss  = np.where(mmiss2 > 0, np.sqrt(mmiss2), -np.sqrt(-mmiss2))

    pN_over_q = pN_mag / np.where(q_mag > 0, q_mag, 1e-30)

    # c.m. momentum of the SRC pair: p_cm = p_lead + p_rec - q = p_miss + p_rec
    pcm_vec = p_lead + p_rec - qvec
    pcm_mag = np.linalg.norm(pcm_vec, axis=1)

    # relative momentum: p_rel = (p_miss - p_rec)/2
    prel_vec = 0.5 * (pmiss_vec - p_rec)
    prel_mag = np.linalg.norm(prel_vec, axis=1)

    # Pcm rotation into the pmiss frame  (same derivation as the sim):
    #   z_hat = pmiss_hat;  x_hat = projection of q onto plane perp to z;
    #   y_hat = z x x.
    pmiss_hat  = pmiss_vec / (pmiss_mag[:, None] + 1e-30)
    z_hat      = pmiss_hat
    q_dot_z    = np.sum(qvec * z_hat, axis=1)
    q_perp     = qvec - q_dot_z[:, None] * z_hat
    q_perp_mag = np.linalg.norm(q_perp, axis=1)
    x_hat      = q_perp / (q_perp_mag[:, None] + 1e-30)
    y_hat      = np.cross(z_hat, x_hat)

    pcm_x = np.sum(pcm_vec * x_hat, axis=1)
    pcm_y = np.sum(pcm_vec * y_hat, axis=1)
    pcm_z = np.sum(pcm_vec * z_hat, axis=1)

    # Angle between p_miss and q (used in Fig 9)
    cos_pmiss_q = np.sum(pmiss_vec * qvec, axis=1) / (
        pmiss_mag * q_mag + 1e-30
    )
    cos_pmiss_q = np.clip(cos_pmiss_q, -1, 1)
    theta_pmiss_q_deg = np.degrees(np.arccos(cos_pmiss_q))

    d = dict(
        pN=pN_mag, prec=prec_mag, q3=q_mag, pmiss=pmiss_mag, pmiss_vec=pmiss_vec,
        theta_pq_deg=pq_angle_lead, mmiss=mmiss, pN_over_q=pN_over_q,
        pcm_vec=pcm_vec, pcm_mag=pcm_mag,
        pcm_x=pcm_x, pcm_y=pcm_y, pcm_z=pcm_z,
        prel_mag=prel_mag, xB=Xb, Q2=Q2, Nu=Nu, qvec=qvec,
        p_lead=p_lead, p_rec=p_rec, Pe=Pe,
        theta_pmiss_q_deg=theta_pmiss_q_deg,
    )

    c = SRC_CUTS
    src_mask_lead = (
        (Xb > c["xB_min"]) &
        (pq_angle_lead < c["theta_pq_max_deg"]) &
        (pN_over_q > c["pN_over_q_min"]) &
        (pN_over_q < c["pN_over_q_max"]) &
        (pmiss_mag > c["pmiss_min"]) &
        (pmiss_mag < c["pmiss_max"]) &
        (mmiss      < c["mmiss_max"])
    )
    d["src_mask_ep"]  = src_mask_lead
    d["src_mask_epp"] = src_mask_lead & (prec_mag > c["precoil_min"])
    return d


def events_to_entry(events, obs_key, mask_key, bins):
    """Return a `data_entry` dict (compatible with `_plot_normalized`) by
    histogramming `events[obs_key][events[mask_key]]` with the given bins.

    Errors are Poisson (sqrt(N)).  Returns None if events is None.
    """
    if events is None:
        return None
    mask = events[mask_key]
    vals = events[obs_key][mask]
    counts, edges = np.histogram(vals, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    err = np.sqrt(np.clip(counts, 0, None))
    return dict(x=centers, y=counts.astype(float),
                yerr_lo=err, yerr_hi=err, edges=edges,
                hist_name=f"(raw events: {obs_key}|{mask_key})")


# ===================================================================
#  Mean-field background
# ===================================================================

# ===================================================================
#  Plotting helpers
# ===================================================================

# Style matching Natalie's paper (ROOT-like appearance):
#   PWIA  = gray dotted       (thin)
#   T+SCX = light blue dashed (medium)
#   FULL  = medium blue solid  (thicker)
#   Data  = black filled circles with error bars
_BLUE = "#5588cc"          # desaturated steel-blue, close to ROOT's kBlue+1
_LBLUE = "#77aadd"         # lighter variant for dashed
_GRAY = "#888888"
_RED  = "#cc3355"          # paper uses a warm color for the MF curve in Fig 8-10

COLORS = dict(pwia=_GRAY, tscx=_LBLUE, full=_BLUE, data="black",
              gcf=_BLUE, mf_full=_RED)
STYLES = dict(pwia=":", tscx="--", full="-",
              gcf="-", mf_full="--")
LABELS = dict(pwia="PWIA", tscx="T+SCX", full="FULL", data="Data",
              gcf="FULL GCF", mf_full="FULL MF")
LINEWIDTHS = dict(pwia=1.2, tscx=1.4, full=1.8,
                  gcf=1.8, mf_full=1.6)

# ROOT-like plot style
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",          # Computer Modern math — matches ROOT's TMathText
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 13,
    "axes.linewidth": 0.8,
    "legend.fontsize": 10,
    "legend.handlelength": 2.5,        # longer legend lines like ROOT
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def make_hist(values, weights, bins):
    """Compute weighted histogram, return (centers, heights)."""
    h, edges = np.histogram(values, bins=bins, weights=weights)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, h


def draw_hist(ax, centers, h, **kwargs):
    """Draw a pre-computed histogram as smooth curve."""
    ax.plot(centers, h, **kwargs)

def area_normalize(histograms, ref_area=None):
    """Normalize a list of (centers, h, style_key) to same area.

    Per Wright et al. (doi:10.48550/arXiv.2104.05090; paper Fig 2 caption): "the calculations are
    individually area normalized to the data". If `ref_area` is provided
    (e.g. the integral of the data histogram), all calculations are
    scaled to match it. Otherwise the first non-zero histogram is used
    as reference (for left panels with no data).
    """
    if ref_area is None:
        # Find the first histogram with nonzero area as reference
        for _, h, _ in histograms:
            a = np.sum(h)
            if a > 0:
                ref_area = a
                break
    if ref_area is None or ref_area == 0:
        return histograms

    normed = []
    for centers, h, key in histograms:
        a = np.sum(h)
        if a > 0:
            normed.append((centers, h * (ref_area / a), key))
        else:
            normed.append((centers, h, key))
    return normed


# (The original `plot_panel` helper has been replaced by `_plot_normalized`
#  which accepts pre-loaded paper histograms via `data_entry`.)


# ===================================================================
#  Figure functions — one per paper figure
# ===================================================================

def _acc_w_ep(d):
    """Per-event CLAS acceptance weight for (e,e'p). 1 if no acceptance applied."""
    return d.get("acc_w_ep", np.ones(len(d["weight"])))


def _acc_w_epp(d):
    """Per-event CLAS acceptance weight for (e,e'pp)."""
    return d.get("acc_w_epp", np.ones(len(d["weight"])))


def _acc_mask_ep(d):
    return d.get("acc_mask_ep", np.ones(len(d["weight"]), dtype=bool))


def _acc_mask_epp(d):
    return d.get("acc_mask_epp", np.ones(len(d["weight"]), dtype=bool))


def _collect_eep(pwia, fsi, var_key, cuts_mask_pwia, cuts_mask_fsi):
    """Collect (values, weights, style_key) for (e,e'p) PWIA/T+SCX/Full.

    Applies CLAS fiducial masks and multiplies by per-event map weights.
    """
    acc_p = _acc_mask_ep(pwia)
    acc_f = _acc_mask_ep(fsi)
    w_p   = _acc_w_ep(pwia)
    w_f   = _acc_w_ep(fsi)

    m_p = (pwia["lead_type"] == PROTON) & cuts_mask_pwia & acc_p
    tscx_w = tscx_weight_eep(pwia)
    m_t = cuts_mask_pwia & acc_p
    m_f = (fsi["lead_type"] == PROTON) & cuts_mask_fsi & acc_f
    return [
        (pwia[var_key][m_p], pwia["weight"][m_p] * w_p[m_p],           "pwia"),
        (pwia[var_key][m_t], pwia["weight"][m_t] * tscx_w[m_t] * w_p[m_t], "tscx"),
        (fsi [var_key][m_f], fsi ["weight"][m_f] * w_f[m_f],           "full"),
    ]


def _collect_eepp(pwia, fsi, var_key, cuts_mask_pwia, cuts_mask_fsi):
    """Collect (values, weights, style_key) for (e,e'pp) PWIA/T+SCX/Full."""
    acc_p = _acc_mask_epp(pwia)
    acc_f = _acc_mask_epp(fsi)
    w_p   = _acc_w_epp(pwia)
    w_f   = _acc_w_epp(fsi)

    rec_pw = cuts_mask_pwia & (pwia["prec"] > SRC_CUTS["precoil_min"]) & acc_p
    rec_fs = cuts_mask_fsi  & (fsi["prec"]  > SRC_CUTS["precoil_min"]) & acc_f

    m_p = (pwia["lead_type"] == PROTON) & (pwia["rec_type"] == PROTON) & rec_pw
    tscx_w = tscx_weight_eepp(pwia)
    m_t = rec_pw
    m_f = (fsi["lead_type"] == PROTON) & (fsi["rec_type"] == PROTON) & rec_fs
    return [
        (pwia[var_key][m_p], pwia["weight"][m_p] * w_p[m_p],           "pwia"),
        (pwia[var_key][m_t], pwia["weight"][m_t] * tscx_w[m_t] * w_p[m_t], "tscx"),
        (fsi [var_key][m_f], fsi ["weight"][m_f] * w_f[m_f],           "full"),
    ]


def _plot_normalized(ax, datasets, bins, xlabel="", title="", logy=False,
                     ylabel="Counts", data_entry=None):
    """Build histograms, area-normalize, draw.

    `data_entry` (optional): dict from load_clas_hists() with keys
        x, y, yerr_lo, yerr_hi, edges.  When supplied:
          * `bins` is overridden by the data's bin edges if they are available,
            so sim and data share identical binning (required for a clean
            area-normalization comparison).
          * All simulation curves are area-normalized to match the data integral.
          * The data is overlaid as black filled circles with (asymmetric) errors.
    """
    # Pin sim binning to the data binning when edges are known
    if data_entry is not None and data_entry.get("edges") is not None:
        bins = data_entry["edges"]

    # Calculation histograms
    hists = []
    for vals, wts, key in datasets:
        c, h = make_hist(vals, wts, bins)
        hists.append((c, h, key))

    ref_area = None
    if data_entry is not None:
        ref_area = float(np.sum(data_entry["y"]))

    hists = area_normalize(hists, ref_area=ref_area)

    for c, h, key in hists:
        draw_hist(ax, c, h, color=COLORS[key], linestyle=STYLES[key],
                  label=LABELS[key], linewidth=LINEWIDTHS.get(key, 1.5))

    # Overlay paper data points (with asymmetric errors if provided)
    if data_entry is not None:
        yerr = np.vstack([data_entry["yerr_lo"], data_entry["yerr_hi"]])
        ax.errorbar(data_entry["x"], data_entry["y"], yerr=yerr,
                    fmt="o", ms=4, color=COLORS["data"], label=LABELS["data"],
                    elinewidth=0.8, capsize=0, zorder=10)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.text(0.95, 0.95, title, transform=ax.transAxes,
                ha="right", va="top", fontsize=11)
    if logy:
        ax.set_yscale("log")
    ax.legend(frameon=False, loc="best")


def _clas(clas, key):
    """Safe lookup in the loaded CLAS-hists dict."""
    return None if clas is None else clas.get(key)


def figure2(pwia, fsi, outdir, clas=None):
    """Fig 2: Leading-proton momentum for 12C(e,e'p).
    Data shown only on right (SRC cuts) panel; sim is area-normalized to
    data (paper Fig 2 caption)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    bins = np.linspace(0.5, 2.5, 31)   # 0.067 GeV/c bins (was 60 -> 30 bins)

    no_cuts = np.ones(len(pwia["weight"]), dtype=bool)
    no_cuts_f = np.ones(len(fsi["weight"]), dtype=bool)
    src_p = apply_src_cuts(pwia)
    src_f = apply_src_cuts(fsi)

    for i, (label, cp, cf, data_key) in enumerate([
        ("All events", no_cuts, no_cuts_f, None),
        ("SRC cuts",   src_p,   src_f,     ("fig2", "src")),
    ]):
        ds = _collect_eep(pwia, fsi, "pN", cp, cf)
        _plot_normalized(axes[i], ds, bins,
                         xlabel=r"$p_{\rm Lead}$ [GeV/c]",
                         title=rf"$^{{12}}$C(e,e'p) — {label}",
                         data_entry=_clas(clas, data_key))
        axes[i].set_xlim(0.5, 2.5)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig2_pLead_eep.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig2_pLead_eep.png")


def figure3(pwia, fsi, outdir, clas=None, clas_ev=None):
    """Fig 3: Lead (top) and recoil (bottom) proton momentum for 12C(e,e'pp).

    Data (SRC-cut column): re-binned from the raw (e,e'pp) event file when
    `clas_ev` is provided; otherwise falls back to the stored Hists-file
    histograms `epp_mom1`/`epp_mom2` (40 / 17 bins).

    Channel: (e,e'pp)."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    bins_lead = np.linspace(0.5, 2.5, 21)   # 0.10 GeV/c bin (was 41)
    bins_rec  = np.linspace(0.35, 1.2, 18)  # 0.05 GeV/c bin (matches epp_mom2)

    no_cuts = np.ones(len(pwia["weight"]), dtype=bool)
    no_cuts_f = np.ones(len(fsi["weight"]), dtype=bool)
    src_p = apply_src_cuts(pwia)
    src_f = apply_src_cuts(fsi)

    for col, (label, cp, cf, use_data) in enumerate([
        ("All events", no_cuts, no_cuts_f, False),
        ("SRC cuts",   src_p,   src_f,     True),
    ]):
        ds_lead = _collect_eepp(pwia, fsi, "pN",   cp, cf)
        ds_rec  = _collect_eepp(pwia, fsi, "prec", cp, cf)

        if use_data:
            # Prefer paper-exact Hists (Natalie Mail 2: "exactly those
            # used in the original paper").  The raw eg2 ntuple has
            # more events with looser cuts and won't match the paper.
            de_lead = (_clas(clas, ("fig3_lead", "src"))
                       or events_to_entry(clas_ev, "pN", "src_mask_epp", bins_lead))
            de_rec  = (_clas(clas, ("fig3_rec",  "src"))
                       or events_to_entry(clas_ev, "prec","src_mask_epp", bins_rec))
        else:
            de_lead = de_rec = None

        _plot_normalized(axes[0, col], ds_lead, bins_lead,
                         xlabel=r"$p_{\rm Lead}$ [GeV/c]",
                         title=rf"$^{{12}}$C(e,e'pp) — {label}",
                         data_entry=de_lead)
        _plot_normalized(axes[1, col], ds_rec, bins_rec,
                         xlabel=r"$p_{\rm Recoil}$ [GeV/c]",
                         data_entry=de_rec)
        axes[0, col].set_xlim(0.5, 2.5)
        axes[1, col].set_xlim(0.38, 1.0)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_pLead_pRec_eepp.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig3_pLead_pRec_eepp.png")


def figure4(pwia, fsi, outdir, clas=None, clas_ev=None):
    """Fig 4: Missing momentum for (e,e'p) and (e,e'pp).

    (e,e'p) data comes from the Hists file (~5600 events, not in the raw
    event ntuple).  (e,e'pp) data is re-binned from the raw event file
    when `clas_ev` is provided."""
    pwia_eep = pwia_eepp = pwia
    fsi_eep = fsi_eepp = fsi
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    bins = np.linspace(0.4, 1.0, 13)   # 0.05 GeV/c bins

    no_cuts_p_eep   = np.ones(len(pwia_eep["weight"]),  dtype=bool)
    no_cuts_p_eepp  = np.ones(len(pwia_eepp["weight"]), dtype=bool)
    no_cuts_f_eep   = np.ones(len(fsi_eep["weight"]),   dtype=bool)
    no_cuts_f_eepp  = np.ones(len(fsi_eepp["weight"]),  dtype=bool)
    src_p_eep  = apply_src_cuts(pwia_eep)
    src_p_eepp = apply_src_cuts(pwia_eepp)
    src_f_eep  = apply_src_cuts(fsi_eep)
    src_f_eepp = apply_src_cuts(fsi_eepp)

    for col, (label, cp_eep, cf_eep, cp_eepp, cf_eepp, use_data) in enumerate([
        ("All events", no_cuts_p_eep, no_cuts_f_eep, no_cuts_p_eepp, no_cuts_f_eepp, False),
        ("SRC cuts",   src_p_eep,     src_f_eep,     src_p_eepp,     src_f_eepp,     True),
    ]):
        ds_eep  = _collect_eep (pwia_eep,  fsi_eep,  "pmiss", cp_eep,  cf_eep)
        ds_eepp = _collect_eepp(pwia_eepp, fsi_eepp, "pmiss", cp_eepp, cf_eepp)

        de_eep = _clas(clas, ("fig4_ep", "src")) if use_data else None
        if use_data:
            de_eepp = (_clas(clas, ("fig4_epp", "src"))
                       or events_to_entry(clas_ev, "pmiss", "src_mask_epp", bins))
        else:
            de_eepp = None

        _plot_normalized(axes[0, col], ds_eep, bins, logy=True,
                         xlabel=r"$p_{\rm Miss}$ [GeV/c]",
                         title=rf"$^{{12}}$C(e,e'p) — {label}",
                         data_entry=de_eep)
        _plot_normalized(axes[1, col], ds_eepp, bins, logy=True,
                         xlabel=r"$p_{\rm Miss}$ [GeV/c]",
                         title=rf"$^{{12}}$C(e,e'pp) — {label}",
                         data_entry=de_eepp)
        for row in range(2):
            axes[row, col].set_xlim(0.4, 1.0)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig4_pmiss.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig4_pmiss.png")


def figure5(pwia, fsi, outdir, clas=None):
    """Fig 5: (e,e'pp)/(e,e'p) ratio vs p_miss after SRC cuts.
    Paper: x=[0.4,1.0], y=[0,0.2], linear. Uses Natalie's pp_to_p TGraph
    for the measured ratio."""
    pwia_eep = pwia_eepp = pwia
    fsi_eep = fsi_eepp = fsi
    fig, ax = plt.subplots(figsize=(6, 5))

    data_entry = _clas(clas, ("fig5_ratio", "src"))
    if data_entry is not None and data_entry.get("edges") is not None:
        bins = data_entry["edges"]
    elif data_entry is not None and data_entry.get("x") is not None:
        xd = data_entry["x"]
        step = xd[1] - xd[0]
        bins = np.concatenate([[xd[0] - step/2], xd + step/2])
    else:
        bins = np.linspace(0.4, 1.0, 14)
    centers = 0.5 * (bins[:-1] + bins[1:])
    src = apply_src_cuts

    def _ratio(d_eep, d_eepp, sel_eep, sel_eepp,
               w_scale_eep, w_scale_eepp):
        """Build weighted pmiss hists for (e,e'p) and (e,e'pp), return ratio."""
        w_ep  = d_eep ["weight"] * w_scale_eep
        w_epp = d_eepp["weight"] * w_scale_eepp
        h_ep,  _ = np.histogram(d_eep ["pmiss"][sel_eep ], bins=bins, weights=w_ep [sel_eep ])
        h_epp, _ = np.histogram(d_eepp["pmiss"][sel_eepp], bins=bins, weights=w_epp[sel_eepp])
        return np.divide(h_epp, h_ep,
                         out=np.zeros_like(h_epp, dtype=float),
                         where=h_ep > 0)

    # ------- PWIA -------
    m_eep_pw  = select_eep(pwia_eep)  & src(pwia_eep)  & _acc_mask_ep(pwia_eep)
    m_eepp_pw = (select_eepp(pwia_eepp) & src(pwia_eepp) & _acc_mask_epp(pwia_eepp)
                 & (pwia_eepp["prec"] > SRC_CUTS["precoil_min"]))
    r_pwia = _ratio(pwia_eep, pwia_eepp, m_eep_pw, m_eepp_pw,
                    _acc_w_ep(pwia_eep), _acc_w_epp(pwia_eepp))
    ax.plot(centers, r_pwia, color=COLORS["pwia"], linestyle=STYLES["pwia"],
            label=LABELS["pwia"], linewidth=LINEWIDTHS["pwia"])

    # ------- Full FSI -------
    m_eep_fs  = select_eep(fsi_eep)  & src(fsi_eep)  & _acc_mask_ep(fsi_eep)
    m_eepp_fs = (select_eepp(fsi_eepp) & src(fsi_eepp) & _acc_mask_epp(fsi_eepp)
                 & (fsi_eepp["prec"] > SRC_CUTS["precoil_min"]))
    r_full = _ratio(fsi_eep, fsi_eepp, m_eep_fs, m_eepp_fs,
                    _acc_w_ep(fsi_eep), _acc_w_epp(fsi_eepp))
    ax.plot(centers, r_full, color=COLORS["full"], linestyle=STYLES["full"],
            label=LABELS["full"], linewidth=LINEWIDTHS["full"])

    # ------- T+SCX -------
    tscx_eep  = tscx_weight_eep (pwia_eep)
    tscx_eepp = tscx_weight_eepp(pwia_eepp)
    m_all_pw  = src(pwia_eep)  & _acc_mask_ep (pwia_eep)
    m_rec_pw  = (src(pwia_eepp) & _acc_mask_epp(pwia_eepp)
                 & (pwia_eepp["prec"] > SRC_CUTS["precoil_min"]))
    r_tscx = _ratio(pwia_eep, pwia_eepp, m_all_pw, m_rec_pw,
                    tscx_eep  * _acc_w_ep (pwia_eep),
                    tscx_eepp * _acc_w_epp(pwia_eepp))
    ax.plot(centers, r_tscx, color=COLORS["tscx"], linestyle=STYLES["tscx"],
            label=LABELS["tscx"], linewidth=LINEWIDTHS["tscx"])

    # ------- data -------
    if data_entry is not None:
        yerr = np.vstack([data_entry["yerr_lo"], data_entry["yerr_hi"]])
        ax.errorbar(data_entry["x"], data_entry["y"], yerr=yerr,
                    fmt="o", ms=4, color=COLORS["data"], label=LABELS["data"],
                    elinewidth=0.8, capsize=0, zorder=10)

    ax.set_xlabel(r"$p_{\rm Miss}$ [GeV/c]")
    ax.set_ylabel(r"A(e,e$'$pp) / A(e,e$'$p)")
    ax.text(0.95, 0.95, r"$^{12}$C — SRC cuts", transform=ax.transAxes,
            ha="right", va="top", fontsize=11)
    ax.set_xlim(0.4, 1.0)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig5_ratio_eepp_eep.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig5_ratio_eepp_eep.png")


def figure5_components(pwia, fsi, outdir):
    """Diagnostic behind Fig 5: the four yield spectra that make the ratio.

    Fig 5 shows (e,e'pp)/(e,e'p) vs p_miss, which folds four numbers per bin
    into one. Here they are drawn separately:

      top    -- (e,e'p) and (e,e'pp) yields vs p_miss for PWIA and FULL, all
                four rescaled to the SAME experimental luminosity (the
                attempts-based scale used by print_effective_transparencies),
                so their absolute heights are directly comparable.
      bottom -- the per-bin effective transparency FULL/PWIA for each channel.

    Reading it: a Fig-5 ratio that is too high can come from the (e,e'pp)
    numerator being too large or the (e,e'p) denominator being too small, and
    only the separated yields say which.
    """
    lo, hi = _ACTIVE_SRC_CUTS["pmiss_min"], _ACTIVE_SRC_CUTS["pmiss_max"]
    bins = np.linspace(lo, hi, 13)
    centers = 0.5 * (bins[:-1] + bins[1:])
    src = apply_src_cuts

    def _n_attempts(d):
        meta = d.get("_meta") or {}
        try:
            return float(meta.get("n_events_attempts"))
        except (TypeError, ValueError):
            return None
    att_p, att_f = _n_attempts(pwia), _n_attempts(fsi)
    if att_p and att_f:
        lumi = att_p / att_f
    else:
        lumi = (pwia.get("weight_gen", pwia["weight"])[_primary(pwia)].sum()
                / fsi.get("weight_gen", fsi["weight"])[_primary(fsi)].sum())

    def _yields(d, channel, wscale):
        if channel == "ep":
            m = select_eep(d) & src(d) & _acc_mask_ep(d)
            w = d["weight"] * _acc_w_ep(d) * wscale
        else:
            m = (select_eepp(d) & src(d) & _acc_mask_epp(d)
                 & (d["prec"] > SRC_CUTS["precoil_min"]))
            w = d["weight"] * _acc_w_epp(d) * wscale
        h, _ = np.histogram(d["pmiss"][m], bins=bins, weights=w[m])
        h2, _ = np.histogram(d["pmiss"][m], bins=bins, weights=(w[m])**2)
        return h, np.sqrt(h2)

    y = {}
    y[("pwia", "ep")]  = _yields(pwia, "ep",  1.0)
    y[("pwia", "epp")] = _yields(pwia, "epp", 1.0)
    y[("full", "ep")]  = _yields(fsi,  "ep",  lumi)
    y[("full", "epp")] = _yields(fsi,  "epp", lumi)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.4, 7.2), sharex=True,
        gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.07))

    spec = [
        (("pwia", "ep"),  _GRAY, ":",  r"PWIA  (e,e$'$p)"),
        (("full", "ep"),  _BLUE, "-",  r"FULL  (e,e$'$p)"),
        (("pwia", "epp"), _GRAY, "-.", r"PWIA  (e,e$'$pp)"),
        (("full", "epp"), _RED,  "--", r"FULL  (e,e$'$pp)"),
    ]
    for key, col, ls, lab in spec:
        h, e = y[key]
        ax.errorbar(centers, h, yerr=e, color=col, linestyle=ls, lw=1.7,
                    label=lab, elinewidth=0.7, capsize=0, marker="")

    ax.set_yscale("log")
    ax.set_ylabel("yield at matched luminosity [arb.]")
    ax.legend(frameon=False, loc="upper right", ncol=1)
    ax.text(0.03, 0.05, r"$^{12}$C — SRC cuts", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=11)

    for chan, col, lab in (("ep", _BLUE, r"(e,e$'$p)"),
                           ("epp", _RED, r"(e,e$'$pp)")):
        hp, ep_ = y[("pwia", chan)]
        hf, ef = y[("full", chan)]
        r = np.divide(hf, hp, out=np.zeros_like(hf), where=hp > 0)
        rerr = np.zeros_like(r)
        ok = (hp > 0) & (hf > 0)
        rerr[ok] = r[ok] * np.sqrt((ef[ok]/hf[ok])**2 + (ep_[ok]/hp[ok])**2)
        axr.errorbar(centers, r, yerr=rerr, color=col, lw=1.7, label=lab,
                     elinewidth=0.7, capsize=0)
    axr.axhline(1.0, color="k", lw=0.6, ls=":")
    axr.axhline(0.58, color=_BLUE, lw=0.8, ls="--", alpha=0.6)
    axr.axhline(0.73, color=_RED, lw=0.8, ls="--", alpha=0.6)
    axr.set_xlabel(r"$p_{\rm Miss}$ [GeV/c]")
    axr.set_ylabel("FULL / PWIA")
    axr.set_xlim(lo, hi)
    axr.set_ylim(bottom=0)
    axr.legend(frameon=False, loc="upper left", ncol=2)
    axr.text(0.98, 0.06,
             "dashed: paper integrated values 0.58 / 0.73",
             transform=axr.transAxes, ha="right", va="bottom", fontsize=8,
             color="0.35")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figD_fig5_components.png"), dpi=200)
    plt.close(fig)
    print("  Saved figD_fig5_components.png")

    print("  Fig5 components (matched luminosity):")
    print(f"    {'pmiss':>8} {'ep PWIA':>11} {'ep FULL':>11} {'T_ep':>7}"
          f" {'epp PWIA':>11} {'epp FULL':>11} {'T_epp':>7}")
    for i, c in enumerate(centers):
        hp_ep = y[("pwia", "ep")][0][i];  hf_ep = y[("full", "ep")][0][i]
        hp_pp = y[("pwia", "epp")][0][i]; hf_pp = y[("full", "epp")][0][i]
        t_ep = hf_ep/hp_ep if hp_ep > 0 else float("nan")
        t_pp = hf_pp/hp_pp if hp_pp > 0 else float("nan")
        print(f"    {c:8.3f} {hp_ep:11.4g} {hf_ep:11.4g} {t_ep:7.3f}"
              f" {hp_pp:11.4g} {hf_pp:11.4g} {t_pp:7.3f}")


def figure6(pwia, fsi, outdir, clas=None, clas_ev=None):
    """Fig 6: CM momentum components for 12C(e,e'pp) in the pmiss frame
    (z-axis along pmiss, q in x-z plane — Ref. [21] convention; Ref. [21] is
    Cohen et al., PRL 121, 092501, doi:10.1103/PhysRevLett.121.092501).


    Paper x-ranges: pcm_x,y in [-0.5,0.5], pcm_z in [-0.2,0.8],
    |pcm| in [0,0.8]. Data (bottom row, SRC cuts) is re-binned from the
    raw (e,e'pp) event file when `clas_ev` is provided.  Binning matches
    the paper (~20 bins per column)."""
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))

    # (axis label, sim-key, xlim, bins)
    components = [
        (r"$p_{c.m.}^x$ [GeV/c]", "pcm_x",   (-0.5, 0.5), np.linspace(-0.5, 0.5, 16)),
        (r"$p_{c.m.}^y$ [GeV/c]", "pcm_y",   (-0.5, 0.5), np.linspace(-0.5, 0.5, 16)),
        (r"$p_{c.m.}^z$ [GeV/c]", "pcm_z",   (-0.2, 0.8), np.linspace(-0.2, 0.8, 16)),
        (r"$|p_{c.m.}|$ [GeV/c]", "pcm_mag", ( 0.0, 0.8), np.linspace( 0.0, 0.8, 16)),
    ]

    no_cuts = np.ones(len(pwia["weight"]), dtype=bool)
    no_cuts_f = np.ones(len(fsi["weight"]), dtype=bool)
    src_p = apply_src_cuts(pwia)
    src_f = apply_src_cuts(fsi)

    for row, (label, cp, cf, use_data) in enumerate([
        ("All events", no_cuts, no_cuts_f, False),
        ("SRC cuts",   src_p,   src_f,     True),
    ]):
        for col_idx, (comp_label, comp_key, xlim, bins) in enumerate(components):
            ax = axes[row, col_idx]
            ds = _collect_eepp(pwia, fsi, comp_key, cp, cf)
            title = rf"$^{{12}}$C(e,e'pp) — {label}" if col_idx == 0 else ""
            # Paper-PDF extracted points (exact points the paper shows);
            # all 4 columns including |pcm| now have real data.  The raw
            # events file is NOT used here — it is a different dataset
            # (Ref [7]-style selection, 605 events with a harder pmiss
            # spectrum), not the Ref [21] c.m.-motion sample.
            clas_key = (("fig6_x","src"), ("fig6_y","src"),
                        ("fig6_z","src"), ("fig6_mag","src"))[col_idx]
            de = _clas(clas, clas_key) if use_data else None
            _plot_normalized(ax, ds, bins,
                             xlabel=comp_label, title=title,
                             data_entry=de)
            ax.set_xlim(xlim)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig6_pcm_components.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig6_pcm_components.png")


def figure7(pwia, fsi, outdir, clas=None, clas_ev=None):
    """Fig 7: Relative momentum for 12C(e,e'pp).

    Data (SRC-cut column) is re-binned from the raw (e,e'pp) event file
    when `clas_ev` is provided; otherwise uses the stored `epp_k` TH1D.
    Channel: (e,e'pp)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    bins = np.linspace(0.2, 1.0, 17)   # 0.05 GeV/c bins (was 21)

    no_cuts = np.ones(len(pwia["weight"]), dtype=bool)
    no_cuts_f = np.ones(len(fsi["weight"]), dtype=bool)
    src_p = apply_src_cuts(pwia)
    src_f = apply_src_cuts(fsi)

    for col, (label, cp, cf, use_data) in enumerate([
        ("All events", no_cuts, no_cuts_f, False),
        ("SRC cuts",   src_p,   src_f,     True),
    ]):
        ds = _collect_eepp(pwia, fsi, "prel_mag", cp, cf)
        if use_data:
            de = (_clas(clas, ("fig7_prel", "src"))
                  or events_to_entry(clas_ev, "prel_mag", "src_mask_epp", bins))
        else:
            de = None
        _plot_normalized(axes[col], ds, bins,
                         xlabel=r"$p_{\rm REL}$ [GeV/c]",
                         title=rf"$^{{12}}$C(e,e'pp) — {label}",
                         data_entry=de)
        axes[col].set_xlim(0.2, 1.0)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig7_prel_eepp.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig7_prel_eepp.png")


# --- Figures 8-10: GCF (p_rel > kF, AV18) vs FULL MF (p_rel < kF, Fermi gas) ---
#
# The paper runs these comparisons with the SRC event selection cuts applied
# and only the FULL (post-FSI) calculation shown for each. For us:
#   * FULL GCF   <- fsi_src  (wf_mode = SRC_only, AV18 for p_rel > kF)
#   * FULL MF    <- fsi_mf   (wf_mode = MF_only,  Fermi gas for p_rel < kF)
# Both go through the same CLAS acceptance + smearing pipeline. Each curve is
# individually area-normalized to the data (where data exists).

def _fsi_vals_wts(d, var_key, channel, cuts_mask):
    """Return (values, weights) for one sample in a channel ('ep'/'epp').

    Applies the flat-list channel selection, CLAS fiducial mask, per-event
    acceptance map weight and -- for 'epp' -- the precoil cut. `cuts_mask`
    should already incorporate apply_src_cuts(d). Identical treatment for
    the GCF and MF samples (both are flat-list resolved upstream).
    """
    if channel == "ep":
        m = select_eep(d) & cuts_mask & _acc_mask_ep(d)
        return d[var_key][m], d["weight"][m] * _acc_w_ep(d)[m]
    m = (select_eepp(d) & cuts_mask & _acc_mask_epp(d)
         & (d["prec"] > SRC_CUTS["precoil_min"]))
    return d[var_key][m], d["weight"][m] * _acc_w_epp(d)[m]


def _draw_gcf_mf(ax, fsi_src, fsi_mf, var_key, bins, *, channel,
                 xlabel="", title="", data_entry=None, logy=False,
                 ylabel="Counts", xlim=None):
    """Plot FULL-GCF and FULL-MF overlays on the same axis after SRC cuts.

    Each curve is individually area-normalized to the data integral (or, if
    `data_entry` is None, to the first non-empty calculation). Data is
    overlaid as black points.
    """
    if data_entry is not None and data_entry.get("edges") is not None:
        bins = data_entry["edges"]

    src_cuts_gcf = apply_src_cuts(fsi_src)
    src_cuts_mf  = apply_src_cuts(fsi_mf)

    v_gcf, w_gcf = _fsi_vals_wts(fsi_src, var_key, channel, src_cuts_gcf)
    v_mf,  w_mf  = _fsi_vals_wts(fsi_mf,  var_key, channel, src_cuts_mf)

    hists = [
        make_hist(v_gcf, w_gcf, bins) + ("gcf",),
        make_hist(v_mf,  w_mf,  bins) + ("mf_full",),
    ]

    ref_area = float(np.sum(data_entry["y"])) if data_entry is not None else None
    hists = area_normalize(hists, ref_area=ref_area)

    for c, h, key in hists:
        draw_hist(ax, c, h, color=COLORS[key], linestyle=STYLES[key],
                  label=LABELS[key], linewidth=LINEWIDTHS[key])

    if data_entry is not None:
        yerr = np.vstack([data_entry["yerr_lo"], data_entry["yerr_hi"]])
        ax.errorbar(data_entry["x"], data_entry["y"], yerr=yerr,
                    fmt="o", ms=4, color=COLORS["data"], label=LABELS["data"],
                    elinewidth=0.8, capsize=0, zorder=10)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.text(0.95, 0.95, title, transform=ax.transAxes,
                ha="right", va="top", fontsize=11)
    if xlim is not None:
        ax.set_xlim(xlim)
    if logy:
        ax.set_yscale("log")
    ax.legend(frameon=False, loc="best")


def figure8(fsi_src, fsi_mf, outdir, clas=None, clas_ev=None):
    """Fig 8: pair c.m. momentum components for 12C(e,e'pp) after SRC cuts,
    comparing FULL GCF (p_rel>kF, AV18) to FULL MF (p_rel<kF, Fermi gas).

    Layout mirrors the bottom row of Fig 6 (4 columns: pcm_x, pcm_y, pcm_z,
    |pcm|). Data for x/y/z comes from the Hists file's sub-binned
    `epp_Pcm*_i_0_0` histograms; |pcm| has no 1D hist in the file, so we
    re-bin it from the raw event file.
    """
    if fsi_mf is None:
        print("  [SKIP] fig8 — no MF sample provided")
        return
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    components = [
        # (label, sim-key, xlim, bins, clas-map-key)  — data keys are the
        # fig6_* paper-PDF extracted points (figs 6 and 8 show the same
        # measured sample).
        (r"$p_{c.m.}^x$ [GeV/c]", "pcm_x",  (-0.5, 0.5), np.linspace(-0.5, 0.5, 21), ("fig6_x", "src")),
        (r"$p_{c.m.}^y$ [GeV/c]", "pcm_y",  (-0.5, 0.5), np.linspace(-0.5, 0.5, 21), ("fig6_y", "src")),
        (r"$p_{c.m.}^z$ [GeV/c]", "pcm_z",  (-0.2, 0.8), np.linspace(-0.2, 0.8, 21), ("fig6_z", "src")),
        (r"$|p_{c.m.}|$ [GeV/c]", "pcm_mag", (0.0, 0.8), np.linspace( 0.0, 0.8, 21), ("fig6_mag", "src")),
    ]

    for col, (xlabel, comp_key, xlim, bins, clas_key) in enumerate(components):
        de = _clas(clas, clas_key)
        title = r"$^{12}$C(e,e'pp) — SRC cuts" if col == 0 else ""
        _draw_gcf_mf(axes[col], fsi_src, fsi_mf, comp_key, bins,
                     channel="epp", xlabel=xlabel, title=title,
                     data_entry=de, xlim=xlim)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig8_pcm_GCF_vs_MF.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig8_pcm_GCF_vs_MF.png")


def figure9(fsi_src, fsi_mf, outdir, clas=None, clas_ev=None):
    """Fig 9: theta(p_miss, q) for 12C(e,e'p) (left) and 12C(e,e'pp) (right),
    after SRC cuts + CLAS acceptance, comparing FULL GCF to FULL MF.

    (e,e'p) data from `ep_Pmq` TH1 (40 bins @ 2 deg), (e,e'pp) from
    `epp_Pmq` TH1 (20 bins @ 4 deg). Paper x-ranges: ~100-180 deg for
    (e,e'p), ~120-180 deg for (e,e'pp)."""
    if fsi_mf is None:
        print("  [SKIP] fig9 — no MF sample provided")
        return

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # (e,e'p) — left
    de_ep = _clas(clas, ("fig9_ep", "src"))
    bins_ep = (de_ep["edges"] if de_ep is not None and de_ep.get("edges") is not None
               else np.linspace(100, 180, 41))
    _draw_gcf_mf(axes[0], fsi_src, fsi_mf, "theta_pmiss_q_deg", bins_ep,
                 channel="ep",
                 xlabel=r"$\theta(p_{\rm Miss},\,q)$ [deg]",
                 title=r"$^{12}$C(e,e'p)",
                 data_entry=de_ep, xlim=(100, 180))

    # (e,e'pp) — right.  Prefer paper-exact Hists.
    de_epp = _clas(clas, ("fig9_epp", "src"))
    if de_epp is not None and de_epp.get("edges") is not None:
        bins_epp = de_epp["edges"]
    elif clas_ev is not None:
        bins_epp = np.linspace(100, 180, 21)
        de_epp = events_to_entry(clas_ev, "theta_pmiss_q_deg",
                                 "src_mask_epp", bins_epp)
    else:
        bins_epp = np.linspace(100, 180, 21)
    _draw_gcf_mf(axes[1], fsi_src, fsi_mf, "theta_pmiss_q_deg", bins_epp,
                 channel="epp",
                 xlabel=r"$\theta(p_{\rm Miss},\,q)$ [deg]",
                 title=r"$^{12}$C(e,e'pp)",
                 data_entry=de_epp, xlim=(120, 180))

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig9_theta_pmiss_q_GCF_vs_MF.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig9_theta_pmiss_q_GCF_vs_MF.png")


def figure10(fsi_src, fsi_mf, outdir, clas=None):
    """Fig 10: (e,e'pp)/(e,e'p) yield ratio vs p_miss after SRC cuts,
    comparing FULL GCF to FULL MF. Same x-range as Fig 5 (0.4-1.0 GeV/c)."""
    if fsi_mf is None:
        print("  [SKIP] fig10 — no MF sample provided")
        return

    fig, ax = plt.subplots(figsize=(6, 5))

    de = _clas(clas, ("fig5_ratio", "src"))
    if de is not None and de.get("edges") is not None:
        bins = de["edges"]
    elif de is not None and de.get("x") is not None:
        xd = de["x"]; step = xd[1] - xd[0]
        bins = np.concatenate([[xd[0] - step/2], xd + step/2])
    else:
        bins = np.linspace(0.4, 1.0, 14)
    centers = 0.5 * (bins[:-1] + bins[1:])

    def _ratio(d):
        """Weighted pmiss histograms for (e,e'p) and (e,e'pp), returns ratio.
        Flat-list selection -- identical for the GCF and MF samples."""
        src_m = apply_src_cuts(d)
        m_ep  = select_eep(d) & src_m & _acc_mask_ep(d)
        m_epp = (select_eepp(d) & src_m & _acc_mask_epp(d)
                 & (d["prec"] > SRC_CUTS["precoil_min"]))
        w_ep  = d["weight"] * _acc_w_ep(d)
        w_epp = d["weight"] * _acc_w_epp(d)
        h_ep,  _ = np.histogram(d["pmiss"][m_ep],  bins=bins, weights=w_ep [m_ep])
        h_epp, _ = np.histogram(d["pmiss"][m_epp], bins=bins, weights=w_epp[m_epp])
        return np.divide(h_epp, h_ep,
                         out=np.zeros_like(h_epp, dtype=float),
                         where=h_ep > 0)

    r_gcf = _ratio(fsi_src)
    r_mf  = _ratio(fsi_mf)

    ax.plot(centers, r_gcf, color=COLORS["gcf"], linestyle=STYLES["gcf"],
            label=LABELS["gcf"], linewidth=LINEWIDTHS["gcf"])
    ax.plot(centers, r_mf,  color=COLORS["mf_full"], linestyle=STYLES["mf_full"],
            label=LABELS["mf_full"], linewidth=LINEWIDTHS["mf_full"])

    if de is not None:
        yerr = np.vstack([de["yerr_lo"], de["yerr_hi"]])
        ax.errorbar(de["x"], de["y"], yerr=yerr,
                    fmt="o", ms=4, color=COLORS["data"], label=LABELS["data"],
                    elinewidth=0.8, capsize=0, zorder=10)

    ax.set_xlabel(r"$p_{\rm Miss}$ [GeV/c]")
    ax.set_ylabel(r"A(e,e$'$pp) / A(e,e$'$p)")
    ax.text(0.95, 0.95, r"$^{12}$C — SRC cuts", transform=ax.transAxes,
            ha="right", va="top", fontsize=11)
    ax.set_xlim(0.4, 1.0)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig10_ratio_GCF_vs_MF.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig10_ratio_GCF_vs_MF.png")


# ===================================================================
#  Effective transparency summary (text result from paper)
# ===================================================================

def figure_q2_stages(pwia, fsi, outdir):
    """Q^2 sanity check: FSI cannot touch Q^2 before cuts.

    Q^2 is fixed entirely by the scattered electron, and the generator's FSI
    never modifies the primary event weight (QEGeneratorFSI.cc: the cascade
    output is a particle list, the weight is untouched). So the PWIA and FULL
    Q^2 spectra must coincide bin for bin BEFORE any selection -- the two
    samples are independent runs of the same primary configuration, matched
    by generation attempts.

    They do separate after the acceptance and after the SRC cuts, because FSI
    changes WHICH events survive those cuts. That separation is the
    effective-transparency effect, not a Q^2 effect.

    If the pre-cut curves fail to overlap, the fault is in the generation
    configuration, not in the FSI -- which is what makes this check useful.
    """
    bins = np.linspace(1.0, 5.0, 41)
    ctr = 0.5 * (bins[1:] + bins[:-1])

    def lumi(d):
        meta = d.get("_meta") or {}
        try:
            return float(meta.get("n_events_attempts"))
        except (TypeError, ValueError):
            return None

    ap, af = lumi(pwia), lumi(fsi)
    scale = (ap / af) if (ap and af) else 1.0

    stages = [
        ("as generated", lambda d: _primary(d),
         lambda d: np.ones(len(d["weight"]))),
        ("in CLAS acceptance", lambda d: select_eep(d) & _acc_mask_ep(d),
         _acc_w_ep),
        ("+ SRC cuts", lambda d: select_eep(d) & _acc_mask_ep(d) & apply_src_cuts(d),
         _acc_w_ep),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (title, mask_fn, wfn) in zip(axes, stages):
        for d, lbl, sc, kw in ((pwia, "PWIA", 1.0, dict(ls=":", color="0.35")),
                               (fsi, "FULL", scale, dict(ls="-", color="C0", lw=2))):
            m = mask_fn(d)
            wbase = (d.get("weight_gen", d["weight"]) if title == "as generated"
                     else d["weight"])
            _, h = make_hist(d["Q2"][m], (wbase * wfn(d))[m] * sc, bins)
            draw_hist(ax, ctr, h, label=lbl, **kw)
        ax.set_xlabel(r"$Q^2$  [GeV$^2$]")
        ax.set_ylabel("Counts (equal luminosity)")
        ax.set_title(title)
        ax.legend(frameon=False)

    # quantify the pre-cut agreement
    pp_, pf_ = _primary(pwia), _primary(fsi)
    wgp, wgf = pwia.get("weight_gen", pwia["weight"]), fsi.get("weight_gen", fsi["weight"])
    _, hp = make_hist(pwia["Q2"][pp_], wgp[pp_], bins)
    _, hf = make_hist(fsi["Q2"][pf_], wgf[pf_] * scale, bins)
    ok = hp > 0
    dev = np.abs(hf[ok] - hp[ok]) / hp[ok]
    axes[0].text(0.97, 0.90,
                 f"max |FULL-PWIA|/PWIA = {100*dev.max():.1f}%\n"
                 f"mean = {100*dev.mean():.1f}%",
                 transform=axes[0].transAxes, ha="right", va="top", fontsize=9)
    print(f"  Q2 pre-cut PWIA vs FULL: mean |diff| {100*dev.mean():.2f}%, "
          f"max {100*dev.max():.2f}%  (must be statistical only)")

    fig.tight_layout()
    out = os.path.join(outdir, "figD_Q2_stages.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved {os.path.basename(out)}")


# Provenance classes for the candidate-origin diagnostic, built from
# part_generation and part_parentRole.
#
#   rescatterCode == NO_INTERACTION_FATE  the input nucleon this particle
#       descends from never interacted, so the particle IS that nucleon and
#       left with exactly its generated momentum.
#   anything else  the input nucleon interacted, so this particle is a cascade
#       product -- either that nucleon after the scatter or a knocked-out
#       medium nucleon.
#
# part_generation cannot make this split: Intranuke re-emits even
# non-interacting hadrons as daughters, so generation is >=1 for everything.
# parentRole cannot either: a knocked-out medium nucleon inherits the role of
# whichever incident nucleon knocked it out.
#
# The fate integer is checked empirically rather than against the GENIE enum
# (two headers define the same names with different values -- see the note in
# GenieFSIHelpers.cc). On a 2M hA sample, fate 1 was the only value for which
# |p - p_pre| < 1 MeV, and it held for 100% of them; every other fate had a
# median displacement of 0.36-2.2 GeV/c.
NO_INTERACTION_FATE = 1

_ORIGIN_LABELS = ["struck, untouched", "product of struck",
                  "partner, untouched", "product of partner",
                  "unattributed"]


def _origin_class(role, fate):
    """Map (parentRole, rescatterCode) onto _ORIGIN_LABELS indices."""
    cls = np.full(len(role), 4, dtype=np.int64)      # default: unattributed
    for r, base in ((0, 0), (1, 2)):
        m = role == r
        cls[m & (fate == NO_INTERACTION_FATE)] = base
        cls[m & (fate != NO_INTERACTION_FATE)] = base + 1
    return cls


def figure_origins(fsi_src, fsi_mf, outdir):
    """Where do the detected lead and recoil actually come from?

    The paper's mean-field argument (Fig 8) is that FULL MF has a much
    BROADER pair c.m. distribution than the GCF. In the MF sample the second
    proton can only be an FSI product, so p_cm = p_miss + p_recoil is set by
    the momentum of whatever the cascade knocked out. This figure breaks the
    selected (e,e'pp) yield down by the provenance of the two detected
    protons and plots |p_cm| separately for each class, which is what tells
    us whether a too-narrow MF p_cm comes from the event definition, from the
    cascade medium, or from the acceptance.
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    bins = np.linspace(0.0, 1.0, 26)
    ctr = 0.5 * (bins[1:] + bins[:-1])

    for col, (d, name) in enumerate(((fsi_src, "FULL GCF"), (fsi_mf, "FULL MF"))):
        if d is None:
            continue
        m = select_eepp(d) & _acc_mask_epp(d) & apply_src_cuts(d) \
            & (d["prec"] > SRC_CUTS["precoil_min"])
        w = (d["weight"] * _acc_w_epp(d))[m]
        lc = _origin_class(d["lead_role"][m], d["lead_fate"][m])
        rc = _origin_class(d["rec_role"][m], d["rec_fate"][m])

        ax = axes[0][col]
        tot = w.sum()
        frac = [[100 * w[(lc == i) & (rc == j)].sum() / tot if tot > 0 else 0
                 for j in range(5)] for i in range(5)]
        im = ax.imshow(frac, cmap="Blues", vmin=0)
        ax.set_xticks(range(5)); ax.set_yticks(range(5))
        ax.set_xticklabels(_ORIGIN_LABELS, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(_ORIGIN_LABELS, fontsize=7)
        ax.set_xlabel("recoil origin"); ax.set_ylabel("lead origin")
        ax.set_title(f"{name} — (e,e'pp) yield share [%]")
        for i in range(5):
            for j in range(5):
                if frac[i][j] >= 0.5:
                    ax.text(j, i, f"{frac[i][j]:.0f}", ha="center", va="center",
                            fontsize=7,
                            color="white" if frac[i][j] > 50 else "black")
        fig.colorbar(im, ax=ax, fraction=0.046)

        ax = axes[1][col]
        pcm = d["pcm_mag"][m]
        draw_hist(ax, ctr, make_hist(pcm, w, bins)[1], color="k", lw=2, label="all")
        for j, lbl in enumerate(_ORIGIN_LABELS):
            sel = rc == j
            if w[sel].sum() <= 0:
                continue
            draw_hist(ax, ctr, make_hist(pcm[sel], w[sel], bins)[1],
                      label=f"recoil: {lbl}")
        ax.set_xlabel(r"$|p_{c.m.}|$  [GeV/c]")
        ax.set_ylabel("Counts")
        ax.set_title(f"{name} — by recoil origin")
        ax.legend(frameon=False, fontsize=7)

    fig.tight_layout()
    out = os.path.join(outdir, "figD_origins.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved {os.path.basename(out)}")


def figure_q2_diag(pwia, fsi, outdir, clas_ev=None):
    """Diagnostic: post-cut Q^2 spectrum, sim vs data.

    The paper's "Q^2 >~ 1.7 GeV^2" is DESCRIPTIVE, not a cut (the measured
    events reach down to Q^2 = 1.27 GeV^2 and 28% sit below 1.7): this
    figure checks that the simulation reproduces that population after the
    (e,e'pp) SRC cuts + acceptance.
    """
    fig, ax = plt.subplots(figsize=(6, 5))
    bins = np.linspace(1.0, 5.0, 21)

    def _sim(d, label_key):
        m = (select_eepp(d) & apply_src_cuts(d) & _acc_mask_epp(d)
             & (d["prec"] > SRC_CUTS["precoil_min"]))
        v = d["Q2_calc"][m]
        w = (d["weight"] * _acc_w_epp(d))[m]
        frac = w[v < 1.7].sum() / w.sum() if w.sum() > 0 else 0.0
        print(f"  Q2 diag {LABELS[label_key]:>5s}: frac(Q2<1.7) = {frac:.3f}")
        return make_hist(v, w, bins) + (label_key,)

    hists = [_sim(pwia, "pwia"), _sim(fsi, "full")]
    data_entry = None
    if clas_ev is not None:
        mask = clas_ev["src_mask_epp"]
        v = clas_ev["Q2"][mask]
        print(f"  Q2 diag  Data: frac(Q2<1.7) = {(v < 1.7).mean():.3f}")
        data_entry = events_to_entry(clas_ev, "Q2", "src_mask_epp", bins)
    ref_area = float(np.sum(data_entry["y"])) if data_entry is not None else None
    hists = area_normalize(hists, ref_area=ref_area)
    for c, h, key in hists:
        draw_hist(ax, c, h, color=COLORS[key], linestyle=STYLES[key],
                  label=LABELS[key], linewidth=LINEWIDTHS[key])
    if data_entry is not None:
        yerr = np.vstack([data_entry["yerr_lo"], data_entry["yerr_hi"]])
        ax.errorbar(data_entry["x"], data_entry["y"], yerr=yerr, fmt="o",
                    ms=4, color=COLORS["data"], label=LABELS["data"])
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_ylabel("Counts")
    ax.set_title(r"$^{12}$C(e,e'pp) — SRC cuts (diagnostic)")
    ax.axvline(1.7, color="gray", linestyle=":", linewidth=1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figD_Q2_diagnostic.png"), dpi=200)
    plt.close(fig)
    print("  Saved figD_Q2_diagnostic.png")


def print_effective_transparencies(pwia, fsi):
    """Compute and print effective transparencies (paper page 4-5).

    Paper definition: ratio of events passing the SRC selection cuts
    **within the CLAS acceptance** for FSI (or T+SCX) to PWIA, assuming
    "equivalent experimental luminosities".  Applying the same fiducial
    mask + map weight to both num and denom is essential for an apples-
    to-apples comparison.

    The two generator outputs in this project have different total
    generator-level weight sums (hA <w> ~ 100x PWIA <w>).  To recover
    the paper's luminosity-matched definition we rescale the FSI
    weights by  Sum(w)_PWIA / Sum(w)_FSI  BEFORE the acceptance + SRC
    cuts are applied.  This is exactly the ratio one would obtain if
    both samples had been thrown to the same integrated luminosity.
    """
    src = apply_src_cuts
    pwia_eep = pwia_eepp = pwia
    fsi_eep = fsi_eepp = fsi

    # Equal-luminosity normalization. The PWIA and FSI samples are two
    # INDEPENDENT generator runs, so per-luminosity yields go as
    # Sum(w_pass)/N_attempts. When both runs share the same generator
    # configuration, FSI never modifies the primary weight, so
    # Sum(w_total)/N_attempts is also equal between them and the Sigma-w
    # ratio equals the attempts ratio; the attempts ratio from the
    # metadata is the robust, auditable choice.
    def _n_attempts(d):
        meta = d.get("_meta") or {}
        try:
            return float(meta.get("n_events_attempts"))
        except (TypeError, ValueError):
            return None
    att_p, att_f = _n_attempts(pwia_eep), _n_attempts(fsi_eep)
    sumw_scale = (pwia_eep.get("weight_gen", pwia_eep["weight"])[_primary(pwia_eep)].sum()
                  / fsi_eep.get("weight_gen", fsi_eep["weight"])[_primary(fsi_eep)].sum())
    if att_p and att_f:
        lumi_scale = att_p / att_f
        ratio = lumi_scale / sumw_scale
        if not (0.9 < ratio < 1.1):
            print(f"  WARNING: attempts-based luminosity scale differs from "
                  f"the Sigma-w ratio by x{ratio:.3g} -- the PWIA and FSI "
                  f"samples were NOT generated with identical configuration "
                  f"(or generator versions). Transparencies are suspect.")
    else:
        lumi_scale = sumw_scale
        print("  WARNING: n_events_attempts missing from metadata; falling "
              "back to Sigma-w luminosity matching (only valid for "
              "identically-configured samples).")

    # PWIA, inside CLAS acceptance (fiducial masks + map weights applied).
    m_eep_pwia  = select_eep(pwia_eep)  & src(pwia_eep)  & _acc_mask_ep(pwia_eep)
    m_eepp_pwia = (select_eepp(pwia_eepp) & src(pwia_eepp) & _acc_mask_epp(pwia_eepp)
                   & (pwia_eepp["prec"] > SRC_CUTS["precoil_min"]))
    w_eep_pwia  = (pwia_eep ["weight"] * _acc_w_ep (pwia_eep))[m_eep_pwia].sum()
    w_eepp_pwia = (pwia_eepp["weight"] * _acc_w_epp(pwia_eepp))[m_eepp_pwia].sum()

    # Full FSI, inside CLAS acceptance.  Scaled to same luminosity as PWIA.
    m_eep_fsi  = select_eep(fsi_eep)  & src(fsi_eep)  & _acc_mask_ep(fsi_eep)
    m_eepp_fsi = (select_eepp(fsi_eepp) & src(fsi_eepp) & _acc_mask_epp(fsi_eepp)
                  & (fsi_eepp["prec"] > SRC_CUTS["precoil_min"]))
    w_eep_fsi  = lumi_scale * (fsi_eep ["weight"] * _acc_w_ep (fsi_eep))[m_eep_fsi].sum()
    w_eepp_fsi = lumi_scale * (fsi_eepp["weight"] * _acc_w_epp(fsi_eepp))[m_eepp_fsi].sum()

    # T+SCX applied on PWIA, inside CLAS acceptance.
    tscx_eep_w  = tscx_weight_eep (pwia_eep)
    tscx_eepp_w = tscx_weight_eepp(pwia_eepp)
    m_eep_tscx  = src(pwia_eep)  & _acc_mask_ep (pwia_eep)
    m_eepp_tscx = (src(pwia_eepp) & _acc_mask_epp(pwia_eepp)
                   & (pwia_eepp["prec"] > SRC_CUTS["precoil_min"]))
    w_eep_tscx  = (pwia_eep ["weight"] * tscx_eep_w  * _acc_w_ep (pwia_eep))[m_eep_tscx].sum()
    w_eepp_tscx = (pwia_eepp["weight"] * tscx_eepp_w * _acc_w_epp(pwia_eepp))[m_eepp_tscx].sum()

    print("\n" + "="*62)
    print("  Effective Transparencies  (within CLAS acceptance, SRC cuts)")
    print("  Both samples rescaled to equivalent experimental luminosity.")
    print("="*62)
    print(f"  {'Channel':<12} {'T+SCX':>10} {'Full':>10}   (paper: T+SCX / Full)")
    print(f"  {'(e,e\'p)':<12} {w_eep_tscx/w_eep_pwia:>10.3f} "
          f"{w_eep_fsi/w_eep_pwia:>10.3f}   (paper: 0.61  / 0.58)")
    print(f"  {'(e,e\'pp)':<12} {w_eepp_tscx/w_eepp_pwia:>10.3f} "
          f"{w_eepp_fsi/w_eepp_pwia:>10.3f}   (paper: 0.83  / 0.73)")
    print("="*62 + "\n")


# ===================================================================
#  Main
# ===================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pwia", default="events/natalie/pwia_SRC.root",
                        help="PWIA SRC ROOT file (AV18 WF, p_rel>=kF; wf_mode=2)")
    parser.add_argument("--fsi", default="events/natalie/fsi_hA_SRC_rho1.root",
                        help="Full FSI SRC ROOT file (AV18 WF, p_rel>=kF; "
                             "wf_mode=2, hA, rho^1 pair vertex)")
    parser.add_argument("--fsi-mf", default="events/natalie/fsi_hA_MF_localFG.root",
                        help=("FULL MF file: single struck nucleon from a Fermi "
                              "gas, GENIE hA FSI in the A-1 residual "
                              "(wf_mode=4). Used only for Figs 8-10. Pass empty "
                              "string to disable the MF comparison."))
    parser.add_argument("--outdir",
                        default="analysis/plots/natalie_paper/hA_acceptance",
                        help="Output directory for plots")
    parser.add_argument("--ebeam", type=float, default=5.01,
                        help="Beam energy [GeV]")
    parser.add_argument("--max-events", type=float, default=None,
                        help=("Cap on events read per SRC sample. The analysis "
                              "holds all samples in memory at once, so this is "
                              "the knob that bounds peak RSS."))
    parser.add_argument("--mf-max-events", type=float, default=None,
                        help="Separate cap for the mean-field sample.")
    parser.add_argument("--convention", choices=("detector", "wright"),
                        default="detector",
                        help=("Event-building convention. 'detector' (default) "
                              "is ours: lead = highest-momentum accepted "
                              "proton, recoil = next highest, nothing ever "
                              "dropped, SRC cuts applied afterwards. 'wright' "
                              "reproduces N. Wright's genie_gcf_loop.C: the "
                              "same lead, but the absorption veto, her "
                              "multi-recoil lottery (events with two "
                              "simultaneously detected recoils are discarded), "
                              "the acceptance map as a binary filter plus a "
                              "single-candidate-only recoil weight, the truth "
                              "p_rel > 0.25 cut, and |p_N|/|q| < 0.96."))
    parser.add_argument("--acc-phi-rotation", type=int, choices=(0, 180),
                        default=180,
                        help=("Evaluate the CLAS fiducials and acceptance map "
                              "at phi + this angle. 180 (default) follows "
                              "genie_gcf_loop.C, which rotates the electron "
                              "and protons by +pi; 0 is the previous "
                              "behaviour."))
    parser.add_argument("--absorption-veto", action=argparse.BooleanOptionalAction,
                        default=True,
                        help=("detector convention: throw the event away "
                              "(weight 0) if either generated nucleon has hA "
                              "fate 4 (absorption). ON by default; "
                              "--no-absorption-veto keeps absorbed events. "
                              "Always on under --convention wright."))
    parser.add_argument("--epp-recoil", choices=("unique", "exactly2", "highest", "all"),
                        default="unique",
                        help=("detector convention: (e,e'pp) recoil rule. "
                              "'unique' (default): exactly one other accepted "
                              "proton above 0.35 GeV/c, else the event is not "
                              "(e,e'pp) (it stays (e,e'p)); softer extras "
                              "allowed. 'exactly2': exactly two accepted "
                              "protons in total, the CLAS data topology. "
                              "'highest': previous rule, second-highest "
                              "proton. 'all': every candidate above 0.35 "
                              "GeV/c gives its own (e,e'pp) entry with that "
                              "proton as the recoil; (e,e'p) counted once."))
    parser.add_argument("--clas-data",
                        default="Acceptance/SRC_eg2_Hists_C.root",
                        help=("Pre-made CLAS paper histograms (Natalie Mail 2). "
                              "Used only for (e,e'p) panels (Fig 2, Fig 4 top, "
                              "Fig 5 ratio) — the (e,e'p) sample is not in the "
                              "raw event file."))
    parser.add_argument("--clas-events",
                        default="data/clas/SRC_e2p_C_GoodRuns_coulomb.root",
                        help=("Raw CLAS (e,e'pp) event ntuple. Used to re-bin "
                              "the (e,e'pp) panels (Fig 3, Fig 4 bottom, Fig 6, "
                              "Fig 7) at any binning we like, matching the paper's "
                              "own figures exactly.  Pass empty string to disable."))
    parser.add_argument("--acc-map",
                        default="Acceptance/map_eg2_adin.root",
                        help="Weighted CLAS acceptance map (protons+electrons)")
    parser.add_argument("--no-acceptance", action="store_true",
                        help="Skip CLAS acceptance (smearing + fiducial + map)")
    parser.add_argument("--seed", type=int, default=None,
                        help=("RNG seed for momentum smearing. "
                              "If omitted, a random seed is drawn from the "
                              "OS entropy pool; the concrete value is "
                              "printed so a specific run can be reproduced "
                              "by passing it back via --seed."))
    parser.add_argument("--map-smoothing", type=int, default=1,
                        help=("Smoothing radius for the acceptance-map "
                              "lookup: 0 = 1^3 bit-faithful AccMap.cpp; "
                              "1 = 3^3 (default); 2 = 5^3 (heavier "
                              "smoothing for noisier maps)."))
    parser.add_argument("--weight-cap-percentile", type=float, default=None,
                        metavar="PCT",
                        help=("Cap per-event generator weights at the given "
                              "percentile (e.g. 99.5), applied per sample. "
                              "Variance reduction at the cost of small bias."))

    args = parser.parse_args()

    global _ACTIVE_SRC_CUTS
    if args.convention == "wright":
        _ACTIVE_SRC_CUTS = SRC_CUTS_WRIGHT
        print("Convention: WRIGHT (genie_gcf_loop.C) — |p_N|/|q| < 0.96, "
              "absorption veto, multi-recoil events discarded, "
              "truth p_rel > 0.25")
    else:
        print(f"Convention: detector-style (default) - absorption veto "
              f"{'ON' if args.absorption_veto else 'off'}, (e,e'pp) recoil "
              f"rule '{args.epp_recoil}'")
    print(f"  acceptance evaluated at phi + {args.acc_phi_rotation} deg")

    os.makedirs(args.outdir, exist_ok=True)

    def _load_and_preprocess(label, path, acc=None, include_neutrons=False,
                             cap=None):
        print(f"Loading {label}: {path}")
        d = load_tree(path, max_events=cap)
        n_raw = len(d["weight"])
        keep = d["weight"] > 0
        if not keep.all():
            for k in list(d.keys()):
                d[k] = d[k][keep]
        print(f"  {n_raw} rows loaded, {int(keep.sum())} with weight > 0 kept")
        meta = load_meta(path)
        if meta is not None:
            print(f"    metadata: wf_mode={meta.get('wf_mode_name')}  "
                  f"N_saved={meta.get('n_events_saved')}  "
                  f"Sigma_w={meta.get('total_weight')}")
        d["_meta"] = meta
        if args.convention == "wright":
            wright_resolve(d, acc, args.ebeam, seed=args.seed,
                           phi_rot=(args.acc_phi_rotation == 180))
            nd = int(d["wright_drop"].sum())
            print(f"    [wright] discarded {nd} of {len(d['weight'])} events "
                  f"({100*nd/len(d['weight']):.1f}%) via absorption veto / "
                  f"no lead / multi-recoil / truth p_rel<0.25")
        else:
            flatlist_resolve(d, acc, args.ebeam,
                             include_neutrons=include_neutrons,
                             phi_rot=(args.acc_phi_rotation == 180),
                             absorption_veto=args.absorption_veto,
                             epp_recoil=args.epp_recoil)
            _p = _primary(d)
            print(f"    [detector] absorption veto {'ON' if args.absorption_veto else 'off'}:"
                  f" {int((d['abs_drop'] & _p).sum())} of {int(_p.sum())} events"
                  f" thrown away (weight 0)   (e,e'pp) ambiguous under "
                  f"'{args.epp_recoil}' (kept as (e,e'p) only): "
                  f"{int((d['epp_ambiguous'] & _p & ~d['abs_drop']).sum())}")
        prim = _primary(d)
        nl = d["n_lead_cand"][prim]
        print(f"    accepted nucleons per event:  0: {int((nl==0).sum())}"
              f"   1: {int((nl==1).sum())}"
              f"   2: {int((nl==2).sum())}"
              f"   >=3: {int((nl>=3).sum())}   (none dropped)")
        print(f"    channels: (e,e'p) {int((select_eep(d) & d['acc_mask_ep']).sum())}"
              f"   (e,e'pp) {int((select_eepp(d) & d['acc_mask_epp']).sum())} entries"
              f" ({int((select_eepp(d) & d['acc_mask_epp'] & prim).sum())} events)")
        if "n_recoil_hard" in d:
            _print_recoil_multiplicity(d)
        return d

    if not args.no_acceptance:
        if args.seed is None:
            args.seed = int(np.random.SeedSequence().entropy) & 0x7FFFFFFF
            seed_src = "random (OS entropy)"
        else:
            seed_src = "user-supplied"
        print(f"CLAS acceptance enabled  (map={args.acc_map}, "
              f"seed={args.seed}  [{seed_src}], map_smoothing_radius={args.map_smoothing})")
        acc = ClasAcceptance(args.acc_map, seed=args.seed,
                             map_smoothing_radius=args.map_smoothing)
    else:
        print("CLAS acceptance disabled (--no-acceptance)")
        acc = None

    # GCF SRC samples (AV18, p_rel>kF): the PWIA / T+SCX / FULL curves of
    # Figs 2-7 and the FULL GCF reference of Figs 8-10. No SRC/MF mixing.
    pwia = _load_and_preprocess("PWIA (SRC)", args.pwia, acc,
                                include_neutrons=True, cap=args.max_events)
    fsi  = _load_and_preprocess("FSI  (SRC)", args.fsi,  acc,
                                cap=args.max_events)

    # Single-nucleon mean-field sample (Figs 8-10 only). FULL only.
    fsi_mf = None
    if args.fsi_mf:
        fsi_mf = _load_and_preprocess("FSI  (MF, single nucleon)", args.fsi_mf,
                                      acc, cap=(args.mf_max_events
                                                or args.max_events))

    # Optional weight cap (variance reduction), applied per sample.
    if args.weight_cap_percentile is not None:
        pct = args.weight_cap_percentile
        if pct <= 0 or pct >= 100:
            parser.error("--weight-cap-percentile must be strictly in (0, 100)")
        print(f"\nWeight cap @ {pct}th percentile:")
        apply_weight_cap(pwia, pct, label="PWIA SRC")
        apply_weight_cap(fsi,  pct, label="FSI  SRC")
        if fsi_mf is not None:
            apply_weight_cap(fsi_mf, pct, label="FSI  MF")

    print(f"Loading CLAS paper hists: {args.clas_data}")
    clas = load_clas_hists(args.clas_data)
    if clas is not None:
        print(f"  Loaded {len(clas)} panel(s): {[k for k in clas.keys()]}")
    else:
        print("  WARNING: paper-hists file not found; (e,e'p) overlays disabled")

    # Paper-PDF extracted points override the coarse stored hists for the
    # figs 6/7/8 (e,e'pp) panels (paper plotted ~0.04-wide bins; the Hists
    # file only stores 0.067-wide ones) and provide the otherwise-missing
    # |p_cm| data column.  Figs 6 and 8 show the same measured points; the
    # fig6_* extraction is the clean one, so it serves both.
    extracted = load_paper_extracted()
    if extracted is not None:
        if clas is None:
            clas = {}
        for panel_key, ex_key in [
            (("fig6_x",   "src"), "fig6_x"),
            (("fig6_y",   "src"), "fig6_y"),
            (("fig6_z",   "src"), "fig6_z"),
            (("fig6_mag", "src"), "fig6_mag"),
            (("fig7_prel","src"), "fig7_prel"),
        ]:
            if ex_key in extracted:
                clas[panel_key] = extracted[ex_key]
        print(f"  Paper-PDF extracted points loaded for figs 6/7/8 "
              f"({len(extracted)} panels)")
    else:
        print("  NOTE: analysis/paper_extracted_data.json not found; "
              "figs 6/7/8 use stored-hist binning")

    if args.clas_events:
        print(f"Loading CLAS raw events: {args.clas_events}")
        clas_ev = load_clas_events(args.clas_events)
        if clas_ev is not None:
            n_tot = len(clas_ev["pN"])
            n_ep  = int(clas_ev["src_mask_ep"].sum())
            n_epp = int(clas_ev["src_mask_epp"].sum())
            print(f"  {n_tot} events; after SRC cuts: {n_ep} (e,e'p), "
                  f"{n_epp} (e,e'pp)  [paper: 363 (e,e'pp)]")
        else:
            print("  WARNING: raw-events file not found; (e,e'pp) overlays "
                  "fall back to the Hists file")
    else:
        clas_ev = None

    print_effective_transparencies(pwia, fsi)

    print(f"Generating plots in {args.outdir}/")
    figure2(pwia, fsi, args.outdir, clas=clas)
    figure3(pwia, fsi, args.outdir, clas=clas, clas_ev=clas_ev)
    figure4(pwia, fsi, args.outdir, clas=clas, clas_ev=clas_ev)
    figure5(pwia, fsi, args.outdir, clas=clas)
    figure5_components(pwia, fsi, args.outdir)
    figure6(pwia, fsi, args.outdir, clas=clas, clas_ev=clas_ev)
    figure7(pwia, fsi, args.outdir, clas=clas, clas_ev=clas_ev)
    figure_q2_diag(pwia, fsi, args.outdir, clas_ev=clas_ev)
    figure_q2_stages(pwia, fsi, args.outdir)
    figure_origins(fsi, fsi_mf, args.outdir)

    # Figs 8-10 compare FULL GCF (AV18 SRC, p_rel>kF) to FULL MF (single
    # struck nucleon, Fermi gas). Both are pure-wavefunction FULL calculations.
    figure8 (fsi, fsi_mf, args.outdir, clas=clas, clas_ev=clas_ev)
    figure9 (fsi, fsi_mf, args.outdir, clas=clas, clas_ev=clas_ev)
    figure10(fsi, fsi_mf, args.outdir, clas=clas)

    print("\nDone.")


if __name__ == "__main__":
    main()
