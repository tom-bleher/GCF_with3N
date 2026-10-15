#!/usr/bin/env python3
"""Shared helper: turn a 3N density table into the proton-odd Jacobi angle.

The tables are stored in the Jacobi set QEGenerator_3N::get_rho() permutes into
(QEGenerator_3N.cc:332-400):

    slot a = the NEUTRON ,  slots 2,3 = the two PROTONS
    k_a  = |p_n| ,  k_b = |p_p1 - p_p2|/2 ,  theta = angle(k_a, k_b)

grid: theta 61 x [0,180] deg, k_a 101 x [0,10] fm^-1, k_b 101 x [0,10] fm^-1,
ordered theta -> k_a -> k_b.

For the proton-odd decomposition we rebuild the three momenta exactly as the
generator does (QEGenerator_3N.cc:205-215, P_cm = 0):

    v1 = (0,0,k_a)              the neutron
    vb = k_b (sin th, 0, cos th)
    v2 = -v1/2 + vb             proton 1
    v3 = -v1/2 - vb             proton 2

and form the two proton-odd pairs.  The remaining (p,n) pair is NOT identical-
particle, so its relative momentum has a well defined sign; we fix it as
(proton - neutron)/2 throughout:

    choice A:  odd = v2,  b_A = (v3 - v1)/2
    choice B:  odd = v3,  b_B = (v2 - v1)/2

theta is rotation invariant, so it is a function of (k_a,k_b,theta) alone and no
Euler-angle sampling is needed.  Each point carries the Cartesian phase-space
weight rho * k_a^2 k_b^2 sin(theta); the 8pi^2 of
d^3p_a d^3p_b = 8pi^2 k_a^2 k_b^2 sin th dk_a dk_b dth is a constant and drops
out of a normalised shape.

Evaluating on the bare grid aliases: the map (th,k_a,k_b) -> theta' piles grid
points into particular output bins.  sample() therefore jitters each grid cell
nrep times and interpolates rho linearly, which integrates over the cell.
"""
import os
import numpy as np

ARR = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/array"
NT, NK = 61, 101
HBARC = 0.1973269804                 # GeV*fm
KMIN = 0.25                          # GeV/c, fkMinInternal (QEGenerator_3N.cc:27)

TABLES = [("AV18", "he3_3Ndensity_AV18.dnst"),
          ("AV18+UIX", "he3_3Ndensity_AV18UIX.dnst"),
          ("AV18+TM", "he3_3Ndensity_AV18TM.dnst"),
          ("CD-Bonn", "he3_3Ndensity_CDBonn.dnst"),
          ("CD-Bonn+TM", "he3_3Ndensity_CDBonnTM.dnst")]
COL = {"AV18": "#2b6cb0", "AV18+UIX": "#2f855a", "AV18+TM": "#b7791f",
       "CD-Bonn": "#c53030", "CD-Bonn+TM": "#6b46c1"}


def table_path(fn, pd_only=False):
    return os.path.join(ARR, fn.replace("he3_3Ndensity_", "he3_3Ndensity_pD_") if pd_only else fn)


def load(fn, pd_only=False):
    return np.loadtxt(table_path(fn, pd_only)).reshape(NT, NK, NK)


def _ang(a, b):
    na, nb = np.linalg.norm(a, axis=-1), np.linalg.norm(b, axis=-1)
    return np.degrees(np.arccos(np.clip(np.sum(a * b, -1) / np.maximum(na * nb, 1e-300), -1, 1)))


def _trilinear(a, idx):
    """linear interpolation of `a` at fractional indices idx (3 x N), clamped"""
    lo = np.floor(idx).astype(np.int64)
    lo = np.minimum(lo, np.array(a.shape)[:, None] - 2)
    lo = np.maximum(lo, 0)
    f = idx - lo
    out = np.zeros(idx.shape[1])
    for c in range(8):
        d = ((c >> 2) & 1, (c >> 1) & 1, c & 1)
        wgt = np.ones(idx.shape[1])
        for k in range(3):
            wgt *= f[k] if d[k] else (1.0 - f[k])
        out += wgt * a[lo[0] + d[0], lo[1] + d[1], lo[2] + d[2]]
    return out


def histograms(rho, bins, nsub=3, chunk=4):
    """Normalised theta distributions from a table, integrated over each grid cell.

    Evaluating on the bare grid aliases badly: the map (th,k_a,k_b) -> theta'
    piles grid points into particular output bins.  Each cell is therefore
    sub-divided at nsub^3 interior points and rho is interpolated linearly.
    The offsets are symmetric about the cell centre, so the theta -> 180-theta
    symmetry of the table (hence the A/B identity) is preserved exactly --
    random jitter would break it at the few-percent level.

    Returns {key: {"h": density, "mean": <theta>}} for key in
    thA / thB / th_native, each computed twice, suffixed "" (full table
    support) and "_sup" (the generator's support, all |p_i| > KMIN).
    """
    off = np.linspace(-0.5, 0.5, nsub + 2)[1:-1]
    o0, o1, o2 = (x.ravel() for x in np.meshgrid(off, off, off, indexing="ij"))
    keys = ["thA", "thB", "th_native"]
    acc = {k + sfx: np.zeros(len(bins) - 1) for k in keys for sfx in ("", "_sup")}

    dth, dk = 180.0 / (NT - 1), (10.0 / (NK - 1)) * HBARC
    ia, ib = np.meshgrid(np.arange(NK), np.arange(NK), indexing="ij")
    for t0 in range(0, NT, chunk):
        it = np.arange(t0, min(t0 + chunk, NT))
        I, A, B = np.meshgrid(it, ia[:, 0], ib[0], indexing="ij")
        base = np.stack([I.ravel(), A.ravel(), B.ravel()]).astype(np.float64)
        idx = np.concatenate([base + np.array([a, b, c])[:, None]
                              for a, b, c in zip(o0, o1, o2)], axis=1)
        idx[0] = np.clip(idx[0], 0, NT - 1)
        idx[1] = np.clip(idx[1], 0, NK - 1)
        idx[2] = np.clip(idx[2], 0, NK - 1)

        r = _trilinear(rho, idx)
        th, ka, kb = np.radians(idx[0] * dth), idx[1] * dk, idx[2] * dk
        z = np.zeros_like(ka)
        v1 = np.stack([z, ka], -1)                         # neutron
        vb = np.stack([kb * np.sin(th), kb * np.cos(th)], -1)
        v2, v3 = -0.5 * v1 + vb, -0.5 * v1 - vb            # the two protons
        n1, n2, n3 = (np.linalg.norm(v, axis=-1) for v in (v1, v2, v3))
        w = r * ka ** 2 * kb ** 2 * np.sin(th)
        sup = (n1 > KMIN) & (n2 > KMIN) & (n3 > KMIN)
        vals = {"thA": _ang(v2, 0.5 * (v3 - v1)), "thB": _ang(v3, 0.5 * (v2 - v1)),
                "th_native": np.degrees(th)}
        for k, v in vals.items():
            acc[k] += np.histogram(v, bins=bins, weights=w)[0]
            acc[k + "_sup"] += np.histogram(v, bins=bins, weights=np.where(sup, w, 0.0))[0]

    c = 0.5 * (bins[1:] + bins[:-1])
    out = {}
    for k, h in acc.items():
        s = max(h.sum(), 1e-300)
        out[k] = {"h": h / s / np.diff(bins), "mean": float(np.sum(c * h) / s)}
    return out


def pair_histograms(rho, bins_prel, bins_ang, nsub=3, chunk=4):
    """Two more pn-pair observables from a table, same cell integration as
    histograms():

      prel      |(p_proton - p_neutron)/2|   relative momentum of the (p,n) pair
      ang_pn    angle(p_proton, p_neutron)   opening angle of the two nucleons

    taking the odd proton to be v2, so the pair is (v3 = proton, v1 = neutron).
    Choice B (odd = v3) gives the identical distribution by pp exchange.

    Returns {key: {"h": density, "mean": <x>}} with keys prel / ang_pn, each
    also with a "_sup" variant restricted to the generator's support.
    """
    off = np.linspace(-0.5, 0.5, nsub + 2)[1:-1]
    o0, o1, o2 = (x.ravel() for x in np.meshgrid(off, off, off, indexing="ij"))
    acc = {k + sfx: np.zeros(len(b) - 1)
           for k, b in (("prel", bins_prel), ("ang_pn", bins_ang)) for sfx in ("", "_sup")}
    BN = {"prel": bins_prel, "ang_pn": bins_ang}

    dth, dk = 180.0 / (NT - 1), (10.0 / (NK - 1)) * HBARC
    ia, ib = np.meshgrid(np.arange(NK), np.arange(NK), indexing="ij")
    for t0 in range(0, NT, chunk):
        it = np.arange(t0, min(t0 + chunk, NT))
        I, A, B = np.meshgrid(it, ia[:, 0], ib[0], indexing="ij")
        base = np.stack([I.ravel(), A.ravel(), B.ravel()]).astype(np.float64)
        idx = np.concatenate([base + np.array([a, b, c])[:, None]
                              for a, b, c in zip(o0, o1, o2)], axis=1)
        idx[0] = np.clip(idx[0], 0, NT - 1)
        idx[1] = np.clip(idx[1], 0, NK - 1)
        idx[2] = np.clip(idx[2], 0, NK - 1)

        r = _trilinear(rho, idx)
        th, ka, kb = np.radians(idx[0] * dth), idx[1] * dk, idx[2] * dk
        z = np.zeros_like(ka)
        v1 = np.stack([z, ka], -1)                         # neutron
        vb = np.stack([kb * np.sin(th), kb * np.cos(th)], -1)
        v2, v3 = -0.5 * v1 + vb, -0.5 * v1 - vb            # the two protons
        n1, n2, n3 = (np.linalg.norm(v, axis=-1) for v in (v1, v2, v3))
        w = r * ka ** 2 * kb ** 2 * np.sin(th)
        sup = (n1 > KMIN) & (n2 > KMIN) & (n3 > KMIN)
        vals = {"prel": np.linalg.norm(0.5 * (v3 - v1), axis=-1),
                "ang_pn": _ang(v3, v1)}
        for k, v in vals.items():
            acc[k] += np.histogram(v, bins=BN[k], weights=w)[0]
            acc[k + "_sup"] += np.histogram(v, bins=BN[k], weights=np.where(sup, w, 0.0))[0]

    out = {}
    for k, h in acc.items():
        b = BN[k.replace("_sup", "")]
        c = 0.5 * (b[1:] + b[:-1]); s = max(h.sum(), 1e-300)
        out[k] = {"h": h / s / np.diff(b), "mean": float(np.sum(c * h) / s)}
    return out


def shape(vals, w, bins, mask=None):
    """normalised density over `bins`, plus its mean"""
    ww = w if mask is None else np.where(mask, w, 0.0)
    h, _ = np.histogram(vals, bins=bins, weights=ww)
    s = max(h.sum(), 1e-300)
    return h / s / np.diff(bins), float(np.sum(0.5 * (bins[1:] + bins[:-1]) * h) / s)
