#!/usr/bin/env python3
"""Acceptance analysis of run_eic_sim.sh runs: e + A -> e' N_lead N_rec (A-2), cf. PRC 105, 034001
(doi:10.1103/PhysRevC.105.034001).

Two steps:

  extract RUN_DIR [RUN_DIR ...] [-j N]
      Reduce each run to RUN_DIR/analysis_events.npz, one row per simulated event:
      generator kinematics (rest frame), cross-section weight, and per nucleon (leading,
      recoil) the far-forward hit-level and reco-level detection flags.

  plots --setting LABEL=RUN_DIR[,RUN_DIR...] [--setting ...] -o OUTDIR
      Figures and a summary table. Q2 slices of one beam setting are listed together;
      each run's events are normalized by that run's own attempts, so slices add.

  band --setting LABEL=RUN_DIR[,RUN_DIR...] [--setting ...] -o OUTDIR [--replicas N]
      GCF contact-uncertainty band (12C PWIA runs): the contacts are varied within their
      uncertainties and every event reweighted exactly by S'(k_rel)/S(k_rel); no new simulation.

Event weight: w = (generator weight) / (attempts of its run), in nb, so sum(w) over a
selection is its cross section; yields for 10 fb^-1 are sigma[nb] * 1e7, times a flat nucleon
transparency (applied only where a figure says so): 0.5 as in the paper, or 1 for GENIE FSI runs
(genQE_FSI, genQE_3N_FSI), which absorb nucleons event by event; plots --transparency overrides.

Detection, hit level (the paper's definition: enough hits to reconstruct a track, or a
deposit in a calorimeter's volume):
  B0 tracker >= 3 stations, off-momentum (OMD) >= 2 stations, Roman pots (RP) >= 2 stations,
  counting hits made by the nucleon itself; ZDC: energy deposited by the nucleon or its
  shower (hits attributed through calorimeter contributions to the nucleon's descendants).
Detection, reco level:
  B0: a B0 CKF track associated to the nucleon; central/B0 charged: an associated
  ReconstructedChargedParticle; ZDC and central/B0 calorimeters: an associated cluster of
  the nucleon's descendants. RP and OMD have no reconstruction for ion beams in EICrecon,
  so they stay at hit level in the reco-level flags.
Angles are measured with respect to the ion beam (rest-frame kinematics boosted head-on,
before the crossing angle), as in the paper's Figs. 3 and 4.
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

TRACKERS = {"B0TrackerHits": 3, "ForwardOffMTrackerHits": 2, "ForwardRomanPotHits": 2}
ZDC_HITS = ["HcalFarForwardZDCHits", "EcalFarForwardZDCHits"]
ZDC_CLUSTERS = ["HcalFarForwardZDCClusterAssociations", "EcalFarForwardZDCClusterAssociations"]
CENTRAL_CLUSTERS = ["LFHCALClusterAssociations", "EcalEndcapPClusterAssociations",
                    "HcalEndcapPInsertClusterAssociations", "B0ECalClusterAssociations",
                    "EcalBarrelClusterAssociations", "HcalBarrelClusterAssociations"]
STATION_GAP_MM = 30.0
P_IRF_MIN = 0.3  # GeV/c, the paper's cut against evaporation nucleons
ROLES = ("lead", "rec")


# ----------------------------------------------------------------------------- extract

def ancestors(st, pb, pe, pidx):
    """Index of the generator-level (status 1) ancestor of every MCParticle, -1 if none."""
    n = len(st)
    parent = np.full(n, -1)
    has = pe > pb
    parent[has] = pidx[pb[has]]
    anc = np.arange(n)
    for _ in range(64):
        step = (st[anc] != 1) & (parent[anc] >= 0)
        if not step.any():
            break
        anc[step] = parent[anc[step]]
    anc[st[anc] != 1] = -1
    return anc


def stations(z):
    if len(z) == 0:
        return 0
    z = np.sort(z)
    return 1 + int(np.sum(np.diff(z) > STATION_GAP_MM))


def final_nucleons(gen_path, nums):
    """Per event: {role: (index among the primaries after the electron, final PDG)}, (-1, 0) when absent.

    genQE trees: the converter writes e', lead, recoil, [remnant]. genDIS: e', recoil, [remnant] (the
    struck system is not written). Particle-list trees (genQE_FSI, genQE_3N_FSI): the highest-momentum
    nucleon descended from the generated lead (parentRole 0), recoil (1) or second recoil (2, 3N);
    none if the cascade absorbed it.
    """
    t = uproot.open(gen_path)["genT"]
    names = set(t.keys())
    lo, hi = int(min(nums)), int(max(nums)) + 1
    if "part_pdg" not in names:
        if "pHadron" in names:
            g = t.arrays(["rec_type"], entry_start=lo, entry_stop=hi, library="np")
            return {n: {"rec": (0, int(g["rec_type"][n - lo]))} for n in nums}
        g = t.arrays(["lead_type", "rec_type"], entry_start=lo, entry_stop=hi, library="np")
        return {n: {"lead": (0, int(g["lead_type"][n - lo])), "rec": (1, int(g["rec_type"][n - lo]))} for n in nums}
    g = t.arrays(["part_pdg", "part_parentRole", "part_px", "part_py", "part_pz"], entry_start=lo, entry_stop=hi, library="np")
    three = "N3_type" in names
    out = {}
    for n in nums:
        i = n - lo
        pdg, role = np.asarray(g["part_pdg"][i]), np.asarray(g["part_parentRole"][i])
        p = np.sqrt(np.asarray(g["part_px"][i]) ** 2 + np.asarray(g["part_py"][i]) ** 2 + np.asarray(g["part_pz"][i]) ** 2)
        res = {}
        for r, name in ((0, "lead"), (1, "rec"), (2, "rec2"))[: 3 if three else 2]:
            m = np.flatnonzero((role == r) & ((pdg == 2212) | (pdg == 2112)))
            j = int(m[np.argmax(p[m])]) if len(m) else -1
            res[name] = (j, int(pdg[j]) if j >= 0 else 0)
        out[n] = res
    return out


def extract_file(sim_path):
    """Per-event detection flags of the leading and recoil nucleon for one sim/reco chunk."""
    sim_path = Path(sim_path)
    reco_path = sim_path.with_name(sim_path.name.replace("sim_", "reco_").replace(".edm4hep.", ".edm4eic."))
    mc = ["MCParticles/MCParticles." + b for b in ("generatorStatus", "PDG", "parents_begin", "parents_end")]
    mc.append("_MCParticles_parents/_MCParticles_parents.index")
    br = ["EventHeader/EventHeader.eventNumber"] + mc
    for c in TRACKERS:
        br += [f"{c}/{c}.position.z", f"_{c}_particle/_{c}_particle.index"]
    for c in ZDC_HITS:
        br += [f"{c}Contributions/{c}Contributions.energy", f"_{c}Contributions_particle/_{c}Contributions_particle.index"]
    evap = (sim_path.parent / "gen_deex.hepmc3.json").exists()
    if evap:  # Geant4 de-excitation products of the A-2 remnant (run_eic_sim.sh --deexcite)
        br += [f"MCParticles/MCParticles.momentum.{c}" for c in "xyz"]
    s = uproot.open(sim_path)["events"].arrays(br, library="ak")
    roles = final_nucleons(sim_path.parent / "gen.root",
                           [int(n[0]) for n in s["EventHeader/EventHeader.eventNumber"]])

    r = None
    if reco_path.exists():
        rt = uproot.open(reco_path)["events"]
        keys = set(rt.keys())
        rbr = ["EventHeader/EventHeader.eventNumber"] + mc
        assoc = ["B0TrackerCKFTrackAssociations", "ReconstructedChargedParticleAssociations"]
        assoc += [c for c in ZDC_CLUSTERS + CENTRAL_CLUSTERS if f"_{c}_sim/_{c}_sim.index" in keys]
        rbr += [f"_{c}_sim/_{c}_sim.index" for c in assoc]
        rbr += ["ScatteredElectronsTruth_objIdx/ScatteredElectronsTruth_objIdx.index",
                "ReconstructedElectrons_objIdx/ReconstructedElectrons_objIdx.index",
                "_ReconstructedChargedParticleAssociations_rec/_ReconstructedChargedParticleAssociations_rec.index"]
        rbr += [f"MCParticles/MCParticles.momentum.{c}" for c in "xyz"]
        rbr += [f"ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.{c}" for c in "xyz"]
        zn = "ReconstructedHcalFarForwardZDCNeutrals"
        if f"{zn}/{zn}.energy" in keys:
            rbr += [f"{zn}/{zn}.energy"] + [f"{zn}/{zn}.momentum.{c}" for c in "xyz"]
        r = rt.arrays(rbr, library="ak")
        r_index = {int(n[0]): k for k, n in enumerate(r["EventHeader/EventHeader.eventNumber"])}

    rows = []
    for k in range(len(s)):
        num = int(s["EventHeader/EventHeader.eventNumber"][k][0])
        st = ak.to_numpy(s[mc[0]][k])
        prim = np.flatnonzero(st == 1)
        anc = ancestors(st, ak.to_numpy(s[mc[2]][k]), ak.to_numpy(s[mc[3]][k]), ak.to_numpy(s[mc[4]][k]))
        row = {"num": num, "n_primaries": len(prim)}
        ev_roles = roles[num]
        pick = lambda pr_, j: pr_[1 + j] if 0 <= j and 1 + j < len(pr_) else -2  # noqa: E731
        idx = {role: pick(prim, j) for role, (j, _) in ev_roles.items()}
        for role, (_, pdg_) in ev_roles.items():
            row[f"{role}_final_pdg"] = pdg_
        for role, j in idx.items():
            for c in TRACKERS:
                z = ak.to_numpy(s[f"{c}/{c}.position.z"][k])[ak.to_numpy(s[f"_{c}_particle/_{c}_particle.index"][k]) == j]
                row[f"{role}_{c}"] = stations(z)
            e = 0.0
            for c in ZDC_HITS:
                ce = ak.to_numpy(s[f"{c}Contributions/{c}Contributions.energy"][k])
                cp = ak.to_numpy(s[f"_{c}Contributions_particle/_{c}Contributions_particle.index"][k])
                e += float(ce[anc[cp] == j].sum()) if len(cp) else 0.0
            row[f"{role}_zdc_edep"] = e
        if evap:
            # primaries after e', lead and recoil are the de-excitation products (appended by tools/deexcite)
            n_known = 1 + max((j for j, _ in ev_roles.values() if j >= 0), default=-1) + 1
            pdg_all = ak.to_numpy(s["MCParticles/MCParticles.PDG"][k])
            mom = [ak.to_numpy(s[f"MCParticles/MCParticles.momentum.{c}"][k]) for c in "xyz"]
            ev = []
            for j in prim[n_known:]:
                hits = []
                for c in TRACKERS:
                    z = ak.to_numpy(s[f"{c}/{c}.position.z"][k])[ak.to_numpy(s[f"_{c}_particle/_{c}_particle.index"][k]) == j]
                    hits.append(stations(z))
                e = 0.0
                for c in ZDC_HITS:
                    ce = ak.to_numpy(s[f"{c}Contributions/{c}Contributions.energy"][k])
                    cp = ak.to_numpy(s[f"_{c}Contributions_particle/_{c}Contributions_particle.index"][k])
                    e += float(ce[anc[cp] == j].sum()) if len(cp) else 0.0
                ev.append((int(pdg_all[j]), float(mom[0][j]), float(mom[1][j]), float(mom[2][j]), *hits, e))
            row["evap"] = ev
        if r is not None and num in r_index:
            q = r_index[num]
            rst = ak.to_numpy(r[mc[0]][q])
            ranc = ancestors(rst, ak.to_numpy(r[mc[2]][q]), ak.to_numpy(r[mc[3]][q]), ak.to_numpy(r[mc[4]][q]))
            rprim = np.flatnonzero(rst == 1)
            ridx = {role: pick(rprim, j) for role, (j, _) in ev_roles.items()} if len(rprim) == len(prim) else None
            row["reco"] = ridx is not None
            row["e_truth_matched"] = len(r["ScatteredElectronsTruth_objIdx/ScatteredElectronsTruth_objIdx.index"][q]) > 0
            row["e_reco"] = len(r["ReconstructedElectrons_objIdx/ReconstructedElectrons_objIdx.index"][q]) > 0
            for role in ev_roles:
                j = ridx[role] if ridx else -2
                sim_of = lambda c: ak.to_numpy(r[f"_{c}_sim/_{c}_sim.index"][q])  # noqa: E731
                row[f"{role}_b0track"] = bool(np.any(sim_of("B0TrackerCKFTrackAssociations") == j))
                row[f"{role}_charged"] = bool(np.any(sim_of("ReconstructedChargedParticleAssociations") == j))
                for group, cols in (("zdc_cluster", ZDC_CLUSTERS), ("central_cluster", CENTRAL_CLUSTERS)):
                    hit = False
                    for c in cols:
                        if f"_{c}_sim/_{c}_sim.index" in r.fields:
                            m = sim_of(c)
                            hit |= bool(np.any(ranc[m[m >= 0]] == j)) if len(m) else False
                    row[f"{role}_{group}"] = hit
                # true lab momentum (after the afterburner) and the reconstructed one: the linked charged
                # particle (B0 or central track), else the ZDC neutron for a ZDC-detected neutron
                if j >= 0:
                    for c in "xyz":
                        row[f"{role}_true_p{c}"] = float(r[f"MCParticles/MCParticles.momentum.{c}"][q][j])
                    cm = np.flatnonzero(sim_of("ReconstructedChargedParticleAssociations") == j)
                    rec = None
                    if len(cm):
                        ri = int(r["_ReconstructedChargedParticleAssociations_rec/_ReconstructedChargedParticleAssociations_rec.index"][q][cm[0]])
                        rec = [float(r[f"ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.{c}"][q][ri]) for c in "xyz"]
                    elif row[f"{role}_zdc_cluster"] and row[f"{role}_final_pdg"] == 2112 and f"{zn}/{zn}.energy" in r.fields:
                        en = ak.to_numpy(r[f"{zn}/{zn}.energy"][q])
                        if len(en):
                            ni = int(np.argmax(en))
                            rec = [float(r[f"{zn}/{zn}.momentum.{c}"][q][ni]) for c in "xyz"]
                    if rec is not None:
                        for c, v in zip("xyz", rec):
                            row[f"{role}_reco_p{c}"] = v
        else:
            row["reco"] = False
        rows.append(row)
    return rows


def lab_to_irf(px, py, pz, a, p_nucleon, m_ion, crossing=0.025, mass=0.93892):
    """|p| in the ion rest frame of lab momenta, undoing the nominal crossing angle (ion beam along
    (-sin, 0, cos) of the crossing angle) and the nominal ion boost, as an analysis would."""
    ca, sa = np.cos(crossing), np.sin(crossing)
    x, z = px * ca + pz * sa, -px * sa + pz * ca
    bg = a * p_nucleon / m_ion
    g = np.hypot(1.0, bg)
    e = np.sqrt(mass**2 + px**2 + py**2 + pz**2)
    zr = g * z - bg * e
    return np.sqrt(x**2 + py**2 + zr**2)


def cc1_2021_ratio(e_beam, k, p, is_p, m_n=0.93892):
    """Per event sigma_cc1(as generated for PRC 105, 034001) / sigma_cc1(current), rest-frame vectors.

    De Forest's cc1 as in eNCrossSection::sigmaccn (n = 1, Kelly form factors, mu_p = 2.79, mu_n = -1.91 as
    in constants.hh). The paper's events were generated with the longitudinal factor Q^4/|q|^2, which
    1a92ca1 (#33) fixed to (Q^2/|q|^2)^2; everything else in the weight is unchanged, so this ratio
    reweights a sample to the paper's cross section. It matches the C++, current and with the old line,
    to 1e-10 (#188). k: scattered electron, p: final leading nucleon (as QEGenerator passes vLead).
    """
    q = np.column_stack([-k[:, 0], -k[:, 1], e_beam - k[:, 2]])
    q2v = (q**2).sum(1)
    q2 = q2v - (e_beam - np.linalg.norm(k, axis=1)) ** 2
    e = np.sqrt(m_n**2 + (p**2).sum(1))
    eb = np.sqrt(m_n**2 + ((p - q) ** 2).sum(1))
    q2b = q2v - (e - eb) ** 2
    tau = q2 / (4 * m_n**2)
    kelly = lambda a1, b1, b2, b3: (1 + a1 * tau) / (1 + b1 * tau + b2 * tau**2 + b3 * tau**3)  # noqa: E731
    ge = np.where(is_p, kelly(-0.24, 10.98, 12.82, 21.97), 1.70 * tau / (1 + 3.3 * tau) / (1 + q2 / 0.71) ** 2)
    gm = np.where(is_p, 2.79 * kelly(0.12, 10.97, 18.86, 6.55), -1.91 * kelly(2.33, 14.72, 24.20, 84.1))
    f1, kf2 = (ge + gm * tau) / (1 + tau), (gm - ge) / (1 + tau)
    x = f1**2 + q2b / (4 * m_n**2) * kf2**2
    w_c = ((e + eb) ** 2 * x - q2v * (f1 + kf2) ** 2) / (4 * e * eb)
    w_t = q2b * (f1 + kf2) ** 2 / (2 * eb * e)
    p_perp = np.linalg.norm(np.cross(p, q), axis=1) / np.sqrt(q2v)
    w_s = p_perp**2 * x / (e * eb)
    w_i = -p_perp * (eb + e) * x / (e * eb)
    t2 = np.tan(0.5 * np.arctan2(np.hypot(k[:, 0], k[:, 1]), k[:, 2])) ** 2
    a_, b_ = np.cross(q, k), np.cross(q, p)
    cphi = np.clip((a_ * b_).sum(1) / (np.linalg.norm(a_, axis=1) * np.linalg.norm(b_, axis=1)), -1, 1)
    rest = ((q2 / (2 * q2v) + t2) * w_t + q2 / q2v * np.sqrt(q2 / q2v + t2) * w_i * cphi
            + (q2 / q2v * cphi**2 + t2) * w_s)
    return (q2**2 / q2v * w_c + rest) / ((q2 / q2v) ** 2 * w_c + rest)


def boost_angles(p_irf, a, p_nucleon, m_ion):
    """Lab polar angle w.r.t. the ion beam [mrad] and lab momentum of rest-frame 3-momenta."""
    bg = a * p_nucleon / m_ion
    g = np.hypot(1.0, bg)
    px, py, pz = -p_irf[:, 0], p_irf[:, 1], -p_irf[:, 2]
    e = np.sqrt(0.93892**2 + (p_irf**2).sum(1))
    pz_lab = g * pz + bg * e
    return np.arctan2(np.hypot(px, py), pz_lab) * 1e3, np.sqrt(px**2 + py**2 + pz_lab**2)


def extract(run, jobs):
    run = Path(run)
    conv = json.loads((run / "gen.hepmc3.json").read_text())
    beams = conv["beams"]
    sims = sorted(run.glob("sim_*.edm4hep.root"))
    if not sims:
        sys.exit(f"{run}: no sim files")
    with ProcessPoolExecutor(jobs) as ex:
        rows = [r for part in ex.map(extract_file, sims) for r in part]
    num = np.array([r["num"] for r in rows])
    e_irf = beams["E_IRF_run_GeV"]
    gt = uproot.open(run / "gen.root")["genT"]
    names = set(gt.keys())
    # generator-level nucleons per role, by run type: 2N SRC (genQE, genQE_FSI incl. mean field),
    # 3N SRC (genQE_3N_FSI), tagged DIS off an SRC pair (genDIS, no final leading nucleon)
    if "p1_gen" in names:
        run_kind = "3N"
        g = gt.arrays(["pe", "p1_gen", "p2_gen", "p3_gen", "N1_type", "N2_type", "N3_type", "weight"], library="np")
        gen_roles = {"lead": ("p1_gen", "N1_type"), "rec": ("p2_gen", "N2_type"), "rec2": ("p3_gen", "N3_type")}
        g["q"] = np.column_stack([-g["pe"][:, 0], -g["pe"][:, 1], e_irf - g["pe"][:, 2]])
    elif "pHadron" in names:
        run_kind = "DIS"
        g = gt.arrays(["pe", "q", "pRec", "rec_type", "lead_type", "p1_onshell", "E1_onshell", "weight"], library="np")
        gen_roles = {"rec": ("pRec", "rec_type")}
    else:
        run_kind = "2N"
        g = gt.arrays(["pe", "pLead", "pRec", "q", "lead_type", "rec_type", "weight"], library="np")
        gen_roles = {"lead": ("pLead", "lead_type"), "rec": ("pRec", "rec_type")}
    pe = g["pe"][num]
    ee = np.linalg.norm(pe, axis=1)
    q2 = 2 * e_irf * (ee - pe[:, 2])
    nu = e_irf - ee
    xb = q2 / (2 * 0.938272 * nu)
    out = {"num": num, "w_nb": g["weight"][num] / conv["attempts"], "q2": q2, "xb": xb}
    if run_kind == "DIS":
        out["lead_type"] = g["lead_type"][num]  # the struck nucleon (no final-state leading nucleon)
    if run_kind == "2N":
        out["cc1_2021_ratio"] = cc1_2021_ratio(e_irf, pe, g["pLead"][num], g["lead_type"][num] == 2212)
    z, a, m = {"C12": (6, 12, 11.174863), "d": (1, 2, 1.875613), "He4": (2, 4, 3.727379)}[beams["ion"]]
    for role, (pb, tb) in gen_roles.items():
        p = g[pb][num]
        out[f"{role}_type"] = g[tb][num]
        out[f"{role}_p_irf"] = np.linalg.norm(p, axis=1)
        out[f"{role}_theta_mrad"], out[f"{role}_p_lab"] = boost_angles(p, a, beams["ion_momentum_per_nucleon_GeV"], m)
    qv = g["q"][num]
    if "lead" in gen_roles:
        # missing momentum and energy of the struck nucleon (PRC 105, 034001 Fig. 7 cuts):
        # p_miss = |p_lead - q|, E_miss = omega - T_lead - T_(A-1), with T_(A-1) = recoil energy of the A-1 system
        pl = g[gen_roles["lead"][0]][num]
        out["pmiss"] = np.linalg.norm(pl - qv, axis=1)
        m_am1 = m - 0.93892
        t_lead = np.sqrt(0.93892**2 + (pl**2).sum(1)) - 0.93892
        out["emiss"] = nu - t_lead - (np.sqrt(m_am1**2 + out["pmiss"] ** 2) - m_am1)
    if run_kind == "DIS":
        # tagged-DIS variables: x' = Q2 / (2 p1.q) with the struck nucleon's (on-shell) four-momentum,
        # and the spectator light-cone fraction alpha_s = (E_s - p_s.q^) / m_N (1 at rest)
        p1, e1 = g["p1_onshell"][num], g["E1_onshell"][num]
        out["xprime"] = q2 / (2 * (e1 * nu - (p1 * qv).sum(1)))
        ps = g["pRec"][num]
        es = np.sqrt(0.93892**2 + (ps**2).sum(1))
        qhat = qv / np.linalg.norm(qv, axis=1)[:, None]
        out["alpha_s"] = (es - (ps * qhat).sum(1)) / 0.93892
    if run_kind == "3N":
        out["rec2_p_irf"] = np.linalg.norm(g["p3_gen"][num], axis=1)
    if any("evap" in r for r in rows):
        # flat per-product arrays, evap_event = row of the event
        ev = [(i, *t) for i, r in enumerate(rows) for t in r.get("evap", [])]
        cols = ["event", "pdg", "px", "py", "pz", "B0TrackerHits", "ForwardOffMTrackerHits", "ForwardRomanPotHits", "zdc_edep"]
        arr = np.array(ev, dtype=float).reshape(-1, len(cols))
        for j, c in enumerate(cols):
            out[f"evap_{c}"] = arr[:, j].astype(int) if c in ("event", "pdg") or c.endswith("Hits") else arr[:, j]
        px, py, pz = out["evap_px"], out["evap_py"], out["evap_pz"]
        ca, sa = np.cos(0.025), np.sin(0.025)  # angle to the ion beam (nominal crossing angle undone)
        out["evap_theta_mrad"] = 1e3 * np.arctan2(np.hypot(px * ca + pz * sa, py), -px * sa + pz * ca)
        nucleon = np.isin(out["evap_pdg"], (2212, 2112))
        out["evap_p_irf"] = np.where(nucleon, lab_to_irf(px, py, pz, a, beams["ion_momentum_per_nucleon_GeV"], m), np.nan)
        for r in rows:
            r.pop("evap", None)
    keys = sorted({k for r in rows for k in r if k != "num"})
    for k in keys:
        mom = k.endswith(("_px", "_py", "_pz"))
        out[k] = np.array([r.get(k, np.nan if mom else False) for r in rows])
    for role in ("lead", "rec", "rec2"):
        for kind in ("true", "reco"):
            if f"{role}_{kind}_px" in out:
                v = [out[f"{role}_{kind}_p{c}"] for c in "xyz"]
                out[f"{role}_{kind}_p_lab_mag"] = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
                out[f"{role}_{kind}_p_irf"] = lab_to_irf(*v, a, beams["ion_momentum_per_nucleon_GeV"], m)
    meta = {"run": str(run), "kind": run_kind, "beams": beams, "attempts": conv["attempts"], "cuts": conv["cuts"],
            "events": len(num), "q2_range": json.loads((run / "run.json").read_text()).get("q2_range_GeV2")
            if (run / "run.json").exists() else None}
    np.savez_compressed(run / "analysis_events.npz", meta=json.dumps(meta), **out)
    need = {"2N": 3, "3N": 4, "DIS": 2}[run_kind]  # e' + nucleons (mean field: e' + one nucleon)
    bad = int(np.sum((out["n_primaries"] < need) & (out.get("rec_type", np.ones(len(num))) != 0)))
    print(f"{run.name}: {len(num)} events -> analysis_events.npz (reco matched {out['reco'].mean():.3f}"
          f"{f', {bad} events with < 3 primaries' if bad else ''})")


# ----------------------------------------------------------------------------- selections

def detected(d, role, level):
    """Far-forward detection of a nucleon (role: lead|rec) at 'hit' or 'reco' level."""
    # the final-state particle decides the detector (FSI can change it, e.g. charge exchange)
    is_p = d.get(f"{role}_final_pdg", d[f"{role}_type"]) == 2212
    if level == "hit":
        p_ok = (d[f"{role}_B0TrackerHits"] >= 3) | (d[f"{role}_ForwardOffMTrackerHits"] >= 2) | (d[f"{role}_ForwardRomanPotHits"] >= 2)
        n_ok = d[f"{role}_zdc_edep"] > 0
    else:
        p_ok = d[f"{role}_b0track"] | (d[f"{role}_ForwardOffMTrackerHits"] >= 2) | (d[f"{role}_ForwardRomanPotHits"] >= 2)
        n_ok = d[f"{role}_zdc_cluster"]
    return np.where(is_p, p_ok, n_ok)


def detected_anywhere(d, role):
    """Far-forward (hit level) or central/B0 reconstruction (charged particle or cluster)."""
    return detected(d, role, "hit") | d[f"{role}_charged"] | d[f"{role}_central_cluster"]


def load(dirs):
    parts = [np.load(Path(r) / "analysis_events.npz") for r in dirs]
    keys = [k for k in parts[0].files if k != "meta"]
    stale = {r: sorted(set(keys) - set(p.files)) for r, p in zip(dirs, parts) if set(keys) - set(p.files)}
    if stale:  # Q2 slices extracted by different versions of this script
        sys.exit("plots: these runs lack columns of " + str(dirs[0]) + "; re-run extract on all runs of the setting:\n"
                 + "\n".join(f"  {r}: {', '.join(m[:5])}{' ...' if len(m) > 5 else ''}" for r, m in stale.items()))
    d = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    if "evap_event" in d:  # per-product event index: offset by the events of the preceding runs
        off = np.cumsum([0] + [len(p["num"]) for p in parts[:-1]])
        d["evap_event"] = np.concatenate([p["evap_event"] + o for p, o in zip(parts, off)])
    if "evap_pdg" in d and not len(d["evap_pdg"]):  # no A-2 remnant (e.g. e+d): no evaporation panels (#202)
        d = {k: v for k, v in d.items() if not k.startswith("evap_")}
    metas = [json.loads(str(p["meta"])) for p in parts]
    return d, metas


def has_fsi(run):
    """True for a GENIE FSI run (genQE_FSI or genQE_3N_FSI without -n), from its run.json; None when the
    run directory has no run.json, so the caller can refuse to guess."""
    rj = Path(run) / "run.json"
    if not rj.exists():
        return None
    meta = json.loads(rj.read_text())
    if "generator" not in meta:  # the first production batch predates the key; all of it is genQE (PWIA)
        print(f"warning: {rj} has no 'generator' entry; taking it as a genQE (PWIA) run", file=sys.stderr)
    gen = meta.get("generator", ["genQE"])
    return gen[0] in ("genQE_FSI", "genQE_3N_FSI") and "-n" not in gen[1:]


# ----------------------------------------------------------------------------- plots

BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"


def style(ax):
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def wstats(x, w, ok, bins):
    gen, _ = np.histogram(x, bins, weights=w)
    acc, _ = np.histogram(x[ok], bins, weights=w[ok])
    s1, _ = np.histogram(x[ok], bins, weights=w[ok] ** 2)
    s0, _ = np.histogram(x[~ok], bins, weights=w[~ok] ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = acc / gen
        err = np.sqrt(s1 * (1 - eff) ** 2 + s0 * eff**2) / gen
    return gen, acc, eff, err


def plots(settings, outdir, lumi_nb_inv, transparency=None, cc1="fixed"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = {lab: load(dirs) for lab, dirs in settings.items()}
    if cc1 == "2021":
        # the paper's cross section (cc1 before the longitudinal fix of #33), for like-for-like comparison
        for lab, (d, _) in data.items():
            if "cc1_2021_ratio" not in d:
                sys.exit(f"plots --cc1 2021: {lab} has no cc1_2021_ratio (2N runs only; re-run extract)")
            d["w_nb"] = d["w_nb"] * d["cc1_2021_ratio"]
        print("plots: weights reweighted to the cc1 cross section PRC 105, 034001 was generated with (before #33)")
    # the flat transparency stands in for FSI in PWIA runs; GENIE FSI runs already absorb nucleons
    # event by event, so they take 1 unless --transparency says otherwise
    if transparency is None:
        # from the directories given on the command line, not the paths stored in the npz at extract time
        fsi = {r: has_fsi(r) for dirs in settings.values() for r in dirs}
        unknown = [r for r, f in fsi.items() if f is None]
        if unknown:
            sys.exit("plots: no run.json in " + ", ".join(unknown) + "; cannot tell GENIE FSI (transparency 1) "
                     "from PWIA (0.5) runs; pass --transparency")
        fsi = set(fsi.values())
        if len(fsi) > 1:
            sys.exit("plots: GENIE FSI and PWIA runs need different transparencies; pass --transparency")
        transparency, why = (1.0, "GENIE FSI runs") if fsi == {True} else (0.5, "PWIA runs, as in PRC 105, 034001")
    else:
        why = "--transparency"
    print(f"plots: nucleon transparency {transparency:g} ({why})")
    scale = lumi_nb_inv * transparency
    summary = {}

    # The figures follow the style of PRC 105, 034001 (Figs. 3-6). The paper's evaporation nucleons (from
    # BeAGLE) exist only in --deexcite runs (BeAGLE/FLUKA or Geant4 de-excitation of the A-2 remnant);
    # there Fig. 3 adds them scaled by 1/3 and Fig. 4 gets their column on the left, as in the paper,
    # and figEvap below shows them in more detail.
    green, red, orange, magenta, brown, blue_ev = "#45c045", "#f0281e", "#ff9800", "#ee5fbf", "#9a5a1c", "#1f3fe0"
    low_energy = {lab: metas[0]["beams"]["ion_momentum_per_nucleon_GeV"] < 75 for lab, (_, metas) in data.items()}
    for lab in data:
        summary.setdefault(lab, {})["transparency"] = transparency
        summary[lab]["cc1"] = ("before the #33 fix, as PRC 105, 034001 was generated" if cc1 == "2021"
                               else "current (longitudinal factor fixed in #33)")

    def place_label(ax, text, color, pts, anchor, taken, bw=0.18, bh=0.13, **kw):
        """Put text at the free spot nearest anchor (axes fractions): no data point inside its box."""
        best = None
        for xf in np.linspace(0.1, 0.9, 41):
            for yf in np.linspace(0.1, 0.78, 35):
                if np.any((np.abs(pts[:, 0] - xf) < bw / 2) & (np.abs(pts[:, 1] - yf) < bh / 2)):
                    continue
                if any(abs(xf - a) < (bw + aw) / 2 and abs(yf - b) < (bh + ah) / 2 for a, b, aw, ah in taken):
                    continue
                dd = (xf - anchor[0]) ** 2 + (yf - anchor[1]) ** 2
                if best is None or dd < best[0]:
                    best = (dd, xf, yf)
        xf, yf = best[1:] if best else anchor
        taken.append((xf, yf, bw, bh))
        ax.text(xf, yf, text, transform=ax.transAxes, color=color, ha="center", va="center", **kw)

    def frame(ax, minor=True):
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=12)
        if minor:
            ax.minorticks_on()

    # Mean-field runs (single-nucleon knockout, no recoil): PRC 105, 034001 Fig. 7 only
    mf = {lab: bool(np.all(d["rec_type"] == 0)) for lab, (d, _) in data.items()}
    if any(mf.values()):
        if not all(mf.values()):
            sys.exit("plots: mean-field and SRC settings cannot be mixed in one call")
        fig, ax = plt.subplots(figsize=(6.6, 3.8))
        bins = np.arange(3, 25.5, 1.0)
        cen7 = 0.5 * (bins[1:] + bins[:-1])
        curves7 = {}
        for lab, (d, _) in data.items():
            w = d["w_nb"] * scale
            kin = (d["q2"] > 3) & (d["pmiss"] < 0.3) & (d["emiss"] < 0.08)
            ok = kin & d["e_truth_matched"] & detected_anywhere(d, "lead")
            h = np.histogram(d["q2"][ok], bins, weights=w[ok])[0]
            summary[lab].update({"mf_sigma_nb_cuts": float(d["w_nb"][kin].sum()),
                                 "mf_e_and_nucleon_detected": float(w[ok].sum() / w[kin].sum()),
                                 "mf_q2_reach_100_per_GeV2": float(bins[1:][h >= 100].max()) if (h >= 100).any() else None})
            c, mk = (red, "o") if low_energy[lab] else ("#1029ff", "s")
            g = h > 0
            ax.plot(cen7[g], h[g], ls="none", marker=mk, ms=5, color=c)
            if "e_reco" in d:  # open markers: the electron identified by EICrecon (ReconstructedElectrons)
                oke = kin & d["e_reco"] & detected_anywhere(d, "lead")
                he = np.histogram(d["q2"][oke], bins, weights=w[oke])[0]
                summary[lab].update({"mf_e_identified_and_nucleon_detected": float(w[oke].sum() / w[kin].sum()),
                                     "mf_q2_reach_100_per_GeV2_e_identified":
                                         float(bins[1:][he >= 100].max()) if (he >= 100).any() else None})
                ax.plot(cen7[he > 0], he[he > 0], ls="none", marker=mk, ms=5, mfc="none", color=c)
            curves7[lab] = (h, c)
        ax.set_yscale("log")
        ax.set_xlim(1, 25)
        order = sorted(curves7, key=lambda l: -float(np.sum(curves7[l][0])))
        for i, lab in enumerate(order):
            h, c = curves7[lab]
            k = 12 if i == 0 else 3
            ax.text(cen7[k], max(h[k], 1e-3) * (9 if i == 0 else 0.2), lab, color=c, fontsize=15, ha="center", va="center")
        if any("e_reco" in d for d, _ in data.values()):
            ax.text(0.97, 0.95, "filled: electron truth-matched\nopen: electron identified", transform=ax.transAxes,
                    fontsize=9, ha="right", va="top", color="#444444")
        frame(ax)
        ax.set_ylabel("Events per bin", fontsize=13, loc="top")
        ax.set_xlabel(r"Q$^2$ [GeV$^2$/c$^4$]", fontsize=13, loc="right")
        fig.tight_layout()
        fig.savefig(outdir / "fig7_q2_mean_field.pdf")
        plt.close(fig)
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        for lab in data:
            print(f"{lab}: {summary[lab]}")
        print(f"wrote {outdir}/fig7_q2_mean_field.pdf summary.json")
        return

    kinds = {lab: metas[0].get("kind", "2N") for lab, (_, metas) in data.items()}
    if len(set(kinds.values())) > 1:
        sys.exit(f"plots: settings of different run types cannot be mixed: {kinds}")
    run_kind = next(iter(kinds.values()))

    def acceptance_grid(var, bins, xlabel, roles, fname, cut_irf=True):
        """Acceptance (hit and reco level) vs a generator variable, rows = settings, cols = p / n."""
        fig, axes = plt.subplots(len(data), 2, figsize=(8.4, 2.8 * len(data) + 0.9), sharex=True, sharey=True,
                                 squeeze=False, gridspec_kw={"wspace": 0, "hspace": 0})
        cen_ = 0.5 * (bins[1:] + bins[:-1])
        for row, (lab, (d, _)) in enumerate(data.items()):
            for col, (code, name) in enumerate(((2212, "proton"), (2112, "neutron"))):
                ax = axes[row, col]
                for j, (level, c) in enumerate((("hit", INK), ("reco", "#1f3fe0"))):
                    xs, ws, oks = [], [], []
                    for role in roles:
                        mm = d[f"{role}_type"] == code
                        if cut_irf:
                            mm &= d[f"{role}_p_irf"] > P_IRF_MIN
                        xs.append(d[var(role)][mm]); ws.append(d["w_nb"][mm]); oks.append(detected(d, role, level)[mm])
                    x, w, ok = np.concatenate(xs), np.concatenate(ws), np.concatenate(oks)
                    if not len(x):
                        continue
                    gen, acc, eff, err = wstats(x, w, ok, bins)
                    g = gen > 0
                    ax.errorbar(cen_[g] + (j - 0.5) * 0.15 * np.diff(bins)[g], eff[g], yerr=err[g], fmt="o", ms=3.5,
                                color=c, elinewidth=1.1, label=f"{level} level")
                ax.text(0.95, 0.93, f"{name}, {lab}", transform=ax.transAxes, ha="right", va="top", fontsize=12,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5), zorder=5)
                ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11)
            axes[row, 0].set_ylabel("Acceptance", fontsize=12)
            axes[row, 0].set_ylim(0, 1.12)
        for ax in axes[-1]:
            ax.set_xlabel(xlabel, fontsize=12)
        axes[0, 1].legend(frameon=False, fontsize=10, loc="center right")
        fig.tight_layout()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(outdir / fname)
        plt.close(fig)

    if run_kind == "3N":
        # 3N SRC: how often 0, 1 or 2 of the recoils are detected far forward (hit level), by recoil-pair type,
        # for events where both recoils carry p_IRF > 0.3 GeV/c
        fig, axes = plt.subplots(1, len(data), figsize=(5.0 * len(data), 3.9), sharey=True, squeeze=False)
        pair_name = {(2212, 2212): "pp", (2212, 2112): "pn", (2112, 2212): "pn", (2112, 2112): "nn"}
        for ax, (lab, (d, _)) in zip(axes[0], data.items()):
            both = (d["rec_p_irf"] > P_IRF_MIN) & (d["rec2_p_irf"] > P_IRF_MIN)
            ndet = detected(d, "rec", "hit").astype(int) + detected(d, "rec2", "hit").astype(int)
            nreco = detected(d, "rec", "reco").astype(int) + detected(d, "rec2", "reco").astype(int)
            cats = sorted({pair_name[(int(a1), int(b1))] for a1, b1 in zip(d["rec_type"][both], d["rec2_type"][both])},
                          key=["pp", "pn", "nn"].index)
            w = d["w_nb"]
            summary[lab]["3N"] = {}
            for i, cat in enumerate(cats):
                mc_ = both & np.array([pair_name[(int(a1), int(b1))] == cat for a1, b1 in zip(d["rec_type"], d["rec2_type"])])
                fr = [float(w[mc_ & (ndet == k)].sum() / w[mc_].sum()) for k in (0, 1, 2)]
                summary[lab]["3N"][cat] = {"share": float(w[mc_].sum() / w[both].sum()), "0 recoils": fr[0], "1 recoil": fr[1],
                                          "2 recoils": fr[2], "2 recoils (reco)": float(w[mc_ & (nreco == 2)].sum() / w[mc_].sum())}
                for k, c in enumerate((MUTED, ORANGE, BLUE)):
                    ax.bar(i + (k - 1) * 0.26, fr[k], width=0.24, color=c, label=["0 detected", "1 detected", "both detected"][k] if i == 0 else None)
                    ax.text(i + (k - 1) * 0.26, fr[k] + 0.015, f"{100 * fr[k]:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK)
            ax.set_xticks(range(len(cats)))
            ax.set_xticklabels([f"recoil pair {c}" for c in cats], fontsize=11)
            ax.set_ylim(0, 1.05)
            ax.text(0.97, 0.95, lab, transform=ax.transAxes, ha="right", va="top", fontsize=13, fontweight="bold")
            frame(ax, minor=False)
        axes[0, 0].set_ylabel("Fraction of 3N events", fontsize=12)
        axes[0, 0].legend(frameon=False, fontsize=10, loc="upper left")
        fig.tight_layout()
        fig.savefig(outdir / "fig3N_recoil_tagging.pdf")
        plt.close(fig)
        acceptance_grid(lambda r: f"{r}_p_irf", np.linspace(0.3, 1.2, 19), r"Recoil Momentum (Ion Rest Frame) [GeV/$c$]",
                        ("rec", "rec2"), "fig3N_recoil_acceptance_vs_pirf.pdf")
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        for lab in data:
            print(f"{lab}: {summary[lab]['3N']}")
        print(f"wrote {outdir}/fig3N_recoil_tagging.pdf fig3N_recoil_acceptance_vs_pirf.pdf summary.json")
        return

    if run_kind == "DIS":
        # Tagged DIS off SRC pairs: recoil acceptance vs the spectator light-cone fraction and momentum,
        # and the tagged yield over (x', Q2) for 10 fb^-1 x transparency
        acceptance_grid(lambda r: "alpha_s", np.linspace(0.5, 1.9, 15), r"Spectator Light-cone Fraction $\alpha_s$",
                        ("rec",), "figDIS_recoil_acceptance_vs_alpha.pdf")
        acceptance_grid(lambda r: f"{r}_p_irf", np.linspace(0.3, 1.2, 19), r"Spectator Momentum (Ion Rest Frame) [GeV/$c$]",
                        ("rec",), "figDIS_recoil_acceptance_vs_pirf.pdf")
        fig, axes = plt.subplots(1, len(data), figsize=(5.0 * len(data), 3.9), sharey=True, squeeze=False)
        for ax, (lab, (d, _)) in zip(axes[0], data.items()):
            tag = (d["rec_p_irf"] > P_IRF_MIN) & detected(d, "rec", "hit")
            H, xe, ye = np.histogram2d(d["xprime"][tag], d["q2"][tag], bins=[np.linspace(0.2, 1.0, 33), np.linspace(3, 25, 45)],
                                       weights=d["w_nb"][tag] * scale)
            im = ax.pcolormesh(xe, ye, np.ma.masked_less(H.T, 1), cmap="Blues", norm=matplotlib.colors.LogNorm(vmin=1), rasterized=True)
            ax.text(0.95, 0.93, lab, transform=ax.transAxes, ha="right", va="top", fontsize=13, fontweight="bold")
            ax.set_xlabel(r"x$'$", fontsize=13, loc="right")
            frame(ax)
            summary[lab]["DIS"] = {"sigma_nb_pIRF>0.3": float(d["w_nb"][d["rec_p_irf"] > P_IRF_MIN].sum()),
                                   "tagged_fraction_hit": float(d["w_nb"][tag].sum() / d["w_nb"][d["rec_p_irf"] > P_IRF_MIN].sum()),
                                   "tagged_events_10fb_x_transparency": float(d["w_nb"][tag].sum() * scale)}
        axes[0, 0].set_ylabel(r"Q$^2$ [GeV$^2$/c$^4$]", fontsize=13, loc="top")
        fig.colorbar(im, ax=axes[0, -1], label=rf"tagged events / bin ({lumi_nb_inv / 1e6:g} fb$^{{-1}}$ $\times$ {transparency:g})")
        fig.tight_layout()
        fig.savefig(outdir / "figDIS_xprime_q2_tagged.pdf")
        plt.close(fig)
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        for lab in data:
            print(f"{lab}: {summary[lab]['DIS']}")
        print(f"wrote {outdir}/figDIS_recoil_acceptance_vs_alpha.pdf figDIS_recoil_acceptance_vs_pirf.pdf figDIS_xprime_q2_tagged.pdf summary.json")
        return

    # Fig. 3: angle w.r.t. the ion beam of recoil (green) and leading (red) nucleons, 10 fb^-1, no p_IRF cut
    fig, axes = plt.subplots(len(data), 1, figsize=(6.6, 3.9 * len(data)), squeeze=False)
    bins = np.arange(0, 51, 1.0)
    cen, hw3 = 0.5 * (bins[1:] + bins[:-1]), 0.5 * np.diff(bins)
    for ax, (lab, (d, _)) in zip(axes[:, 0], data.items()):
        w = d["w_nb"] * lumi_nb_inv
        curves, errs = {}, {}
        for role, c, mk in (("rec", green, "o"), ("lead", red, "^")):
            x = d[f"{role}_theta_mrad"]
            h = np.histogram(x, bins, weights=w)[0] / 1e3
            e = np.sqrt(np.histogram(x, bins, weights=w**2)[0]) / 1e3
            g = h > 0
            ax.errorbar(cen[g], h[g], xerr=hw3[g], yerr=e[g], fmt=mk, ms=4, color=c, elinewidth=1, capsize=0)
            curves[role], errs[role] = h, e
        if "evap_pdg" in d:  # evaporation nucleons, scaled by 1/3 as in the paper
            nuc = np.isin(d["evap_pdg"], (2212, 2112))
            we = w[d["evap_event"][nuc]] / 3
            h = np.histogram(d["evap_theta_mrad"][nuc], bins, weights=we)[0] / 1e3
            e = np.sqrt(np.histogram(d["evap_theta_mrad"][nuc], bins, weights=we**2)[0]) / 1e3
            g = h > 0
            ax.errorbar(cen[g], h[g], xerr=hw3[g], yerr=e[g], fmt="s", ms=3.5, color=blue_ev, elinewidth=1, capsize=0)
            curves["evap"], errs["evap"] = h, e
        ymax = 1.15 * max(float(np.max(v)) for v in curves.values())
        ax.set_xlim(0, 50)
        ax.set_ylim(0, ymax)
        frame(ax)
        for x0 in (20, 30):
            ax.axvline(x0, color=orange if x0 == 20 else magenta, lw=1, ls=":")
        y_arrow = 0.965 * ymax
        ax.annotate("", xy=(0.3, y_arrow), xytext=(19.7, y_arrow), arrowprops=dict(arrowstyle="<->", color=orange, lw=1.8))
        ax.annotate("", xy=(49.5, y_arrow), xytext=(30.3, y_arrow), arrowprops=dict(arrowstyle="->", color=magenta, lw=1.8))
        ax.text(4.2, 1.01, "far-forward region", transform=ax.get_xaxis_transform(), ha="left", va="bottom",
                color=orange, fontsize=12, fontweight="bold")
        ax.text(40, 1.01, "central hadron endcap", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=magenta, fontsize=12, fontweight="bold")
        ax.text(0.0, 1.01, r"$\times10^3$", transform=ax.transAxes, ha="left", va="bottom", fontsize=12)
        zdc = 4.0
        ax.plot([zdc, zdc], [-0.07, 0.04], transform=ax.get_xaxis_transform(), color=brown, lw=2, clip_on=False)
        ax.text(zdc, -0.2, "ZDC limit", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=brown, fontsize=12, fontweight="bold")
        ax.text(0.92, 0.84, lab, transform=ax.transAxes, ha="right", fontsize=12, fontweight="bold")
        # direct labels in free space next to each curve (data points incl. error bars kept clear)
        pts = []
        for role in curves:
            h, e = curves[role], errs[role]
            for y in (h, h + e, np.maximum(h - e, 0)):
                g = h > 0
                pts.append(np.column_stack([cen[g] / 50, y[g] / ymax]))
        pts = np.vstack(pts + [np.array([[0.92, 0.86], [0.85, 0.86]])])  # keep the setting tag clear
        taken = []
        labels = [("rec", green, "Recoil\nNucleon", lambda k: (cen[k] / 50 + 0.14, 0.62 * curves["rec"][k] / ymax)),
                  ("lead", red, "Leading\nNucleon", lambda k: (cen[k] / 50, 0.25))]
        if "evap" in curves:
            labels.append(("evap", blue_ev, "Evaporation\nNucleon\n" + r"$\times$1/3", lambda k: (cen[k] / 50 + 0.08, 0.5)))
        for role, c, text, anchor in labels:
            k = int(np.argmax(curves[role]))
            place_label(ax, text, c, pts, anchor(k), taken, fontsize=12, fontweight="bold")
        ax.set_ylabel("Events per bin", fontsize=13)
        ax.set_xlabel(r"$\theta$ [mrad]", fontsize=13, loc="right")
    fig.tight_layout(h_pad=2.5)
    fig.savefig(outdir / "fig3_angles.pdf")
    plt.close(fig)

    # Fig. 4: momentum in the ion rest frame vs opening angle to the ion beam (jet map, empty bins white):
    # recoil nucleons, and evaporation nucleons on the left when the runs have them
    def fig4_panel(ax, theta, p_irf, w, text, lab, xmax):
        H, xe, ye = np.histogram2d(theta, p_irf, bins=[np.linspace(0, xmax, 61), np.linspace(0, 0.8, 41)], weights=w)
        ax.pcolormesh(xe, ye, np.ma.masked_less_equal(H, 0).T, cmap="jet", vmin=0, shading="flat", rasterized=True)
        ax.axhline(P_IRF_MIN, color="red", lw=2, ls=(0, (5, 3)))
        ax.text(0.05, 0.9, lab, transform=ax.transAxes, fontsize=15, va="top", color="#333333",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
        ax.text(0.95, 0.08, text, transform=ax.transAxes, fontsize=13, ha="right", va="bottom", color="#333333")
        ax.set_xticks(np.arange(0, xmax + 1, 5 if xmax >= 10 else 2))
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 0.8)
        ax.set_xlabel("Opening Angle [mrad]", fontsize=12, color="#333333")
        ax.set_ylabel(r"Momentum (IRF) [GeV/$c$]", fontsize=12, color="#333333")
        ax.tick_params(labelsize=11, colors="#333333")

    evap4 = all("evap_pdg" in d for d, _ in data.values())  # the paper's left column: evaporation nucleons
    fig, axes = plt.subplots(len(data), 2 if evap4 else 1, figsize=(9.2 if evap4 else 4.6, 2.9 * len(data)), squeeze=False)
    for row, (lab, (d, _)) in enumerate(data.items()):
        w = d["w_nb"] * lumi_nb_inv
        panels = [(axes[row, -1], d["rec_theta_mrad"], d["rec_p_irf"], w, "Recoil\nNucleons")]
        if evap4:
            nuc = np.isin(d["evap_pdg"], (2212, 2112))
            panels.insert(0, (axes[row, 0], d["evap_theta_mrad"][nuc], d["evap_p_irf"][nuc], w[d["evap_event"][nuc]],
                              "Evaporation\nNucleons"))
        for ax, th, pirf, wt, text in panels:
            fig4_panel(ax, th, pirf, wt, text, lab, 15 if low_energy[lab] else 8)
    fig.tight_layout()
    fig.savefig(outdir / "fig4_pirf_vs_theta.pdf")
    plt.close(fig)

    # Fig. 5: Q2 of events with a recoil detected far forward, 10 fb^-1 x transparency, one point per GeV^2
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    bins = np.arange(3, 25.5, 1.0)
    cen5 = 0.5 * (bins[1:] + bins[:-1])
    curves5 = {}
    for lab, (d, _) in data.items():
        m = (d["rec_p_irf"] > P_IRF_MIN) & detected(d, "rec", "hit")
        h, _ = np.histogram(d["q2"][m], bins, weights=d["w_nb"][m] * scale)
        summary[lab]["q2_reach_100_per_GeV2"] = float(bins[1:][h >= 100].max()) if (h >= 100).any() else None
        c, mk = (red, "o") if low_energy[lab] else ("#1029ff", "s")
        g = h > 0
        ax.plot(cen5[g], h[g], ls="none", marker=mk, ms=5, color=c)
        curves5[lab] = (h, c)
        if "e_reco" in d:  # open markers: also the electron identified by EICrecon (ReconstructedElectrons)
            me = m & d["e_reco"]
            he = np.histogram(d["q2"][me], bins, weights=d["w_nb"][me] * scale)[0]
            summary[lab]["q2_reach_100_per_GeV2_e_identified"] = float(bins[1:][he >= 100].max()) if (he >= 100).any() else None
            summary[lab]["recoil_tagged_fraction_e_identified"] = float(he.sum() / h.sum()) if h.sum() else None
            ax.plot(cen5[he > 0], he[he > 0], ls="none", marker=mk, ms=5, mfc="none", color=c)
    ax.set_yscale("log")
    ax.set_xlim(1, 25)
    lo = min(float(np.min(h[h > 0])) for h, _ in curves5.values())
    hi = max(float(np.max(h)) for h, _ in curves5.values())
    ax.set_ylim(lo / 2, hi * 2)
    # direct labels: the higher curve labeled above itself, the lower below
    order = sorted(curves5, key=lambda l: -float(np.sum(curves5[l][0])))
    for i, lab in enumerate(order):
        h, c = curves5[lab]
        k = 12 if i == 0 else 3
        ax.text(cen5[k], h[k] * (9 if i == 0 else 0.2), lab, color=c, fontsize=15, ha="center", va="center")
    if any("e_reco" in d for d, _ in data.values()):
        ax.text(0.97, 0.95, "filled: no electron requirement (as PRC 105)\nopen: electron identified", transform=ax.transAxes,
                fontsize=9, ha="right", va="top", color="#444444")
    frame(ax)
    ax.set_ylabel("Events per bin", fontsize=13, loc="top")
    ax.set_xlabel(r"Q$^2$ [GeV$^2$/c$^4$]", fontsize=13, loc="right")
    fig.tight_layout()
    fig.savefig(outdir / "fig5_q2.pdf")
    plt.close(fig)

    # Fig. 6 in the style of PRC 105, 034001: rows are beam settings, columns recoil p and n, panels
    # share edges. Generated recoils over the full p_IRF range (black crosses), accepted at hit level
    # above the p_IRF cut (filled), accepted at reco level (blue crosses, in the slot of the paper's
    # smeared points), all for 10 fb^-1 x transparency. fig6b: acceptance vs p_IRF.
    pink, pblue = "#f4c4c4", "#1f3fe0"
    lbins = np.linspace(0.0, 1.3, 40)
    lc, hw = 0.5 * (lbins[1:] + lbins[:-1]), 0.5 * np.diff(lbins)
    hists = {}
    for lab, (d, _) in data.items():
        for code, name in ((2212, "proton"), (2112, "neutron")):
            m = d["rec_type"] == code
            x, w = d["rec_p_irf"][m], d["w_nb"][m] * scale
            cut = x > P_IRF_MIN
            sel = {"gen": np.ones(len(x), bool), "hit": cut & detected(d, "rec", "hit")[m],
                   "reco": cut & detected(d, "rec", "reco")[m], "cut": cut}
            hists[lab, name] = {k: (np.histogram(x[v], lbins, weights=w[v])[0],
                                    np.sqrt(np.histogram(x[v], lbins, weights=w[v] ** 2)[0])) for k, v in sel.items()}
            hists[lab, name]["w2"] = {k: np.histogram(x[cut & ~v], lbins, weights=w[cut & ~v] ** 2)[0] for k, v in sel.items()
                                      if k in ("hit", "reco")}
            # "smeared" spectrum: reconstructed momentum boosted back to the ion rest frame, for the
            # reco-level accepted events (PRC 105, 034001 Fig. 6 blue points; there a parametrized smearing)
            if "rec_reco_p_irf" in d:
                xr = d["rec_reco_p_irf"][m]
                ok = sel["reco"] & np.isfinite(xr)
                hists[lab, name]["smear"] = (np.histogram(xr[ok], lbins, weights=w[ok])[0],
                                             np.sqrt(np.histogram(xr[ok], lbins, weights=w[ok] ** 2)[0]))
                hists[lab, name]["smear_coverage"] = float(w[ok].sum() / max(w[sel["reco"]].sum(), 1e-300))
            summary[lab][f"recoil {name[0]}"] = {k: float(w[sel[k]].sum() / w[cut].sum()) for k in ("hit", "reco")}

    def paper_grid(ylabel):
        n = len(data)
        fig, axes = plt.subplots(n, 2, figsize=(8.4, 2.8 * n + 0.9), sharey="row", squeeze=False,
                                 gridspec_kw={"wspace": 0, "hspace": 0})
        for row in range(n):
            for col in range(2):
                ax = axes[row, col]
                ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11)
                ax.set_xlim(0, 1.3)
                ax.set_xticks([0.0, 0.5, 1.0] if col == 0 else [0.5, 1.0])
                if row < n - 1:
                    ax.tick_params(labelbottom=False)
            axes[row, 0].set_ylabel(ylabel, fontsize=12)
        fig.subplots_adjust(left=0.13, right=0.98, top=1 - 0.45 / fig.get_figheight(), bottom=0.55 / fig.get_figheight())
        fig.supxlabel(r"Spectator Momentum (Ion Rest Frame) [GeV/$c$]", fontsize=12, y=0.1 / fig.get_figheight())
        return fig, axes

    def labels(ax, letter, text):
        ax.text(0.05, 0.93, f"({letter})", transform=ax.transAxes, va="top", fontsize=13)
        ax.text(0.95, 0.93, text, transform=ax.transAxes, va="top", ha="right", fontsize=12)

    def spectra(log, letters, fname):
        fig, axes = paper_grid("Events per bin")
        top = max(float(np.max(hists[k]["gen"][0] + hists[k]["gen"][1])) for k in hists)
        k10 = 0 if log else int(np.floor(np.log10(top)))
        f = 10.0 ** k10
        it = iter(letters)
        for row, lab in enumerate(data):
            ymax = 0.0
            for col, name in enumerate(("proton", "neutron")):
                ax, H = axes[row, col], hists[lab, name]
                ax.stairs(H["hit"][0] / f, lbins, fill=True, color=pink, lw=0, label="Accepted", zorder=1)
                blue = "smear" if "smear" in H else "reco"
                for key, c, lbl in (("gen", INK, "Generated"), (blue, pblue, "Reconstructed")):
                    v, e = H[key]
                    g = v > 0
                    ax.errorbar(lc[g], v[g] / f, xerr=hw[g], yerr=e[g] / f, fmt="none", ecolor=c, elinewidth=1.4,
                                label=lbl, zorder=2)
                    if g.any():
                        ymax = max(ymax, float(np.max(v[g] + e[g])) / f)
                labels(ax, next(it), f"{name}, {lab}")
                cov = H.get("smear_coverage", 1.0)
                if cov < 0.8:
                    why = "no OMD/RP reconstruction\nfor ion beams" if name == "proton" else "no reconstructed\nZDC neutron"
                    ax.text(0.95, 0.55, f"reco momentum for {100 * cov:.0f}%\nof accepted ({why})",
                            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=pblue)
            if log:
                axes[row, 0].set_yscale("log")
                axes[row, 0].set_ylim(1, ymax * 8)
            else:
                axes[row, 0].set_ylim(0, ymax * 1.3)
                if row > 0:  # rows touch: drop the top tick label so it does not collide with the row above
                    from matplotlib.ticker import MaxNLocator
                    axes[row, 0].yaxis.set_major_locator(MaxNLocator(nbins=6, prune="upper"))
        if not log and k10:
            axes[0, 0].text(0.0, 1.02, rf"$\times10^{{{k10}}}$", transform=axes[0, 0].transAxes, fontsize=12, va="bottom")
        h, l_ = axes[0, 1].get_legend_handles_labels()
        order = [l_.index(x) for x in ("Generated", "Accepted", "Reconstructed")]
        h, l_ = [h[i] for i in order], [l_[i] for i in order]
        if log:  # the paper's log panels carry no legend; put it in one row above the grid
            fig.legend(h, l_, frameon=False, fontsize=11, ncol=3, loc="upper center", bbox_to_anchor=(0.555, 1.0))
        else:
            axes[0, 1].legend(h, l_, frameon=False, fontsize=11, loc="upper right", bbox_to_anchor=(1.0, 0.85))
        fig.savefig(outdir / fname)
        plt.close(fig)

    spectra(True, "abcdefgh", "fig6_recoil_acceptance.pdf")
    spectra(False, "efghijkl", "fig6_recoil_acceptance_linear.pdf")

    # Fig. 2: momentum resolution of the reconstructed recoil nucleons vs lab momentum (B0 tracks for
    # protons, ZDC for neutrons); the paper used Yellow-Report parametrizations here
    if all("rec_reco_p_lab_mag" in d for d, _ in data.values()):
        fig, ax = plt.subplots(figsize=(6.6, 4.0))
        pb = np.arange(10, 205, 10.0)
        for lab, (d, _) in data.items():
            # one colour per setting and species, so the two settings' neutrons stay apart (#204)
            for code, name, col, mk in ((2112, "Neutron", "#1029ff" if low_energy[lab] else "#00a0c0", "s"),
                                        (2212, "Proton", red if low_energy[lab] else "#d2a679", "o")):
                m = (d["rec_final_pdg"] == code) & np.isfinite(d["rec_reco_p_lab_mag"]) & (d["rec_p_irf"] > P_IRF_MIN)
                pt, rel = d["rec_true_p_lab_mag"][m], (d["rec_reco_p_lab_mag"][m] - d["rec_true_p_lab_mag"][m]) / d["rec_true_p_lab_mag"][m]
                xs, ys, es = [], [], []
                for lo, hi in zip(pb[:-1], pb[1:]):
                    k = (pt >= lo) & (pt < hi)
                    if k.sum() >= 30:
                        q16, q84 = np.percentile(rel[k], [16, 84])
                        sig = 0.5 * (q84 - q16)
                        xs.append(0.5 * (lo + hi)); ys.append(sig); es.append(sig / np.sqrt(2 * k.sum()))
                if xs:
                    ax.errorbar(xs, ys, yerr=es, fmt=mk, ms=5, color=col, label=f"{name} ({lab})", capsize=0)
                    summary[lab][f"recoil {name.lower()} dp/p (68% half-width)"] = float(np.median(ys))
                    summary[lab][f"recoil {name.lower()} dp/p bias (median)"] = float(np.median(rel)) if len(rel) else None
        ax.set_xlim(0, 210)
        ax.set_ylim(0, 1.6 * ax.get_ylim()[1])  # headroom: the note and the legend sit above the points
        frame(ax)
        ax.set_xlabel(r"p [GeV/$c$]", fontsize=13, loc="right")
        ax.set_ylabel(r"$\Delta$p/p (68% half-width)", fontsize=13, loc="top")
        ax.text(0.03, 0.95, r"recoil nucleons, $p_\mathrm{IRF} > 0.3$ GeV/$c$; protons: B0 tracks, neutrons: ZDC",
                transform=ax.transAxes, fontsize=9, va="top", color=MUTED)
        ax.legend(frameon=False, fontsize=10, loc="upper right", bbox_to_anchor=(1.0, 0.9))
        fig.tight_layout()
        fig.savefig(outdir / "fig2_resolution.pdf")
        plt.close(fig)

    # fig6b: acceptance vs p_IRF (accepted / generated above the p_IRF cut), hit and reco level
    fig, axes = paper_grid("Acceptance")
    it = iter("abcdefgh")
    for row, lab in enumerate(data):
        for col, name in enumerate(("proton", "neutron")):
            ax, H = axes[row, col], hists[lab, name]
            gen = H["cut"][0]
            ok = gen > 0
            for j, (key, c, lbl) in enumerate((("hit", INK, "hit level"), ("reco", pblue, "reco level"))):
                acc = H[key][0]
                eff = np.divide(acc, gen, out=np.zeros_like(acc), where=ok)
                err = np.sqrt(H[key][1] ** 2 * (1 - eff) ** 2 + H["w2"][key] * eff ** 2) / np.where(ok, gen, 1)
                ax.errorbar(lc[ok] + (j - 0.5) * 0.25 * hw[ok], eff[ok], yerr=err[ok], fmt="o", ms=3.5, color=c,
                            elinewidth=1.2, label=lbl)
            ax.axvline(P_IRF_MIN, color=MUTED, lw=0.8, ls=":")
            labels(ax, next(it), f"{name}, {lab}: {100 * summary[lab]['recoil ' + name[0]]['hit']:.0f}%")
        axes[row, 0].set_ylim(0, 1.15)
    axes[0, 1].legend(frameon=False, fontsize=11, loc="upper right", bbox_to_anchor=(1.0, 0.85))
    fig.savefig(outdir / "fig6b_acceptance_vs_pirf.pdf")
    plt.close(fig)

    # Evaporation (runs with the Geant4 de-excitation of the A-2 remnant): Fig. 3 with the evaporation
    # nucleons as in PRC 105, 034001 (there from BeAGLE), their rest-frame momenta, where they land in
    # the far-forward detectors, and how often they share the ZDC with a recoil neutron
    if all("evap_pdg" in d for d, _ in data.values()):
        blue_ev = "#1f3fe0"
        fig, axes = plt.subplots(len(data), 2, figsize=(11.0, 3.6 * len(data)), squeeze=False)
        for row, (lab, (d, _)) in enumerate(data.items()):
            w = d["w_nb"] * lumi_nb_inv
            ev, pdg = d["evap_event"], d["evap_pdg"]
            nuc = np.isin(pdg, (2212, 2112))
            we = w[ev]
            ax = axes[row, 0]
            tb = np.arange(0, 20.25, 0.25) if low_energy[lab] else np.arange(0, 8.1, 0.1)
            tc = 0.5 * (tb[1:] + tb[:-1])
            for x, wt, c, name in ((d["rec_theta_mrad"], w, green, "Recoil nucleon"), (d["lead_theta_mrad"], w, red, "Leading nucleon"),
                                   (d["evap_theta_mrad"][nuc & (pdg == 2212)], we[nuc & (pdg == 2212)], blue_ev, "Evaporation protons"),
                                   (d["evap_theta_mrad"][pdg == 2112], we[pdg == 2112], "#00a0c0", "Evaporation neutrons")):
                h = np.histogram(x, tb, weights=wt)[0]
                ax.step(tc, h, where="mid", color=c, lw=1.6, label=name)
            ax.axvline(4.0, color=brown, lw=1, ls=":")
            ax.text(4.0, 1.01, "ZDC limit", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color=brown, fontsize=10)
            ax.set_yscale("log")
            ax.set_xlim(tb[0], tb[-1])
            frame(ax)
            ax.set_xlabel(r"$\theta$ w.r.t. ion beam [mrad]", fontsize=12, loc="right")
            ax.set_ylabel("Particles per bin (10 fb$^{-1}$)", fontsize=11)
            ax.text(0.03, 0.95, lab, transform=ax.transAxes, ha="left", va="top", fontsize=12, fontweight="bold")
            lo_, hi_ = ax.get_ylim()
            ax.set_ylim(lo_, hi_ * 60)  # headroom for the legend
            ax.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
            ax = axes[row, 1]
            pb = np.linspace(0, 0.6, 61)
            for code, c, name in ((2212, blue_ev, "Evaporation protons"), (2112, "#00a0c0", "Evaporation neutrons")):
                m = pdg == code
                ax.stairs(np.histogram(d["evap_p_irf"][m], pb, weights=we[m])[0], pb, color=c, lw=1.6, label=name)
            for code, c, name in ((2212, green, "Recoil protons"), (2112, "#7d9a3a", "Recoil neutrons")):
                m = d["rec_type"] == code
                ax.stairs(np.histogram(d["rec_p_irf"][m], pb, weights=w[m])[0], pb, color=c, lw=1.2, ls="--", label=name)
            ax.axvline(P_IRF_MIN, color=MUTED, lw=1, ls=":")
            ax.set_yscale("log")
            frame(ax)
            ax.set_xlabel(r"Momentum (Ion Rest Frame) [GeV/$c$]", fontsize=12, loc="right")
            lo_, hi_ = ax.get_ylim()
            ax.set_ylim(lo_, hi_ * 60)
            ax.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
            # summary: multiplicities, detection, and the overlap with the recoil sample
            n_ev = float(w.sum())
            det_p = (d["evap_B0TrackerHits"] >= 3) | (d["evap_ForwardOffMTrackerHits"] >= 2) | (d["evap_ForwardRomanPotHits"] >= 2)
            det_n = d["evap_zdc_edep"] > 0
            sp, sn = pdg == 2212, pdg == 2112
            frac = lambda m_, ok: float(we[m_ & ok].sum() / max(we[m_].sum(), 1e-300))  # noqa: E731
            zdc_evap = np.zeros(len(w), bool)
            np.logical_or.at(zdc_evap, ev[sn & det_n], True)
            recn = (d["rec_type"] == 2112) & (d["rec_p_irf"] > P_IRF_MIN) & detected(d, "rec", "hit")
            summary[lab]["evaporation"] = {
                "neutrons_per_event": float(we[sn].sum() / n_ev), "protons_per_event": float(we[sp].sum() / n_ev),
                "products_per_event_by_pdg": {str(int(c_)): float(we[pdg == c_].sum() / n_ev) for c_ in np.unique(pdg)},
                "evap_protons_B0": frac(sp, d["evap_B0TrackerHits"] >= 3), "evap_protons_OMD": frac(sp, d["evap_ForwardOffMTrackerHits"] >= 2),
                "evap_protons_RP": frac(sp, d["evap_ForwardRomanPotHits"] >= 2), "evap_protons_any": frac(sp, det_p),
                "evap_neutrons_ZDC": frac(sn, det_n),
                "evap_nucleons_above_pIRF_cut": frac(nuc, d["evap_p_irf"] > P_IRF_MIN),
                "events_with_evap_neutron_in_ZDC": float(w[zdc_evap].sum() / n_ev),
                "tagged_recoil_neutrons_sharing_ZDC_with_evap_neutron": float(w[recn & zdc_evap].sum() / max(w[recn].sum(), 1e-300)),
            }
        fig.tight_layout()
        fig.savefig(outdir / "figEvap_angles_momenta.pdf")
        plt.close(fig)
        for lab in data:
            print(f"{lab}: evaporation {summary[lab]['evaporation']}")


    # Pair-type table: recoil-tagged and triple-coincidence fractions
    lines = []
    for lab, (d, metas) in data.items():
        m0 = d["rec_p_irf"] > P_IRF_MIN
        w = d["w_nb"]
        rec_ok = detected(d, "rec", "hit")
        lead_ok = detected_anywhere(d, "lead")
        e_ok = d["e_truth_matched"] if "e_truth_matched" in d else np.ones(len(w), bool)
        sigma = float(w[m0].sum())
        summary[lab]["sigma_nb_cuts_pIRF"] = sigma
        summary[lab]["electron_truth_matched"] = float(w[m0 & e_ok].sum() / w[m0].sum())
        summary[lab]["runs"] = [m_["run"] for m_ in metas]
        for lt, rt, pair in ((2212, 2212, "pp"), (2212, 2112, "pn"), (2112, 2212, "np"), (2112, 2112, "nn")):
            m = m0 & (d["lead_type"] == lt) & (d["rec_type"] == rt)
            if not m.any():
                continue
            f = lambda ok: float(w[m & ok].sum() / w[m].sum())  # noqa: E731
            summary[lab][f"pair_{pair}"] = {"weight_fraction": float(w[m].sum() / sigma), "recoil": f(rec_ok),
                                            "lead": f(lead_ok), "triple_e_lead_rec": f(rec_ok & lead_ok & e_ok)}
            lines.append(f"{lab:8s} {pair}  share {w[m].sum() / sigma:5.3f}  recoil {f(rec_ok):5.2f}  "
                         f"lead {f(lead_ok):5.2f}  e+lead+recoil {f(rec_ok & lead_ok & e_ok):5.2f}")
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (outdir / "pairs.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for lab in data:
        print(f"{lab}: recoil p {summary[lab]['recoil p']}, recoil n {summary[lab]['recoil n']}, "
              f"Q2 reach (>=100/GeV2) {summary[lab]['q2_reach_100_per_GeV2']}, "
              f"electron {summary[lab]['electron_truth_matched']:.2f}")
    print(f"wrote {outdir}/fig3_angles.pdf fig4_pirf_vs_theta.pdf fig5_q2.pdf fig2_resolution.pdf fig6_recoil_acceptance.pdf fig6_recoil_acceptance_linear.pdf fig6b_acceptance_vs_pirf.pdf summary.json pairs.txt")


# ----------------------------------------------------------------------------- contact uncertainty band

# GCF k-space contacts from the Scale and Scheme paper (Cruz-Torres et al., Nature Phys. 17, 306,
# doi:10.1038/s41567-020-01053-7) with their uncertainties, as set in
# src/libraries/nucleus/gcfSRC.cc (set_Contacts_SS_k, AV18), per A/2; gcfSRC::randomize_Contacts
# varies each one independently with a Gaussian of this width
CONTACTS_SS_K = {"C12": {"Cpp0": (1.140, 0.210), "Cnn0": (1.140, 0.210), "Cpn0": (1.244, 0.319), "Cpn1": (15.876, 1.770)}}
GEV_FM = 0.1973  # src/include/constants.hh


def av18_phi_sq():
    """AV18 pn spin-0 and spin-1 |phi|^2 tables (src/libraries/nucleus/universal_functions/AV18.hh)."""
    import re
    text = (Path(__file__).resolve().parent.parent / "src/libraries/nucleus/universal_functions/AV18.hh").read_text()
    tab = {}
    for name in ("AV18_pn0", "AV18_pn1"):
        body = re.search(name + r"\[100\]\s*=\s*\{([^}]*)\}", text).group(1)
        tab[name] = np.array([float(x) for x in body.replace("\n", "").split(",") if x.strip()])
        assert len(tab[name]) == 100
    return tab["AV18_pn0"], tab["AV18_pn1"]


def phi_sq(table, k):
    """gcfSRC::get_phiSq without its constant normalization (ratios at equal k only)."""
    b = np.asarray(k) / GEV_FM / 0.1
    i = np.clip(np.floor(b).astype(int), 1, 99)
    x = b - i
    v = x * table[i] + (1 - x) * table[i - 1]
    return np.where((b >= 1) & (b <= 100), v, 0.0)


def contact_basis(d, k):
    """Per event, the coefficients b_c of the contact-reweighting factor f = sum_c C'_c b_c (f = 1 at the
    nominal contacts): the weight goes with S(k_rel) = 2 C_pp0 phi_pp0^2, 2 C_nn0 phi_nn0^2 or
    C_pn0 phi_pn0^2 + C_pn1 phi_pn1^2 of the struck pair (gcfSRC::get_S)."""
    c = {n: v[0] for n, v in CONTACTS_SS_K["C12"].items()}
    t0, t1 = av18_phi_sq()
    f0, f1 = phi_sq(t0, k), phi_sq(t1, k)
    lt, rt = d["lead_type"], d["rec_type"]
    pp, nn = (lt == 2212) & (rt == 2212), (lt == 2112) & (rt == 2112)
    pn = ~pp & ~nn
    b = np.zeros((len(k), 4))  # columns: Cpp0, Cnn0, Cpn0, Cpn1
    b[pp, 0] = 1 / c["Cpp0"]
    b[nn, 1] = 1 / c["Cnn0"]
    den = c["Cpn0"] * f0 + c["Cpn1"] * f1
    ok = pn & (den > 0)
    b[ok, 2] = f0[ok] / den[ok]
    b[ok, 3] = f1[ok] / den[ok]
    return b


def band(settings, outdir, lumi_nb_inv, transparency, n_rep, seed):
    """Contact-uncertainty band of the 12C PWIA results by exact reweighting of the simulated events."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scale = lumi_nb_inv * transparency
    names = list(CONTACTS_SS_K["C12"])
    nom = np.array([CONTACTS_SS_K["C12"][n][0] for n in names])
    sig = np.array([CONTACTS_SS_K["C12"][n][1] for n in names])
    rng = np.random.default_rng(seed)
    reps = np.clip(nom + sig * rng.standard_normal((n_rep, 4)), 0, None)
    pct = lambda a: np.percentile(a, [16, 50, 84], axis=0)  # noqa: E731
    summary = {"contacts_per_A_over_2": {n: {"value": v[0], "sigma": v[1]} for n, v in CONTACTS_SS_K["C12"].items()},
               "replicas": n_rep, "seed": seed, "method": "independent Gaussian contacts (gcfSRC::randomize_Contacts), "
               "exact per-event reweighting S'(k_rel)/S(k_rel); 16th/50th/84th percentiles over replicas"}
    red, blue = "#f0281e", "#1029ff"
    bins = np.arange(3, 25.5, 1.0)
    cen = 0.5 * (bins[1:] + bins[:-1])
    lbins = np.linspace(0.0, 1.3, 40)
    lc = 0.5 * (lbins[1:] + lbins[:-1])
    fig5, ax5 = plt.subplots(figsize=(6.6, 3.8))
    fig6, ax6 = plt.subplots(len(settings), 2, figsize=(8.4, 2.8 * len(settings) + 0.9), sharey="row", squeeze=False,
                             gridspec_kw={"wspace": 0, "hspace": 0})
    figr, axr = plt.subplots(1, 2, figsize=(10.0, 3.8), sharey=True)
    for row, (lab, dirs) in enumerate(settings.items()):
        d, metas = load(dirs)
        if any(m_.get("kind", "2N") != "2N" for m_ in metas) or metas[0]["beams"]["ion"] != "C12":
            sys.exit("band: 12C SRC-pair (genQE) runs only")
        if any(has_fsi(r) for r in dirs):
            sys.exit("band: PWIA runs only (the reweighting needs pRel; GENIE FSI runs take no flat transparency)")
        k = np.concatenate([np.linalg.norm(uproot.open(Path(r) / "gen.root")["genT"]["pRel"].array(library="np")
                                           [np.load(Path(r) / "analysis_events.npz")["num"]], axis=1) for r in dirs])
        B = contact_basis(d, k)
        assert np.allclose(B @ nom, 1.0), "reweighting factor is not 1 at the nominal contacts"
        w = d["w_nb"]
        lin = lambda sel, x=None, bb=None: (np.stack([np.histogram(x[sel], bb, weights=(w * B[:, c])[sel])[0] for c in range(4)])  # noqa: E731
                                            if x is not None else np.array([(w * B[:, c])[sel].sum() for c in range(4)]))
        cut = d["rec_p_irf"] > P_IRF_MIN
        tag = cut & detected(d, "rec", "hit")
        out = {}
        # cross section after the cuts and the p_IRF cut, and the recoil-tagged yields for 10 fb^-1 x transparency
        out["sigma_nb_cuts_pIRF"] = pct(reps @ lin(cut))
        for code, name in ((2212, "p"), (2112, "n")):
            m = d["rec_type"] == code
            out[f"recoil {name} tagged events"] = pct(reps @ lin(tag & m) * scale)
            out[f"recoil {name} acceptance (hit)"] = pct((reps @ lin(tag & m)) / (reps @ lin(cut & m)))
        for lt, rt, pair in ((2212, 2212, "pp"), (2212, 2112, "pn"), (2112, 2212, "np"), (2112, 2112, "nn")):
            m = cut & (d["lead_type"] == lt) & (d["rec_type"] == rt)
            out[f"share {pair}"] = pct((reps @ lin(m)) / (reps @ lin(cut)))
        # Fig. 5 with the band, and the Q2 reach per replica
        H = reps @ lin(tag, d["q2"], bins) * scale
        h_nom = nom @ lin(tag, d["q2"], bins) * scale
        reach = np.array([bins[1:][h >= 100].max() if (h >= 100).any() else np.nan for h in H])
        out["q2_reach_100_per_GeV2"] = np.nanpercentile(reach, [16, 50, 84])
        lo, hi = np.percentile(H, [16, 84], axis=0)
        c, mk = (red, "o") if metas[0]["beams"]["ion_momentum_per_nucleon_GeV"] < 75 else (blue, "s")
        g = h_nom > 0
        ax5.fill_between(cen[g], lo[g], hi[g], step=None, color=c, alpha=0.2, lw=0)
        ax5.plot(cen[g], h_nom[g], ls="none", marker=mk, ms=5, color=c, label=lab)
        ok = h_nom > 0
        axr[0].fill_between(cen[ok], (lo / h_nom - 1)[ok], (hi / h_nom - 1)[ok], color=c, alpha=0.25, lw=0, label=lab)
        # Fig. 6: accepted recoil spectra (hit level) with the band
        for col, (code, name) in enumerate(((2212, "proton"), (2112, "neutron"))):
            ax = ax6[row, col]
            m = tag & (d["rec_type"] == code)
            Hs = reps @ lin(m, d["rec_p_irf"], lbins) * scale
            hn = nom @ lin(m, d["rec_p_irf"], lbins) * scale
            l6, h6 = np.percentile(Hs, [16, 84], axis=0)
            ok = hn > 0
            ls_ = "-" if code == 2212 else "--"
            for y in (l6, h6):
                axr[1].plot(lc[ok], (y / hn - 1)[ok], ls=ls_, color=c, lw=1.6,
                            label=f"recoil {name}, {lab}" if y is l6 else None)
            ax.stairs(hn, lbins, color="#f4c4c4", fill=True, lw=0, label="Accepted (nominal contacts)")
            ax.fill_between(lc, l6, h6, step="mid", color=INK, alpha=0.25, lw=0, label="contact uncertainty (68%)")
            ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11)
            ax.set_xlim(0, 1.3)
            ax.set_xticks([0.0, 0.5, 1.0] if col == 0 else [0.5, 1.0])
            rel = (h6.sum() - l6.sum()) / 2 / max(hn.sum(), 1e-300)
            ax.text(0.95, 0.93, f"{name}, {lab}\n$\\pm${100 * rel:.0f}% (68%)", transform=ax.transAxes, va="top", ha="right", fontsize=11)
            if row < len(settings) - 1:
                ax.tick_params(labelbottom=False)
        ax6[row, 0].set_yscale("log")
        ax6[row, 0].set_ylim(1, None)
        ax6[row, 0].set_ylabel("Events per bin", fontsize=12)
        summary[lab] = {key: {"p16": float(v[0]), "median": float(v[1]), "p84": float(v[2])} for key, v in out.items()}
        summary[lab]["runs"] = [m_["run"] for m_ in metas]
        print(f"{lab}: " + ", ".join(f"{key} {v[1]:.4g} [{v[0]:.4g}, {v[2]:.4g}]" for key, v in out.items()))
    ax5.set_yscale("log")
    ax5.set_xlim(1, 25)
    ax5.tick_params(direction="in", which="both", top=True, right=True, labelsize=12)
    ax5.set_ylabel("Events per bin", fontsize=13, loc="top")
    ax5.set_xlabel(r"Q$^2$ [GeV$^2$/c$^4$]", fontsize=13, loc="right")
    ax5.axhline(100, color=MUTED, lw=0.8, ls=":")
    ax5.legend(frameon=False, fontsize=12, title="band: GCF contact uncertainty (68%)", title_fontsize=10)
    fig5.tight_layout()
    fig5.savefig(outdir / "figband_q2.pdf")
    for ax, xl in zip(axr, (r"Q$^2$ [GeV$^2$/c$^4$]", r"Spectator Momentum (Ion Rest Frame) [GeV/$c$]")):
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11)
        ax.set_xlabel(xl, fontsize=12, loc="right")
    axr[0].set_xlim(3, 25)
    axr[1].set_xlim(P_IRF_MIN, 1.2)
    axr[0].set_ylabel("relative contact uncertainty (68%)", fontsize=12)
    axr[0].text(0.03, 0.95, "recoil-tagged yield", transform=axr[0].transAxes, va="top", fontsize=11)
    axr[1].text(0.03, 0.95, "accepted recoils (hit level)", transform=axr[1].transAxes, va="top", fontsize=11)
    axr[0].legend(frameon=False, fontsize=10, loc="lower right")
    axr[1].legend(frameon=False, fontsize=9, loc="lower right", ncol=2)
    figr.tight_layout()
    figr.savefig(outdir / "figband_relative.pdf")
    fig6.subplots_adjust(left=0.13, right=0.98, top=1 - 0.45 / fig6.get_figheight(), bottom=0.55 / fig6.get_figheight())
    fig6.supxlabel(r"Spectator Momentum (Ion Rest Frame) [GeV/$c$]", fontsize=12, y=0.1 / fig6.get_figheight())
    h_, l_ = ax6[0, 1].get_legend_handles_labels()
    fig6.legend(h_, l_, frameon=False, fontsize=10, ncol=2, loc="upper center", bbox_to_anchor=(0.555, 1.0))
    fig6.savefig(outdir / "figband_recoil_spectra.pdf")
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {outdir}/figband_q2.pdf figband_recoil_spectra.pdf figband_relative.pdf summary.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract", help="reduce run directories to analysis_events.npz")
    e.add_argument("runs", nargs="+")
    e.add_argument("-j", "--jobs", type=int, default=16)
    p = sub.add_parser("plots", help="figures and summary from extracted runs")
    p.add_argument("--setting", action="append", required=True, help="LABEL=RUN_DIR[,RUN_DIR...] (Q2 slices of one setting)")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--lumi-fb", type=float, default=10.0, help="integrated luminosity per nucleus [fb^-1] (default 10)")
    p.add_argument("--transparency", type=float, help="flat nucleon transparency factor (default 0.5, or 1 for GENIE FSI "
                   "runs, from run.json)")
    p.add_argument("--cc1", choices=("fixed", "2021"), default="fixed",
                   help="cross section: fixed (current cc1, default) or 2021 (cc1 before the #33 fix, as PRC 105, "
                        "034001 was generated; per-event reweighting, 2N runs)")
    b = sub.add_parser("band", help="GCF contact-uncertainty band of 12C runs by reweighting (no new simulation)")
    b.add_argument("--setting", action="append", required=True, help="LABEL=RUN_DIR[,RUN_DIR...] (Q2 slices of one setting)")
    b.add_argument("-o", "--outdir", required=True)
    b.add_argument("--lumi-fb", type=float, default=10.0)
    b.add_argument("--transparency", type=float, default=0.5)
    b.add_argument("--replicas", type=int, default=1000)
    b.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.cmd == "extract":
        for r in a.runs:
            extract(r, a.jobs)
    else:
        settings = {}
        for s in a.setting:
            lab, dirs = s.split("=", 1)
            settings[lab] = dirs.split(",")
        if a.cmd == "band":
            band(settings, a.outdir, a.lumi_fb * 1e6, a.transparency, a.replicas, a.seed)
        else:
            plots(settings, a.outdir, a.lumi_fb * 1e6, a.transparency, a.cc1)


if __name__ == "__main__":
    main()
