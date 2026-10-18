#!/usr/bin/env python3
"""Convert GCF generator trees (target rest frame) to collider-frame HepMC3 for the EIC.

The generators work in the target rest frame with the electron beam along +z.
An EIC event with an electron of energy E_e colliding head-on with an ion of
momentum A * p_N (p_N per nucleon) is the same event in the rest frame with an
electron of energy

    E_IRF = gamma * (1 + beta) * E_e,    gamma * beta = A * p_N / M_A,

so run the generator at E_IRF (about 440 GeV for e+12C at 5x41 and 2362 GeV at
10x110; this script prints the value) and convert its output here. The
conversion rotates by pi about y, so that the electron moves along -z (ePIC
convention), and boosts along +z with the ion velocity. The output is head-on,
with the electron and the ion as status-4 beams. Feed it to the afterburner
(abconv) for the crossing angle and beam effects, then to npsim:

    tools/gcf2hepmc.py gen.root gen.hepmc3 --beams 5x41 --ion C12 --paper-cuts
    abconv -p ip6_eAu_41x5 -f hepmc3 -o ab gen.hepmc3
    npsim --compactFile geometry.xml --inputFiles ab.hepmc --outputFile sim.edm4hep.root -N -1

Supported inputs (detected from the tree):
  * genQE (genQE_FSI/build/genQE), tree genT: the scattered electron, the
    leading and recoil nucleons and the A-2 remnant nucleus. The remnant needs
    the mAm2 branch; use --no-remnant for files written before it existed.
    The remnant is skipped when it is not a bound nucleus (e.g. nn, pp).
  * genDIS, tree genT with pHadron: tagged DIS off an SRC pair. The electron, the recoil and the A-2
    remnant are written; the struck system is a single unhadronized four-vector and is left out.
  * genQE_FSI and genQE_3N_FSI, tree genT with the part_* list, and
    SRC_analysis_2N, tree events: the scattered electron plus the final-state
    particle list (post-FSI, or the unscattered nucleons with FSI off). GENIE
    drops the residual nucleus from that list, so no remnant is written.

The HepMC event number is the source tree entry, so simulated events can be
joined back to the generator tree, and the event weight is the generator
weight (nb per attempt). A JSON file next to the output records the source,
beams, cuts, counts and the cross section of the written sample.
"""
import argparse
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

M_E = 0.51099895e-3
M_P = 0.93827208816
M_N = 0.93956542052
M_N_GCF = 0.93892  # nucleon mass the GCF generators put on shell (particleMasses.hh)
M_U = 0.9314941024  # atomic mass unit [GeV], as in nuclearMasses.hh

# Target nuclei: (Z, A, nuclear mass [GeV]); masses as in nuclearMasses.hh
IONS = {
    "d": (1, 2, 2.01410178 * M_U - M_E),
    "He3": (2, 3, 3.0160293 * M_U - 2 * M_E),
    "He4": (2, 4, 4.00260325415 * M_U - 2 * M_E),
    "C12": (6, 12, 12.0 * M_U - 6 * M_E),
}
PDG_MASS = {11: M_E, 2212: M_P, 2112: M_N}


def ion_pdg(z, a):
    return 1000000000 + 10000 * z + 10 * a


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="generator ROOT file")
    ap.add_argument("output", nargs="?", help="output HepMC3 ASCII file")
    ap.add_argument("--print-ebeam-irf", action="store_true",
                    help="print the rest-frame electron energy to run the generator at, and exit")
    ap.add_argument("--beams", required=True,
                    help="electron energy x ion momentum per nucleon, in GeV, e.g. 5x41 or 10x110")
    ap.add_argument("--ion", required=True, choices=sorted(IONS), help="target nucleus of the generator run")
    ap.add_argument("--ebeam-irf", type=float,
                    help="rest-frame electron energy of the generator run [GeV]; needed only for "
                         "trees that do not store q (genQE_3N_FSI), and checked against --beams")
    ap.add_argument("--rtol", type=float, default=1e-3,
                    help="allowed relative mismatch between the run's and the beams' E_IRF (default 1e-3)")
    ap.add_argument("--paper-cuts", action="store_true", help="Q2 > 3 GeV^2 and xB > 1.2 (PRC 105, 034001)")
    ap.add_argument("--q2-min", type=float, help="keep Q2 > this [GeV^2]")
    ap.add_argument("--q2-max", type=float, help="keep Q2 < this [GeV^2]")
    ap.add_argument("--xb-min", type=float, help="keep xB > this (xB = Q2 / 2 m_p nu)")
    ap.add_argument("--no-remnant", action="store_true", help="do not write the A-2 remnant nucleus (genQE input)")
    ap.add_argument("--first", type=int, default=0, help="first tree entry to consider")
    ap.add_argument("-n", "--max-events", type=int, help="write at most this many events")
    args = ap.parse_args()
    if not args.print_ebeam_irf and (args.input is None or args.output is None):
        ap.error("input and output are required")
    if args.paper_cuts:
        args.q2_min = 3.0 if args.q2_min is None else args.q2_min
        args.xb_min = 1.2 if args.xb_min is None else args.xb_min
    return args


class Frame:
    """Rest frame (electron along +z) -> head-on collider frame (electron along -z, ion along +z)."""

    def __init__(self, e_beam, p_nucleon, ion):
        z, a, m_ion = IONS[ion]
        self.e_beam, self.p_nucleon, self.z, self.a, self.m_ion = e_beam, p_nucleon, z, a, m_ion
        self.bg = a * p_nucleon / m_ion
        self.g = math.hypot(1.0, self.bg)
        self.e_irf = self.g * (1.0 + self.bg / self.g) * e_beam

    def __call__(self, px, py, pz, e):
        """Transform rest-frame four-momenta (arrays or scalars)."""
        px, pz = -np.asarray(px), -np.asarray(pz)  # rotation by pi about y
        return px, np.asarray(py), self.g * pz + self.bg * e, self.g * e + self.bg * pz


def read_tree(path):
    f = uproot.open(path)
    keys = {k.split(";")[0] for k in f.keys()}
    if "genT" in keys:
        t = f["genT"]
    elif "events" in keys:
        t = f["events"]
    else:
        sys.exit(f"{path}: no genT or events tree")
    names = set(t.keys())
    if "part_pdg" in names:
        kind = "particles"
    elif {"pHadron", "pRec", "pe"} <= names:
        kind = "genDIS"
    elif {"pLead", "pRec", "pe"} <= names:
        kind = "genQE"
    else:
        sys.exit(f"{path}: tree {t.name} has neither part_pdg nor pLead/pRec/pe")
    seed = int(f["seed"].member("fVal")) if "seed" in keys else None
    attempts = None
    if "hAttempts" in keys:
        attempts = float(f["hAttempts"].values()[0])
    meta = Path(str(path) + ".meta.txt")
    if attempts is None and meta.exists():
        for line in meta.read_text().splitlines():
            if line.startswith("n_events_attempts="):
                attempts = float(line.split("=", 1)[1])
    return t, names, kind, attempts, seed


def electron_and_beam(t, names, sl):
    """Scattered electron (px, py, pz, E) per entry and the run's rest-frame beam energy (or None)."""
    if "electron" in names:  # SRC_analysis_2N: electron[4] and q[4] as (px, py, pz, E)
        k = t["electron"].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
        q = t["q"].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
        beam = k + q
    else:  # genT: pe[3] (massless) and, except genQE_3N_FSI, q[3]
        pe = t["pe"].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
        k = np.column_stack([pe, np.linalg.norm(pe, axis=1)])
        beam = None
        if "q" in names:
            q = t["q"].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
            b = pe + q
            beam = np.column_stack([b, np.linalg.norm(b, axis=1)])
    e_beam = None
    if beam is not None and len(beam):
        if np.max(np.abs(beam[:, :2])) > 1e-6 * np.max(beam[:, 3]):
            sys.exit("the run's beam (electron + q) is not along +z: not a GCF rest-frame tree?")
        e_beam = float(np.median(beam[:, 2]))
        if np.max(np.abs(beam[:, 2] / e_beam - 1)) > 1e-6:
            sys.exit("the beam energy (electron + q) varies between entries: radiative (-O) or Coulomb (-C) "
                     "runs are not supported")
    return k, e_beam


def kinematics(k, e_beam):
    nu = e_beam - k[:, 3]
    q2 = 2.0 * e_beam * (k[:, 3] - k[:, 2])  # massless electron along +z
    with np.errstate(divide="ignore", invalid="ignore"):
        xb = q2 / (2.0 * M_P * nu)
    return q2, xb


def remnant_species(z, a, lead, rec):
    """PDG of the A-2 remnant (a nucleon, or a nucleus), or None when it is not bound; and its (Z, A)."""
    z2 = z - int(lead == 2212) - int(rec == 2212)
    a2 = a - 2
    if a2 == 1:
        return (2212 if z2 == 1 else 2112), (z2, a2)
    if a2 >= 2 and 0 < z2 < a2:
        return ion_pdg(z2, a2), (z2, a2)
    return None, (z2, a2)


def fmt(*x):
    return " ".join(f"{v:.10e}" for v in x)


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "describe", "--always", "--dirty"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    args = parse_args()
    e_beam, p_nucleon = map(float, args.beams.lower().split("x"))
    frame = Frame(e_beam, p_nucleon, args.ion)
    if args.print_ebeam_irf:
        print(f"{frame.e_irf:.6f}")
        return
    t, names, kind, attempts, seed = read_tree(args.input)
    stop = t.num_entries
    sl = slice(args.first, stop)

    k, e_run = electron_and_beam(t, names, sl)
    if e_run is None:
        if args.ebeam_irf is None:
            sys.exit(f"{args.input}: the tree does not store q, so pass --ebeam-irf (the generator's beam energy)")
        e_run = args.ebeam_irf
    if abs(e_run / frame.e_irf - 1) > args.rtol:
        sys.exit(f"the generator ran at E_IRF = {e_run:.4f} GeV, but {args.beams} with {args.ion} needs "
                 f"{frame.e_irf:.4f} GeV (gamma = {frame.g:.4f})")

    w = t["weight"].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
    q2, xb = kinematics(k, e_run)
    sel = w > 0
    if args.q2_min is not None:
        sel &= q2 > args.q2_min
    if args.q2_max is not None:
        sel &= q2 < args.q2_max
    if args.xb_min is not None:
        sel &= xb > args.xb_min
    idx = np.flatnonzero(sel)[: args.max_events]

    # Final-state particles per selected entry, as (pdg, px, py, pz, E, m) in the rest frame.
    events = []
    balance = 0.0
    skipped_remnant = {}
    if kind == "genQE":
        a = {b: t[b].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
             for b in ("pLead", "pRec", "pAm2", "lead_type", "rec_type")}
        has_m = "mAm2" in names
        if not args.no_remnant and not has_m:
            sys.exit(f"{args.input} has no mAm2 branch (genQE older than this script): pass --no-remnant")
        m_am2 = t["mAm2"].array(library="np", entry_start=sl.start, entry_stop=sl.stop) if has_m else None
        z, a_num, m_a = IONS[args.ion]
        for i in idx:
            parts = [(11, *k[i, :3], math.hypot(np.linalg.norm(k[i, :3]), M_E), M_E)]
            for pdg, p in ((int(a["lead_type"][i]), a["pLead"][i]), (int(a["rec_type"][i]), a["pRec"][i])):
                m = PDG_MASS[pdg]
                parts.append((pdg, *p, math.sqrt(m * m + p @ p), m))
            if has_m:  # four-momentum balance with the generator's own masses
                e_tot = (k[i, 3] + math.sqrt(M_N_GCF**2 + a["pLead"][i] @ a["pLead"][i])
                         + math.sqrt(M_N_GCF**2 + a["pRec"][i] @ a["pRec"][i])
                         + math.sqrt(m_am2[i]**2 + a["pAm2"][i] @ a["pAm2"][i]))
                p_tot = k[i, :3] + a["pLead"][i] + a["pRec"][i] + a["pAm2"][i]
                d = np.append(p_tot - [0.0, 0.0, e_run], e_tot - e_run - m_a)
                balance = max(balance, float(np.max(np.abs(d))))
            if not args.no_remnant:
                pdg, (z2, a2) = remnant_species(z, a_num, int(a["lead_type"][i]), int(a["rec_type"][i]))
                if pdg is None:
                    if a2 > 0:
                        key = f"Z={z2},A={a2}"
                        skipped_remnant[key] = skipped_remnant.get(key, 0) + 1
                else:
                    p, m = a["pAm2"][i], float(m_am2[i])
                    parts.append((pdg, *p, math.sqrt(m * m + p @ p), m))
            events.append(parts)
    elif kind == "genDIS":
        # Tagged DIS off an SRC pair: the electron, the recoil and the A-2 remnant. The struck system
        # (pHadron, EHadron) is a single unhadronized four-vector, so it is not written; the recoil
        # tagging does not depend on it.
        a = {b: t[b].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
             for b in ("pRec", "pAm2", "EAm2", "lead_type", "rec_type")}
        z, a_num, m_a = IONS[args.ion]
        for i in idx:
            parts = [(11, *k[i, :3], math.hypot(np.linalg.norm(k[i, :3]), M_E), M_E)]
            pdg, p = int(a["rec_type"][i]), a["pRec"][i]
            m = PDG_MASS[pdg]
            parts.append((pdg, *p, math.sqrt(m * m + p @ p), m))
            if not args.no_remnant:
                rpdg, (z2, a2) = remnant_species(z, a_num, int(a["lead_type"][i]), pdg)
                if rpdg is None:
                    if a2 > 0:
                        key = f"Z={z2},A={a2}"
                        skipped_remnant[key] = skipped_remnant.get(key, 0) + 1
                else:
                    p2, e2 = a["pAm2"][i], float(a["EAm2"][i])
                    parts.append((rpdg, *p2, e2, math.sqrt(max(e2 * e2 - p2 @ p2, 0.0))))
            events.append(parts)
    else:
        b = {n: t[n].array(library="np", entry_start=sl.start, entry_stop=sl.stop)
             for n in ("part_pdg", "part_px", "part_py", "part_pz", "part_E")}
        for i in idx:
            parts = [(11, *k[i, :3], math.hypot(np.linalg.norm(k[i, :3]), M_E), M_E)]
            for pdg, px, py, pz, e in zip(*(b[n][i] for n in ("part_pdg", "part_px", "part_py", "part_pz", "part_E"))):
                m2 = e * e - (px * px + py * py + pz * pz)
                parts.append((int(pdg), px, py, pz, e, math.sqrt(max(m2, 0.0))))
            events.append(parts)

    # Beams: the rest-frame beams transformed with the same map as the final state.
    ebx, eby, ebz, ebe = frame(0.0, 0.0, e_run, math.hypot(e_run, M_E))
    ibx, iby, ibz, ibe = frame(0.0, 0.0, 0.0, frame.m_ion)
    pdg_ion = ion_pdg(frame.z, frame.a)
    beams = (f"P 1 0 11 {fmt(ebx, eby, ebz, ebe, M_E)} 4\n"
             f"P 2 0 {pdg_ion} {fmt(ibx, iby, ibz, ibe, frame.m_ion)} 4\n"
             "V -1 0 [1,2]\n")

    n_part = {}
    with open(args.output, "w") as out:
        out.write("HepMC::Version 3.02.02\nHepMC::Asciiv3-START_EVENT_LISTING\n")
        for i, parts in zip(idx, events):
            out.write(f"E {i + sl.start} 1 {len(parts) + 2}\nU GEV MM\nW {w[i]:.10e}\n{beams}")
            for j, (pdg, px, py, pz, e, m) in enumerate(parts):
                lx, ly, lz, le = frame(px, py, pz, e)
                out.write(f"P {j + 3} -1 {pdg} {fmt(lx, ly, lz, le, m)} 1\n")
                n_part[pdg] = n_part.get(pdg, 0) + 1
        out.write("HepMC::Asciiv3-END_EVENT_LISTING\n")

    info = {
        "input": str(Path(args.input).resolve()),
        "tree": t.name,
        "input_kind": kind,
        "output": str(Path(args.output).resolve()),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gcf_commit": git_commit(),
        "generator_seed": seed,
        "beams": {"electron_GeV": e_beam, "ion": args.ion, "ion_pdg": pdg_ion,
                  "ion_momentum_per_nucleon_GeV": p_nucleon, "ion_momentum_GeV": ibz,
                  "gamma": frame.g, "E_IRF_beams_GeV": frame.e_irf, "E_IRF_run_GeV": e_run},
        "frame": "head-on: electron along -z, ion along +z; no crossing angle or beam effects",
        "cuts": {"q2_min": args.q2_min, "q2_max": args.q2_max, "xb_min": args.xb_min,
                 "xb_definition": "Q2 / (2 m_p nu), m_p = 0.938272 GeV"},
        "entries_scanned": int(len(w)),
        "entries_passing": int(sel.sum()),
        "events_written": int(len(idx)),
        "first_entry": args.first,
        "attempts": attempts,
        "sum_weight_written": float(w[idx].sum()),
        "sigma_written_nb": float(w[idx].sum() / attempts) if attempts else None,
        "remnant": None if kind not in ("genQE", "genDIS") or args.no_remnant else "A-2 nucleus with generator mass (incl. E*)",
        "struck_system": "not written (genDIS hadronic system is unhadronized)" if kind == "genDIS" else None,
        "remnant_skipped_unbound": skipped_remnant,
        "particles_written_by_pdg": {str(p): n for p, n in sorted(n_part.items())},
        "max_rest_frame_4momentum_imbalance_GeV": balance if kind == "genQE" and "mAm2" in names else None,
    }
    if args.max_events is not None and len(idx) == args.max_events:
        info["note"] = "event count capped by --max-events: sigma_written_nb covers the written events only"
    Path(args.output + ".json").write_text(json.dumps(info, indent=2) + "\n")
    sigma = f"{info['sigma_written_nb']:.4g} nb" if info["sigma_written_nb"] is not None else "n/a (attempts unknown)"
    print(f"{args.input}: {sel.sum()} of {len(w)} entries pass, wrote {len(idx)} events to {args.output} "
          f"(E_IRF = {e_run:.3f} GeV, gamma = {frame.g:.3f}, sigma = {sigma})")
    if skipped_remnant:
        print(f"  unbound A-2 remnants not written: {skipped_remnant}")


if __name__ == "__main__":
    main()
