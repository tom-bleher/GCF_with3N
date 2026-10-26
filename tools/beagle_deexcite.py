#!/usr/bin/env python3
"""De-excite the A-2 remnant of GCF events with BeAGLE (FLUKA evaporation, fission, Fermi break-up).

    beagle_deexcite.py IN.hepmc3 OUT.hepmc3 [--seed S] [--jobs N] [--json OUT.json]

This is the evaporation stage of PRC 105, 034001 (doi:10.1103/PhysRevC.105.034001): GCF events
"processed through BeAGLE, run without any rescattering of the outgoing nucleons". It is a drop-in
alternative to tools/deexcite (Geant4) with the same input, output layout and JSON summary:

  * IN is gcf2hepmc.py output (head-on collider frame, beams as status 4), so the event selection
    and order are the converter's. Events without a remnant (A >= 3) are copied unchanged.
  * Each event with a remnant is boosted to the ion rest frame (a boost along z, which keeps the
    electron along -z as BeAGLE expects) and written as a BeAGLE MODEL GCF-FT record: the nine
    tracks DT_GCFEVNTQE (src/dpm_gcf.f) requires. The remnant's excitation is E* = m - M_gs(Z, A),
    with m its generated mass (GCF writes the A-2 mass including the sampled E*) and M_gs from the
    GCF mass table. The struck system is taken from energy-momentum conservation, so genDIS
    events (which do not carry it) work too.
  * BeAGLE runs at nuclear-reaction level 2: the remnant goes to FLUKA (DT_FLUKAIT), nothing is
    rescattered. FLUKA conserves the remnant's energy to the keV level BeAGLE writes, and its
    momentum to about 1e-4 of itself (up to ~0.2 MeV/c), where Geant4 is exact; both are far below
    the detector resolution. The products are boosted back and appended as status 1 after all other
    particles; the remnant becomes a status 3 entry with the products as daughters, exactly as
    tools/deexcite does (see its header for why status 3).

BeAGLE needs FLUKA's data tables in its working directory, so every job runs in its own directory
OUT.beagle/job_K with links to $FLUPRO/data. Large samples are split into contiguous blocks, one
BeAGLE process each (about 30 s of initialization per process); block K seeds both BeAGLE and
FLUKA's own generator (FSEED card) with S*1000 + K.

Environment (defaults in parentheses):
  BEAGLESYS  BeAGLE checkout with the BeAGLE executable (~/src/BeAGLE)
  FLUPRO     CERN FLUKA 4 installation (~/src/beagle-deps/fluka)
  LHAPATH    LHAPDF 5 PDF sets BeAGLE initializes (~/src/beagle-deps/lhapdf5/share/lhapdf/PDFsets)
Run it inside eic-shell, where BeAGLE was built (see BeAGLE's build_local.sh).
"""
import argparse
import json
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

M_E = 0.51099895e-3
M_P = 0.93827208816
M_U = 0.9314941024  # atomic mass unit [GeV], as in nuclearMasses.hh

# Ground-state nuclear masses [GeV] by (Z, A), as in src/include/constants/nuclearMasses.hh
M_GS = {
    (1, 3): 3.01604928199 * M_U - M_E,
    (2, 3): 3.0160293 * M_U - 2 * M_E,
    (2, 4): 4.00260325415 * M_U - 2 * M_E,
    (3, 6): 6.015122795 * M_U - 3 * M_E,
    (4, 8): 8.00530510 * M_U - 4 * M_E,
    (4, 9): 9.0121831 * M_U - 4 * M_E,
    (4, 10): 10.013534 * M_U - 4 * M_E,
    (5, 10): 10.0129370 * M_U - 5 * M_E,
    (5, 11): 11.0093054 * M_U - 5 * M_E,
    (6, 10): 10.016853 * M_U - 6 * M_E,
    (6, 11): 11.0114336 * M_U - 6 * M_E,
    (6, 12): 12.0 * M_U - 6 * M_E,
    (7, 14): 14.0030740048 * M_U - 7 * M_E,
    (8, 16): 15.99491461956 * M_U - 8 * M_E,
}
BEAGLE_ION = 80000  # DPMJET's code for a nucleus; A and Z are in the IDRES/IDXRES columns
MATCH_TOL = 2e-5    # GeV: BeAGLE writes momenta with 1 keV resolution


def ion_za(pdg):
    return (pdg // 10000) % 1000, (pdg // 10) % 1000


def ion_pdg(z, a):
    return 1000000000 + 10000 * z + 10 * a


def charge(pdg):
    if pdg > 1000000000:
        return ion_za(pdg)[0]
    return {2212: 1, 2112: 0, 11: -1, 22: 0}[pdg]


class Event:
    """One HepMC3 ASCII event as written by gcf2hepmc.py: raw lines plus the parsed particles."""

    def __init__(self, lines):
        self.lines = lines
        self.parts = []  # (line index, id, pdg, [px, py, pz, e], m, status)
        for k, line in enumerate(lines):
            if line.startswith("P "):
                f = line.split()
                self.parts.append((k, int(f[1]), int(f[3]), [float(x) for x in f[4:8]], float(f[8]), int(f[9])))


def read_hepmc(path):
    header, events, cur = [], [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("E "):
                cur = []
                events.append(cur)
            elif line.startswith("HepMC::Asciiv3-END"):
                cur = None
                continue
            if cur is None:
                if not events:
                    header.append(line)
                continue
            cur.append(line)
    return header, [Event(e) for e in events]


def boost_z(p, beta):
    """Boost (px, py, pz, E) along z with velocity beta."""
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    return [p[0], p[1], g * (p[2] + beta * p[3]), g * (p[3] + beta * p[2])]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def mass_signed(p):
    m2 = p[3] ** 2 - p[0] ** 2 - p[1] ** 2 - p[2] ** 2
    return math.copysign(math.sqrt(abs(m2)), m2)


def gcf_record(ievent, sel):
    """BeAGLE GCF-FT header and nine tracks for one event, in the ion rest frame."""
    e_beam, ion, ep, rec, rem, beta = sel["e_beam"], sel["ion"], sel["ep"], sel["rec"], sel["rem"], sel["beta"]
    irf = lambda p: boost_z(p, -beta)
    k, kp, pr, pa = irf(e_beam[3]), irf(ep[3]), irf(rec[3]), irf(rem[3])
    z_t, a_t = ion_za(ion[2])
    m_t = ion[4]
    z_r, a_r = ion_za(rem[2])
    # struck system after absorbing q, from energy-momentum conservation (the leading nucleon for genQE)
    lead = [k[0] - kp[0] - pr[0] - pa[0], k[1] - kp[1] - pr[1] - pa[1],
            k[2] - kp[2] - pr[2] - pa[2], k[3] + m_t - kp[3] - pr[3] - pa[3]]
    z_lead = z_t - charge(rec[2]) - z_r
    if z_lead not in (0, 1):
        raise ValueError(f"event {ievent}: struck system has charge {z_lead}")
    lt = 2212 if z_lead else 2112
    q = sub(k, kp)
    nu = q[3]
    q2 = q[0] ** 2 + q[1] ** 2 + q[2] ** 2 - nu * nu
    p1 = sub(lead, q)  # struck nucleon before absorbing q (off shell)
    # remnant with its generated mass (incl. E*); the boosted energy only carries rounding
    m_r = rem[4]
    pa[3] = math.sqrt(m_r * m_r + pa[0] ** 2 + pa[1] ** 2 + pa[2] ** 2)
    estar = m_r - M_GS[(z_r, a_r)]
    phi = math.atan2(kp[1], kp[0])
    zz = lambda pdg: int(pdg == 2212)
    head = (f"0 {ievent} 11 {a_t} {z_t} {k[2]:.12e} {rec[2]} {lt} {q2:.12e} {q2 / (2 * M_P * nu):.12e} "
            f"{nu:.12e} {phi:.12e} {p1[0]:.12e} {p1[1]:.12e} {p1[2]:.12e} {max(estar, 0.0):.12e} "
            f"{sel['w']:.12e} {sel['w']:.12e} 9")
    # I ISTHKK IDHKK JMOHKK(2) JMOHKK(1) JDAHKK(1) JDAHKK(2) px py pz E m A Z NOBAM
    tracks = [(12, lt, 0, 0, 7, 0, p1, 1, zz(lt)),
              (12, rec[2], 0, 0, 8, 0, pr, 1, zz(rec[2])),
              (12, BEAGLE_ION, 0, 0, 9, 0, pa, a_r, z_r),
              (21, 11, 0, 0, 5, 6, k, 0, -1),
              (21, 22, 0, 4, 7, 0, q, 0, 0),
              (1, 11, 0, 4, 0, 0, kp, 0, -1),
              (1, lt, 5, 1, 0, 0, lead, 1, zz(lt)),
              (1, rec[2], 0, 2, 0, 0, pr, 1, zz(rec[2])),
              (1000, BEAGLE_ION, 0, 3, 0, 0, pa, a_r, z_r)]
    lines = [head]
    for i, (st, pdg, mo2, mo1, da1, da2, p, a, z) in enumerate(tracks, 1):
        m = m_r if pdg == BEAGLE_ION else mass_signed(p)
        lines.append(f"{i} {st} {pdg} {mo2} {mo1} {da1} {da2} {p[0]:.12e} {p[1]:.12e} {p[2]:.12e} "
                     f"{p[3]:.12e} {m:.12e} {a} {z} 0")
    return lines, {"irf": [kp, lead, pr], "rem_irf": pa, "estar": estar, "za": (z_r, a_r)}


def select(ev):
    """Beams, e', recoil and remnant of an event, or None when it has no remnant to de-excite."""
    beams = [p for p in ev.parts if p[5] == 4]
    final = [p for p in ev.parts if p[5] == 1]
    rems = [p for p in final if p[2] > 1000000000 and ion_za(p[2])[1] >= 3]
    if not rems:
        return None
    if len(rems) > 1:
        raise ValueError("more than one remnant nucleus in an event")
    e_beam = next(p for p in beams if p[2] == 11)
    ion = next(p for p in beams if p[2] > 1000000000)
    if abs(ion[3][0]) > 1e-9 * ion[3][3] or abs(ion[3][1]) > 1e-9 * ion[3][3]:
        raise ValueError("the ion beam is not along z: run before the afterburner, on gcf2hepmc.py output")
    ep = next(p for p in final if p[2] == 11)
    rem = rems[0]
    nucleons = [p for p in final if p[2] in (2212, 2112) and p[0] < rem[0]]
    w = next((float(l.split()[1]) for l in ev.lines if l.startswith("W ")), 1.0)
    return {"e_beam": e_beam, "ion": ion, "ep": ep, "rec": nucleons[-1], "rem": rem,
            "beta": ion[3][2] / ion[3][3], "w": w}


def card(code, what=(), sdum=""):
    """BeAGLE/DPMJET input card: A10, up to six numbers in columns 11-70 (read list-directed), A8."""
    nums = "".join(f" {int(x)}" if float(x).is_integer() else f" {x}" for x in what)
    if len(nums) > 60:
        raise ValueError(f"{code} card: numbers do not fit in columns 11-70")
    return f"{code:<10}{nums:<60}{sdum}".rstrip()


def run_job(job, todo, workdir, env, seed, z_t, a_t, e_e, p_n):
    """Write one BeAGLE input block, run BeAGLE on it, return its output path."""
    workdir.mkdir(parents=True, exist_ok=True)
    for f in Path(env["FLUPRO"], "data").iterdir():
        link = workdir / f.name
        if not link.exists():
            link.symlink_to(f)
    with open(workdir / "gcf.txt", "w") as f:
        f.write(f"GCF events for BeAGLE, written by beagle_deexcite.py\n{len(todo)} {a_t} {z_t}\n")
        f.write("\n".join(["header"] * 5) + "\n")
        for lines in todo:
            f.write("\n".join(lines) + "\n")
    steer = "\n".join([
        card("PROJPAR", sdum="ELECTRON"),
        card("TARPAR", (a_t, z_t, 0.0)),
        card("MOMENTUM", (e_e, p_n)),
        # GCF-FT: seed, nuclear-reaction level 2 (hand the (A-2)* to FLUKA, no rescattering)
        card("MODEL", (seed, 2.0, 0.0, 0.0, 0.0), "GCF-FT"),
        # FLUKA's own generator (FLRNDM), which does the break-up: RNINIT(13, seed, 3, 3), BeAGLE's
        # tested values with the sequence seed set; without it every job draws FLUKA's default sequence
        card("FSEED", (13, seed, 3, 3)),
        f"{'INPFILE':<10}gcf.txt",
        f"{'OUTFILE':<10}beagle.txt",
        card("START", (0.0, 0.0)),
        "STOP"]) + "\n"
    (workdir / "steer.inp").write_text(steer)
    for old in workdir.glob("fort.*"):
        old.unlink()
    with open(workdir / "steer.inp") as fin, open(workdir / "beagle.log", "w") as log:
        rc = subprocess.run([str(Path(env["BEAGLESYS"]) / "BeAGLE")], stdin=fin, stdout=log, stderr=subprocess.STDOUT,
                            cwd=workdir, env=env).returncode
    if rc != 0:
        raise RuntimeError(f"BeAGLE job {job} failed (exit {rc}); see {workdir}/beagle.log")
    return workdir / "beagle.txt"


def read_beagle(path):
    """Final-state particles per event of a BeAGLE GCF output file, in order."""
    events, cur = [], None
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) > 20 and t[0] == "0":
                cur = []
                events.append(cur)
            elif cur is not None and len(t) == 18 and t[0].isdigit() and t[1] == "1":
                pdg, a, z = int(t[2]), int(t[15]), int(t[16])
                if pdg == BEAGLE_ION:
                    pdg = (2212 if z else 2112) if a == 1 else ion_pdg(z, a)
                cur.append((pdg, [float(x) for x in t[7:11]], float(t[11])))
    return events


def git_describe(path):
    try:
        return subprocess.run(["git", "-C", str(path), "describe", "--always", "--dirty"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="gcf2hepmc.py output (HepMC3 ASCII)")
    ap.add_argument("output", help="output HepMC3 ASCII file")
    ap.add_argument("--seed", type=int, default=1, help="seed; BeAGLE job K uses seed*1000 + K (default 1)")
    ap.add_argument("--jobs", type=int, default=8, help="parallel BeAGLE processes (default 8)")
    ap.add_argument("--min-block", type=int, default=2000, help="smallest number of events per BeAGLE process")
    ap.add_argument("--json", help="summary file (default OUTPUT.json)")
    ap.add_argument("--keep", action="store_true", help="keep BeAGLE's input and output files")
    args = ap.parse_args()

    home = Path.home()
    env = dict(os.environ)
    env.setdefault("BEAGLESYS", str(home / "src/BeAGLE"))
    env.setdefault("FLUPRO", str(home / "src/beagle-deps/fluka"))
    env.setdefault("LHAPATH", str(home / "src/beagle-deps/lhapdf5/share/lhapdf/PDFsets"))
    for key, need in (("BEAGLESYS", "BeAGLE"), ("FLUPRO", "data/nuclear.bin"), ("LHAPATH", "cteq6ll.LHpdf")):
        if not Path(env[key], need).exists():
            sys.exit(f"beagle_deexcite: {key}={env[key]} has no {need}")
    # job K seeds BeAGLE's RANLUX and FLUKA's FLRNDM (RNINIT, below 9e8) with seed*1000 + K;
    # a BeAGLE seed of 0 would mean a time-based seed
    if not 1 <= args.seed or args.seed * 1000 + args.jobs >= 9 * 10**8:
        sys.exit("beagle_deexcite: --seed must be between 1 and 899999")

    header, events = read_hepmc(args.input)
    if not events:
        sys.exit(f"beagle_deexcite: no events in {args.input}")
    sels = [select(ev) for ev in events]
    todo = [i for i, s in enumerate(sels) if s]
    beams = sels[todo[0]] if todo else None

    out_events = {}
    stats = Counter()
    products = Counter()
    mult = {"n": Counter(), "p": Counter()}
    estars = []
    max_de_rest = max_dp_rest = max_dp = max_negative = 0.0
    n_jobs = 0
    if todo:
        z_t, a_t = ion_za(beams["ion"][2])
        e_e = -beams["e_beam"][3][2]
        p_n = beams["ion"][3][2] / a_t
        n_jobs = max(1, min(args.jobs, len(todo) // args.min_block))
        blocks = [todo[k * len(todo) // n_jobs:(k + 1) * len(todo) // n_jobs] for k in range(n_jobs)]
        recs = {}
        inputs = []
        for block in blocks:
            lines = []
            for j, i in enumerate(block, 1):
                rec_lines, recs[i] = gcf_record(j, sels[i])
                lines.append(rec_lines)
            inputs.append(lines)
        work = Path(args.output + ".beagle")
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            outs = list(pool.map(lambda k: run_job(k, inputs[k], work / f"job_{k}", env, args.seed * 1000 + k,
                                                   z_t, a_t, e_e, p_n), range(n_jobs)))
        for block, path in zip(blocks, outs):
            bev = read_beagle(path)
            if len(bev) != len(block):
                sys.exit(f"beagle_deexcite: {path} has {len(bev)} events, expected {len(block)}")
            for i, final in zip(block, bev):
                r, s = recs[i], sels[i]
                # BeAGLE copies e', the struck system and the recoil through unchanged; the rest is the remnant's
                rest = list(final)
                for want in r["irf"]:
                    hit = [x for x in rest if max(abs(a - b) for a, b in zip(x[1][:3], want[:3])) < MATCH_TOL]
                    if len(hit) != 1:
                        sys.exit(f"beagle_deexcite: event {events[i].lines[0]}: e'/lead/recoil not passed through")
                    rest.remove(hit[0])
                z_r, a_r = r["za"]
                sum_a = sum(1 if p in (2212, 2112) else ion_za(p)[1] if p > 1000000000 else 0 for p, _, _ in rest)
                sum_z = sum(charge(p) for p, _, _ in rest)
                if (sum_a, sum_z) != (a_r, z_r):
                    stats["violations"] += 1
                tot = [sum(x[1][c] for x in rest) for c in range(4)]
                max_de_rest = max(max_de_rest, abs(tot[3] - r["rem_irf"][3]))
                max_dp_rest = max(max_dp_rest, *(abs(a - b) for a, b in zip(tot[:3], r["rem_irf"][:3])))
                lab = [(pdg, boost_z(p, s["beta"]), m) for pdg, p, m in rest]
                tot_lab = [sum(x[1][c] for x in lab) for c in range(4)]
                max_dp = max(max_dp, *(abs(a - b) for a, b in zip(tot_lab, s["rem"][3])))
                if r["estar"] < 0:
                    stats["negative"] += 1
                    max_negative = max(max_negative, -r["estar"])
                estars.append(max(r["estar"], 0.0))
                for pdg, _, _ in lab:
                    products[pdg] += 1
                mult["n"][sum(p == 2112 for p, _, _ in lab)] += 1
                mult["p"][sum(p == 2212 for p, _, _ in lab)] += 1
                out_events[i] = lab
        if not args.keep:
            for k in range(n_jobs):
                for f in ("gcf.txt", "beagle.txt"):
                    (work / f"job_{k}" / f).unlink(missing_ok=True)
    for i, s in enumerate(sels):
        if not s:
            mult["n"][0] += 1
            mult["p"][0] += 1

    with open(args.output, "w") as out:
        out.write("\n".join(header) + "\n")
        for i, ev in enumerate(events):
            lines = list(ev.lines)
            if i in out_events:
                rem, prods = sels[i]["rem"], out_events[i]
                f = lines[rem[0]].split()
                f[9] = "3"  # history only; the products are the final state (as tools/deexcite)
                lines[rem[0]] = " ".join(f)
                e = lines[0].split()
                n_vtx, n_part = int(e[2]), int(e[3])
                vid = -(n_vtx + 1)
                e[2], e[3] = str(n_vtx + 1), str(n_part + len(prods))
                lines[0] = " ".join(e)
                lines.append(f"V {vid} 0 [{rem[1]}]")
                for j, (pdg, p, m) in enumerate(prods, n_part + 1):
                    lines.append(f"P {j} {vid} {pdg} {p[0]:.10e} {p[1]:.10e} {p[2]:.10e} {p[3]:.10e} {m:.10e} 1")
            out.write("\n".join(lines) + "\n")
        out.write("HepMC::Asciiv3-END_EVENT_LISTING\n")

    n_rem = len(out_events)
    info = {
        "input": str(Path(args.input).resolve()),
        "output": str(Path(args.output).resolve()),
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tool": "BeAGLE MODEL GCF-FT, nuclear-reaction level 2 (FLUKA de-excitation, no rescattering)",
        "beagle": {"dir": env["BEAGLESYS"], "commit": git_describe(env["BEAGLESYS"])},
        "fluka": {"dir": env["FLUPRO"],
                  "version": (Path(env["FLUPRO"], "Version.tag").read_text().split("\n")[0]
                              if Path(env["FLUPRO"], "Version.tag").exists() else None)},
        "seed": args.seed,
        "beagle_jobs": n_jobs,
        "events": len(events),
        "remnants": n_rem,
        "mean_Estar_GeV": sum(estars) / n_rem if n_rem else 0.0,
        "max_Estar_GeV": max(estars, default=0.0),
        "remnants_below_ground_state": stats["negative"],
        "max_below_ground_state_GeV": max_negative,
        "max_4momentum_imbalance_rest_frame_MeV": 1e3 * max(max_de_rest, max_dp_rest),
        "max_energy_imbalance_ion_rest_frame_MeV": 1e3 * max_de_rest,
        "max_momentum_imbalance_ion_rest_frame_MeV": 1e3 * max_dp_rest,
        "rest_frame": "ion rest frame (BeAGLE's frame), not the remnant's",
        "max_4momentum_imbalance_lab_GeV": max_dp,
        "A_or_Z_not_conserved": stats["violations"],
        "products_by_pdg": {str(p): n for p, n in sorted(products.items())},
        "events_by_multiplicity": {k: {str(m): n for m, n in sorted(v.items())} for k, v in mult.items()},
    }
    Path(args.json or args.output + ".json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"deexcite: {len(events)} events, {n_rem} remnants, mean E* {1e3 * info['mean_Estar_GeV']:.3f} MeV, "
          f"max imbalance {1e3 * max_de_rest:.4f} MeV in energy, {1e3 * max_dp_rest:.4f} MeV in momentum (ion rest frame), "
          f"A/Z not conserved in {stats['violations']}; {stats['negative']} below the ground state "
          f"(max {1e3 * max_negative:.3f} MeV) [BeAGLE/FLUKA]")
    # FLUKA conserves the remnant's energy to the keV level BeAGLE writes, but its break-up
    # conserves momentum only to about 1e-4 of the remnant's momentum (up to ~0.2 MeV/c here)
    if stats["violations"] or max_de_rest > 5e-5 or max_dp_rest > 1e-3:
        print("beagle_deexcite: conservation check failed", file=sys.stderr)
        sys.exit(4)


if __name__ == "__main__":
    main()
