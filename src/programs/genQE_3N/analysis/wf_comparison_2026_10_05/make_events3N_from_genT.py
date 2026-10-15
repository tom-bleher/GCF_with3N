#!/usr/bin/env python3
"""genT -> 3N-style `events` tree (lead / recoil2 / recoil3), the format
plotting/plot_fig4_xB_by_config.py and the other 3N figure scripts read.

The repo's plotting/make_events_from_genT.py writes the detection-based 2N-style
tree (lead + a single recoil), which those 3N scripts cannot use. This keeps the
three generated nucleons: lead = p1_gen (outgoing struck nucleon, as written by
the generator), recoil2 = p2_gen, recoil3 = p3_gen. The figure scripts form the
initial struck momentum themselves as lead - q.

Handles both schemas: current (p1_gen/p2_gen/p3_gen) and the pre-replacement
code (pLead_pre/p2_pre/p3_pre).

Usage: make_events3N_from_genT.py --input genT.root --output events.root [--ebeam 6.0]
"""
import argparse
import numpy as np
import uproot

MN = 0.93892


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--ebeam", type=float, default=6.0)
    a = ap.parse_args()

    f = uproot.open(a.input)
    t = f["genT"]
    have = set(t.keys())
    gen = ("pLead_pre", "p2_pre", "p3_pre") if "p1_gen" not in have else ("p1_gen", "p2_gen", "p3_gen")
    arr = t.arrays(["weight", "pe", *gen, "N1_type", "N2_type", "N3_type"], library="np")

    w = arr["weight"]
    pe = arr["pe"]
    lead3, r2_3, r3_3 = arr[gen[0]], arr[gen[1]], arr[gen[2]]
    E = a.ebeam
    q3 = np.column_stack([-pe[:, 0], -pe[:, 1], E - pe[:, 2]])
    pe_mag = np.linalg.norm(pe, axis=1)
    nu = E - pe_mag
    Q2 = np.sum(q3**2, axis=1) - nu**2
    xB = Q2 / (2 * MN * np.maximum(nu, 1e-12))
    theta_e = np.degrees(np.arccos(np.clip(pe[:, 2] / np.maximum(pe_mag, 1e-12), -1, 1)))
    pmiss = np.linalg.norm(lead3 - q3, axis=1)

    def p4(p3v, m=MN):
        return np.column_stack([p3v, np.sqrt(np.sum(p3v**2, axis=1) + m*m)])

    out = {
        "weight": w,
        "N1_type": arr["N1_type"].astype(np.int32),
        "N2_type": arr["N2_type"].astype(np.int32),
        "N3_type": arr["N3_type"].astype(np.int32),
        "electron": np.column_stack([pe, pe_mag]),
        "lead": p4(lead3), "recoil2": p4(r2_3), "recoil3": p4(r3_3),
        "q": np.column_stack([q3, nu]),
        "Q2": Q2, "xB": xB, "nu": nu, "pmiss": pmiss,
        "scattering_angle": theta_e,
    }
    with uproot.recreate(a.output) as fo:
        fo.mktree("events", {
            "weight": "f8", "N1_type": "i4", "N2_type": "i4", "N3_type": "i4",
            "electron": ("f8", (4,)), "lead": ("f8", (4,)),
            "recoil2": ("f8", (4,)), "recoil3": ("f8", (4,)), "q": ("f8", (4,)),
            "Q2": "f8", "xB": "f8", "nu": "f8", "pmiss": "f8", "scattering_angle": "f8",
        })
        fo["events"].extend(out)
    print(f"wrote {a.output}: {len(w)} events")


if __name__ == "__main__":
    main()
