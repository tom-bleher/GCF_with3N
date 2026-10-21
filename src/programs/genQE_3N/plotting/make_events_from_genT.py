#!/usr/bin/env python3
"""genT -> events converter: detection-based (e,e'pp) event building.

Reads a generator tree (``genT``) whose post-FSI final state is a flat
particle list (nPart, part_pdg, part_px/py/pz/E) and writes the flat
``events`` tree the paper figure scripts read.

DETECTION-BASED (e,e'pp) EVENT BUILDING — identical rule to
SRC_analysis_2N.cpp:
  * detected proton = final-state particle with pdg 2212 and |p| > kF
  * lead   = the detected proton closest in angle to q
  * recoil = the highest-momentum other detected proton
  * pmiss  = |p_lead - q| (from the DETECTED lead)
The paper's lead-identification cuts (theta_pq < 8 deg, p_lead > 0.75 |q|),
the (e,e'pp) requirement nProtonAboveKF == 2, and all other kinematic cuts
are applied downstream in the figure scripts.

Works for any genT with the flat-list branches (3N and 2N generators alike).

Usage:
    .venv/bin/python plotting/make_events_from_genT.py \
        --input events/3N_FSI_....root --output events_3N.root \
        --ebeam 6.0 --kF 0.22
"""
import argparse

import numpy as np
import uproot

MN = 0.93892  # GeV
P_CODE = 2212
N_CODE = 2112


GENT_BRANCHES = ["weight", "pe", "part_pdg", "part_px", "part_py", "part_pz",
                 "part_E"]


def build_events(arrays, ebeam, kf):
    """Apply the detection-based (e,e'pp) event building to one chunk of genT
    arrays (GENT_BRANCHES, awkward) and return the `events` columns as a dict
    of numpy arrays."""
    import awkward as ak

    w = np.asarray(arrays["weight"])
    pe = np.asarray(arrays["pe"])           # (N,3) scattered electron
    pdg = arrays["part_pdg"]
    px, py, pz = arrays["part_px"], arrays["part_py"], arrays["part_pz"]
    pE = arrays["part_E"]

    n = len(w)
    Ee = np.linalg.norm(pe, axis=1)
    electron = np.column_stack([pe, Ee])

    qx = -pe[:, 0]
    qy = -pe[:, 1]
    qz = ebeam - pe[:, 2]
    qE = ebeam - Ee
    qvec = np.column_stack([qx, qy, qz])
    qmag = np.linalg.norm(qvec, axis=1)

    Q2 = qmag**2 - qE**2
    nu = qE
    with np.errstate(divide="ignore", invalid="ignore"):
        xB = np.where(nu > 0, Q2 / (2 * MN * nu), -1.0)
    scattering_angle = np.degrees(
        np.arccos(np.clip(pe[:, 2] / (Ee + 1e-30), -1, 1)))

    # --- detection-based selection over the particle list ---
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    is_p_det = (pdg == P_CODE) & (pmag > kf)
    is_n_det = (pdg == N_CODE) & (pmag > kf)
    nProton = ak.to_numpy(ak.sum(is_p_det, axis=1)).astype(np.int32)
    nNucleon = (nProton
                + ak.to_numpy(ak.sum(is_n_det, axis=1))).astype(np.int32)

    # angle of each detected proton to q; +inf for everything else
    dot_q = (px * qx[:, None] + py * qy[:, None] + pz * qz[:, None])
    cosang = dot_q / (pmag * qmag[:, None] + 1e-30)
    # np.clip is not a ufunc, so it chokes on jagged arrays; min/max are.
    ang = np.arccos(np.minimum(np.maximum(cosang, -1.0), 1.0))
    ang_masked = ak.where(is_p_det, ang, np.inf)
    lead_idx = ak.argmin(ang_masked, axis=1, keepdims=True)

    def gather(field, idx, default=0.0):
        return ak.to_numpy(ak.fill_none(ak.firsts(field[idx]), default))

    lead = np.column_stack([gather(px, lead_idx), gather(py, lead_idx),
                            gather(pz, lead_idx), gather(pE, lead_idx)])
    has_lead = ak.to_numpy(
        ak.fill_none(ak.firsts(ang_masked[lead_idx]), np.inf)) < np.inf
    lead[~has_lead] = 0.0

    # recoil: highest-|p| detected proton that is not the lead
    li = ak.to_numpy(ak.fill_none(ak.firsts(lead_idx), -1))
    local_idx = ak.local_index(pmag, axis=1)
    is_rec_cand = is_p_det & (local_idx != li[:, None])
    p_masked = ak.where(is_rec_cand, pmag, -1.0)
    rec_idx = ak.argmax(p_masked, axis=1, keepdims=True)
    recoil = np.column_stack([gather(px, rec_idx), gather(py, rec_idx),
                              gather(pz, rec_idx), gather(pE, rec_idx)])
    has_rec = ak.to_numpy(
        ak.fill_none(ak.firsts(p_masked[rec_idx]), -1.0)) > 0.0
    recoil[~has_rec] = 0.0

    p1 = lead[:, :3] - qvec
    pmiss = np.where(has_lead, np.linalg.norm(p1, axis=1), -1.0)

    out = {
        "weight": w,
        "electron": electron,
        "q": np.column_stack([qvec, qE]),
        "lead": lead,
        "recoil": recoil,
        "Q2": Q2,
        "xB": xB,
        "nu": nu,
        "pmiss": pmiss,
        "scattering_angle": scattering_angle,
        "nProtonAboveKF": nProton,
        "nNucleonAboveKF": nNucleon,
    }
    return out


def convert(input_path, output_path, ebeam, kf, chunk="200 MB"):
    fin = uproot.open(input_path)
    tree = fin["genT"]
    n_total = tree.num_entries
    print(f"Input : {input_path}  ({n_total} events)")
    print(f"Output: {output_path}  (Ebeam={ebeam} GeV, kF={kf} GeV/c)")

    fout = uproot.recreate(output_path)
    first = True

    for arrays in tree.iterate(GENT_BRANCHES, library="ak", step_size=chunk):
        out = build_events(arrays, ebeam, kf)
        if first:
            # mktree, not plain assignment: assignment writes an RNTuple in
            # recent uproot, which the TTree-based figure scripts cannot read.
            fout.mktree("events", {
                "weight": "f8",
                "electron": ("f8", (4,)),
                "q": ("f8", (4,)),
                "lead": ("f8", (4,)),
                "recoil": ("f8", (4,)),
                "Q2": "f8", "xB": "f8", "nu": "f8", "pmiss": "f8",
                "scattering_angle": "f8",
                "nProtonAboveKF": "i4", "nNucleonAboveKF": "i4",
            })
            first = False
        fout["events"].extend(out)

    # carry the attempts histogram through for per-attempt normalization
    if "hAttempts" in fin:
        h = fin["hAttempts"]
        fout["hAttempts"] = h.to_numpy()
        print(f"hAttempts: {h.values()[0]:.0f} generation attempts")
    fout.close()
    print("Done.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="genT ROOT file")
    ap.add_argument("--output", required=True, help="events ROOT file")
    ap.add_argument("--ebeam", type=float, default=6.0, help="beam energy (GeV)")
    ap.add_argument("--kF", type=float, default=0.22,
                    help="detection threshold |p| > kF (GeV/c)")
    args = ap.parse_args()
    convert(args.input, args.output, args.ebeam, args.kF)


if __name__ == "__main__":
    main()
