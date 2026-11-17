#!/usr/bin/env python3
"""Check the Scale-and-Scheme contacts in gcfSRC.cc against the published table.

Parses set_Contacts_SS_k() and set_Contacts_SS_r() in
src/libraries/nucleus/gcfSRC.cc and compares every C_pn^{s=1} and C_pp^{s=0}
value and uncertainty with data/cruztorres2021_contacts.csv, the numbers of
Cruz-Torres et al., Nature Phys. 17, 306 (2021), Source Data MOESM6.

It fails (exit 1) when
  - a published value or uncertainty differs from the code,
  - a published cell has no code entry (except the known gaps in KNOWN_GAPS),
  - the code has a non-zero Cpn1/Cpp0/Cnn0 the table does not publish,
  - a code block does not assign all eight contacts exactly once, or assigns
    one without the "* ( A * 0.5)" factor,
  - Cnn0 != Cpp0 for an N = Z nucleus.
Cpn0 (C_pn^{s=0}) is not in the published Source Data and is not checked.

  check_contacts.py [--src gcfSRC.cc] [--table cruztorres2021_contacts.csv]
"""

import argparse
import csv
import re
import sys
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SRC = HERE.parents[1] / "src" / "libraries" / "nucleus" / "gcfSRC.cc"
DEFAULT_TABLE = HERE / "data" / "cruztorres2021_contacts.csv"

NUCLEI = {"2H": (1, 1), "3H": (1, 2), "3He": (2, 1), "4He": (2, 2), "6Li": (3, 3),
          "12C": (6, 6), "16O": (8, 8), "40Ca": (20, 20)}
NAMES = {zn: n for n, zn in NUCLEI.items()}
FIELDS = ["Cpp0", "d_Cpp0", "Cnn0", "d_Cnn0", "Cpn0", "d_Cpn0", "Cpn1", "d_Cpn1"]

# Published cells with no code set: the paper gives NV2+3-Ia* pp0 for 6Li but no
# pn1, so gcfSRC.cc has no NV2_1a 6Li block.
KNOWN_GAPS = {("pp0", "r", "NV2_1a", "6Li"), ("pp0", "k", "NV2_1a", "6Li")}

FUNC_RE = re.compile(r"bool\s+gcfSRC::set_Contacts_(\w+)\s*\(")
NUC_RE = re.compile(r"\(\s*Z\s*==\s*(\d+)\s*\)\s*&&\s*\(\s*N\s*==\s*(\d+)\s*\)")
U_RE = re.compile(r"\bu\s*==\s*(\w+)\s*\)")
ASSIGN_RE = re.compile(r"^\s*(d_C(?:pp0|nn0|pn0|pn1)|C(?:pp0|nn0|pn0|pn1))\s*=")
VALUE_RE = re.compile(r"^\s*(\w+)\s*=\s*([0-9]*\.?[0-9]+)\s*\*\s*\(\s*A\s*\*\s*0\.5\s*\)\s*;")


def parse_code(path, errors):
    """Return {(space, Z, N, u): {field: (Decimal, line)}} for the SS sets."""
    blocks = {}
    space = zn = u = None
    for i, line in enumerate(Path(path).read_text().splitlines(), 1):
        m = FUNC_RE.search(line)
        if m:
            name = m.group(1)
            space = {"SS_k": "k", "SS_r": "r"}.get(name)
            zn = u = None
            continue
        if space is None:
            continue
        m = NUC_RE.search(line)
        if m:  # the 40Ca branch also matches 40Ar; the first pair names it
            zn, u = (int(m.group(1)), int(m.group(2))), None
            continue
        m = U_RE.search(line)
        if m:
            u = m.group(1)
            continue
        if not ASSIGN_RE.match(line):
            continue
        m = VALUE_RE.match(line)
        if m is None or zn is None or u is None:
            errors.append(f"gcfSRC.cc:{i}: cannot parse contact assignment: {line.strip()}")
            continue
        block = blocks.setdefault((space, zn[0], zn[1], u), {})
        if m.group(1) in block:
            errors.append(f"gcfSRC.cc:{i}: {m.group(1)} assigned twice in {space}-space "
                          f"{NAMES.get(zn, zn)} {u}")
        block[m.group(1)] = (Decimal(m.group(2)), i)
    return blocks


def read_table(path):
    rows = (l for l in Path(path).read_text().splitlines() if l and not l.startswith("#"))
    return list(csv.DictReader(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--src", default=DEFAULT_SRC, help="gcfSRC.cc to check")
    ap.add_argument("--table", default=DEFAULT_TABLE, help="published contacts (CSV)")
    args = ap.parse_args()

    errors = []
    blocks = parse_code(args.src, errors)
    table = read_table(args.table)

    for key, block in sorted(blocks.items()):
        missing = [f for f in FIELDS if f not in block]
        if missing:
            sp, z, n, u = key
            errors.append(f"{sp}-space {NAMES.get((z, n), (z, n))} {u}: "
                          f"no assignment of {', '.join(missing)}")

    published = set()
    nok = 0
    for r in table:
        cell = (r["contact"], r["space"], r["interaction"], r["nucleus"])
        z, n = NUCLEI[r["nucleus"]]
        key = (r["space"], z, n, r["interaction"])
        published.add((key, r["contact"]))
        if key not in blocks:
            if cell in KNOWN_GAPS:
                print(f"note: published {' '.join(cell)} = {r['value']} +- {r['unc']} "
                      "has no code set (known gap)")
            else:
                errors.append(f"published {' '.join(cell)} = {r['value']} +- {r['unc']}: "
                              "no code set")
            continue
        field = "Cpn1" if r["contact"] == "pn1" else ("Cnn0" if z == 1 else "Cpp0")
        for f, want in ((field, r["value"]), ("d_" + field, r["unc"])):
            if f not in blocks[key]:
                continue  # already reported as a missing assignment
            have, line = blocks[key][f]
            if have == Decimal(want):
                nok += 1
            else:
                errors.append(f"gcfSRC.cc:{line}: {r['space']}-space {r['nucleus']} "
                              f"{r['interaction']} {f} = {have}, published {want}")

    for key, block in sorted(blocks.items()):
        sp, z, n, u = key
        nuc = NAMES.get((z, n), f"Z={z},N={n}")
        for contact, fields in (("pn1", ["Cpn1"]), ("pp0", ["Cpp0", "Cnn0"])):
            if (key, contact) in published:
                continue
            for f in fields:
                if f in block and block[f][0] != 0:
                    errors.append(f"gcfSRC.cc:{block[f][1]}: {sp}-space {nuc} {u} {f} = "
                                  f"{block[f][0]} is not in the published table")
        if z == n:
            for f in ("Cpp0", "d_Cpp0"):
                g = f.replace("pp", "nn")
                if f in block and g in block and block[f][0] != block[g][0]:
                    errors.append(f"gcfSRC.cc:{block[g][1]}: {sp}-space {nuc} {u} "
                                  f"{g} = {block[g][0]} != {f} = {block[f][0]} (N = Z)")

    for e in errors:
        print("FAIL:", e)
    print(f"{nok} of {2 * len(table)} published values and uncertainties match "
          f"({len(blocks)} code sets parsed); {len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
