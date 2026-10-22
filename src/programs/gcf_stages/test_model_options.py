#!/usr/bin/env python3
"""Every gcf_generate model accepts its own options and rejects the others'.

For each model-specific option in tools/gcf/schema.py (an Opt whose `models`
is not all of them) and each model, runs `gcf_generate --check` with that
option and checks that it is accepted (exit 0) when the schema lists the model
and rejected with an "applies to" error (exit 1) when it does not, and that
the options for every model (--cs, --q2) are accepted by each. Before that it
compares the binary's option table (kModelOptions in gcf_generate.cc) with the
schema directly. Together this keeps the binary and the front end's schema in
agreement.

    python3 src/programs/gcf_stages/test_model_options.py [path/to/gcf_generate]

Needs the environment gcf_generate runs in (ROOT, GENIE libraries and
GXMLPATH, e.g. from build/gcf_env.sh). Default binary: build/gcf_generate
next to this file.
"""

import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools", "gcf"))
import schema  # noqa: E402

# Beam and target per model for which --check finds an accepted event.
BASE = {
    "3N": ["--beam", "6", "-A", "12", "-Z", "6"],
    "2N": ["--beam", "5.01", "-A", "12", "-Z", "6"],
    "MF": ["--beam", "5.01", "-A", "12", "-Z", "6"],
}


def table_mismatches():
    """kModelOptions in gcf_generate.cc against the model-specific Opts of the
    schema: the same flags with the same models (no binary needed)."""
    src = open(os.path.join(HERE, "gcf_generate.cc")).read()
    body = re.search(r"kModelOptions\[\]\s*=\s*\{(.*?)\n\};", src, re.S)
    if not body:
        return ["kModelOptions not found in gcf_generate.cc"]
    binary = {flag: set(models.split("/"))
              for flag, models in re.findall(r'\{"(--[\w-]+)",\s*"([\w/]+)"', body.group(1))}
    schema_tab = {o.flag: set(o.models) for o in schema.GENERATE
                  if o.flag is not None and set(o.models) != set(schema.MODELS)}
    out = []
    for flag in sorted(set(binary) | set(schema_tab)):
        b, s = binary.get(flag), schema_tab.get(flag)
        if b != s:
            out.append(f"{flag}: gcf_generate.cc {sorted(b) if b else 'every model'}, "
                       f"schema.py {sorted(s) if s else 'every model'}")
    return out


def sample_values(phase_space_file):
    """A valid, non-default command-line value for every option with a flag
    (other than model, beam and events)."""
    return {
        "cs": ["--cs", "cc2"],
        "q2": ["--q2", "2:8"],
        "sigma_cm": ["--sigma-cm", "0.12"],
        "wave_function": ["--wave-function", "AV18+UIX"],
        "theta_e": ["--theta-e", "10:30"],
        "kF": ["--kF", "0.3"],
        "cm": ["--no-cm"],
        "nn": ["--nn", "N2LO"],
        "estar": ["--estar", "0.02"],
        "kcut": ["--kcut", "0.3"],
        "onset": ["--onset", "local"],
        "prel_is": ["--prel-is"],
        "fermi_gas": ["--fermi-gas", "local"],
        "kF_fg": ["--kF-fg", "0.22"],
        "rad": ["--rad"],
        "coulomb": ["--coulomb"],
        "phase_space": ["--phase-space", phase_space_file],
    }


def run(binary, args):
    out = subprocess.run([binary] + args + ["--check"], capture_output=True, text=True,
                         timeout=300)
    return out.returncode, out.stdout + out.stderr


def main():
    binary = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "build", "gcf_generate")
    failures, checked = table_mismatches(), 0
    for f in failures:
        print(f"FAIL table: {f}")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as ps:
        ps.write("QSq 1 10\n")
    try:
        samples = sample_values(ps.name)
        flagged = [o for o in schema.GENERATE
                   if o.flag is not None and o.key not in ("model", "beam", "events")]
        missing = [o.key for o in flagged if o.key not in samples]
        if missing:
            print(f"FAIL: no test value for {missing}; add them to sample_values()")
            return 1
        cases = [(o.key, samples[o.key], o.models) for o in flagged]
        # Other spellings of the same options.
        cases += [("wave_function (-u)", ["-u", "AV18+UIX"], ["3N"]),
                  ("wave_function (--interaction)", ["--interaction", "2"], ["3N"]),
                  ("nn (default value)", ["--nn", "AV18"], ["2N"])]
        for model in schema.MODELS:
            rc, text = run(binary, ["--model", model] + BASE[model])
            if rc != 0:
                failures.append(f"{model} without options: exit {rc}: {text.strip()[-300:]}")
            for key, args, models in cases:
                checked += 1
                rc, text = run(binary, ["--model", model] + BASE[model] + args)
                if model in models:
                    ok = rc == 0
                    want = "accepted"
                else:
                    ok = rc == 1 and "applies to" in text
                    want = "rejected"
                status = "ok  " if ok else "FAIL"
                print(f"{status} {model:2s} {' '.join(args):38s} {want}")
                if not ok:
                    failures.append(f"{model} {' '.join(args)}: expected {want}, exit {rc}: "
                                    f"{text.strip()[-300:]}")
    finally:
        os.unlink(ps.name)
    print(f"{checked} model x option cases, {len(failures)} failures")
    for f in failures:
        print("  " + f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
