"""The `gcf` command: run pipeline stages from a config file or the command line.

  gcf run CONFIG.toml            run the stages a config describes
  gcf generate [options]         PWIA events -> event file   (gcf_generate)
  gcf fsi [options]              event file -> + GENIE FSI   (gcf_fsi)
  gcf quicklook FILE             a page of control plots
  gcf gui                        the web interface
  gcf example [3N|2N|MF|fsi]     print an example config
  gcf env                        print the environment the stages run with
"""

import argparse
import shlex
import sys
from pathlib import Path

from . import config, runner, schema

EXAMPLES = {
    "3N": """\
# 3N SRC events on 12C at 6 GeV, then GENIE hN FSI.
[run]
name = "c12_3N_6GeV"
output_dir = "runs"
stages = ["generate", "fsi"]
jobs = 8
seed = 0              # 0: draw one; the saved config records it

[generate]
model = "3N"
beam = 6.0
events = 100000
target = [12, 6]
wave_function = "AV18"  # 3He Faddeev density; AV18+UIX, CDBonn, ... or AV18-pD, ...
theta_e = [5.0, 40.0]

[fsi]
model = "hN"            # hN (full cascade) or hA (effective)
transport = "shared"
off = []                # e.g. ["pion.absorption", "nucleon.inelastic"]; see gcf_fsi --list-channels
""",
    "2N": """\
# 2N SRC pair events on 12C at 4.5 GeV, then GENIE hN FSI.
[run]
name = "c12_2N_4p5GeV"
output_dir = "runs"
stages = ["generate", "fsi"]
jobs = 8

[generate]
model = "2N"
beam = 4.5
events = 100000
target = [12, 6]
onset = "fixed"

[fsi]
model = "hN"
""",
    "MF": """\
# Mean-field (single nucleon, local Fermi gas) events on 12C, PWIA only.
[run]
name = "c12_MF"
output_dir = "runs"
stages = ["generate"]
jobs = 4

[generate]
model = "MF"
beam = 4.5
events = 50000
target = [12, 6]
fermi_gas = "local"
""",
    "fsi": """\
# GENIE FSI on an existing event file (from gcf_generate or any generator
# writing the GCF event-file format, see docs/PIPELINE.md).
[run]
name = "my_events_hA"
output_dir = "runs"
stages = ["fsi"]
jobs = 8

[fsi]
input = "my_events.root"
model = "hA"
target = [12, 6]      # only needed if the file has no target_A/target_Z
""",
}


def add_schema_args(parser, opts, skip=()):
    for o in opts:
        if o.key in skip:
            continue
        help_text = o.help + (f" [{o.unit}]" if o.unit else "")
        if o.key == "target":
            continue                      # handled by the caller
        if o.key == "off":
            parser.add_argument("--off", dest="off", action="append", default=None, metavar="CHANNEL",
                                help="switch a cascade channel off (repeatable; gcf_fsi --list-channels)")
            continue
        if o.key == "transport":
            parser.add_argument("--independent", dest="transport", action="store_const",
                                const="independent", default=None,
                                help="one GENIE record per nucleon (default: shared remnant)")
            continue
        if o.flag is None:
            continue
        if o.kind == "bool":
            parser.add_argument(o.flag, dest=o.key, default=None,
                                action="store_false" if o.invert else "store_true",
                                help=help_text or o.label)
        elif o.kind == "range":
            parser.add_argument(o.flag, dest=o.key, default=None, metavar="MIN:MAX",
                                type=_range, help=help_text or o.label)
        else:
            typ = {"float": float, "int": int, "str": str}[o.kind]
            parser.add_argument(o.flag, dest=o.key, default=None, type=typ,
                                choices=o.choices, help=help_text or o.label)


def _range(text):
    try:
        lo, hi = (float(x) for x in text.split(":"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected MIN:MAX, got {text!r}") from None
    return [lo, hi]


def _target(text):
    try:
        a, z = (int(x) for x in text.split(":"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected A:Z, got {text!r}") from None
    return [a, z]


def section_from_args(ns, opts):
    return {o.key: getattr(ns, o.key) for o in opts
            if getattr(ns, o.key, None) is not None}


def add_run_args(parser):
    parser.add_argument("-j", "--jobs", type=int, default=1, help="parallel jobs (default 1)")
    parser.add_argument("-S", "--seed", type=int, default=0,
                        help="base seed; job i uses seed + i (default: drawn and printed)")
    parser.add_argument("--keep-chunks", action="store_true", help="keep the per-job files")


def work_dir_for(output):
    output = Path(output).resolve()
    return output.parent / f".{output.name}.gcf"


def run_plan(cfg, run_dir, final_outputs=None, quiet=False):
    cfg = runner.resolve_seed(cfg)
    plan = runner.plan_from_config(cfg, run_dir, final_outputs)
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    (Path(run_dir) / "config.toml").write_text(
        config.dumps(cfg, comment="Resolved configuration of this run (re-run with: gcf run config.toml)"))
    if not quiet:
        print(f"Run '{cfg['run']['name']}': {', '.join(cfg['run']['stages'])}, "
              f"{len(plan.jobs)} job(s), base seed {cfg['run']['seed']}, logs in {run_dir}/logs")
    return runner.execute(plan, quiet=quiet)


def cmd_run(args):
    cfg = config.load(args.config)
    if args.jobs:
        cfg["run"]["jobs"] = args.jobs
    run_dir = Path(cfg["run"]["output_dir"]) / cfg["run"]["name"]
    if (run_dir / "status.json").exists() and not args.force:
        raise config.ConfigError(f"{run_dir} already holds a run; pick another [run] name "
                                 f"or pass --force to overwrite it")
    return run_plan(cfg, run_dir, quiet=args.quiet)


def cmd_generate(args):
    gen = section_from_args(args, schema.GENERATE)
    gen["target"] = [args.target_A, args.target_Z]
    stages = ["generate"] + (["fsi"] if args.fsi_output else [])
    raw = {"run": {"name": Path(args.output).stem, "stages": stages, "jobs": args.jobs,
                   "seed": args.seed, "keep_chunks": args.keep_chunks},
           "generate": gen}
    if args.fsi_output:
        raw["fsi"] = {"model": args.fsi_model}
    cfg = config.validate(raw)
    finals = {"generate": args.output}
    if args.fsi_output:
        finals["fsi"] = args.fsi_output
    return run_plan(cfg, work_dir_for(args.output), finals)


def cmd_fsi(args):
    fsi = section_from_args(args, schema.FSI)
    if args.target:
        fsi["target"] = args.target
    raw = {"run": {"name": Path(args.output).stem, "stages": ["fsi"], "jobs": args.jobs,
                   "seed": args.seed, "keep_chunks": args.keep_chunks},
           "fsi": fsi}
    cfg = config.validate(raw)
    return run_plan(cfg, work_dir_for(args.output), {"fsi": args.output})


def cmd_quicklook(args):
    from . import quicklook
    out = quicklook.make(args.file, args.output)
    print(f"Wrote {out}")
    return True


def cmd_gui(args):
    from . import gui
    gui.serve(args.host, args.port, Path(args.runs_dir))
    return True


def cmd_example(args):
    print(EXAMPLES[args.kind], end="")
    return True


def cmd_env(args):
    env = runner.environment()
    for key in ("PATH", "LD_LIBRARY_PATH", "ROOTSYS", "GENIE", "GXMLPATH"):
        if key in env:
            print(f"export {key}={shlex.quote(env[key])}")
    print(f"# stage binaries: {runner.BIN_DIR}")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gcf", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split("\n\n", 1)[0],
        epilog="Each command has --help. Docs: docs/PIPELINE.md")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="run the stages described by a TOML config")
    p.add_argument("config")
    p.add_argument("-j", "--jobs", type=int, default=0, help="override [run] jobs")
    p.add_argument("--force", action="store_true", help="overwrite an existing run")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("generate", help="PWIA events -> event file (optionally + FSI)")
    p.add_argument("-o", "--output", required=True, help="output event file")
    p.add_argument("-A", "--target-A", dest="target_A", type=int, default=12)
    p.add_argument("-Z", "--target-Z", dest="target_Z", type=int, default=6)
    p.add_argument("--fsi-output", help="also run FSI (hN unless --fsi-model) into this file")
    p.add_argument("--fsi-model", default="hN", choices=["hN", "hA"])
    add_schema_args(p, schema.GENERATE)
    add_run_args(p)
    p.set_defaults(func=cmd_generate)

    p = sub.add_parser("fsi", help="add GENIE FSI to an event file")
    p.add_argument("-o", "--output", required=True, help="output file")
    p.add_argument("--target", type=_target, metavar="A:Z",
                   help="target for files without target_A/target_Z")
    add_schema_args(p, schema.FSI)
    add_run_args(p)
    p.set_defaults(func=cmd_fsi)

    p = sub.add_parser("quicklook", help="control plots of an event file (PNG)")
    p.add_argument("file")
    p.add_argument("-o", "--output", help="PNG path (default: next to the file)")
    p.set_defaults(func=cmd_quicklook)

    p = sub.add_parser("gui", help="web interface (open it through an SSH tunnel)")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--host", default="127.0.0.1",
                   help="interface to listen on (default: localhost only)")
    p.add_argument("--runs-dir", default=str(runner.REPO / "runs"),
                   help="where GUI runs are written (default: <repo>/runs)")
    p.set_defaults(func=cmd_gui)

    p = sub.add_parser("example", help="print an example config")
    p.add_argument("kind", nargs="?", default="3N", choices=sorted(EXAMPLES))
    p.set_defaults(func=cmd_example)

    p = sub.add_parser("env", help="print the environment the stages run with")
    p.set_defaults(func=cmd_env)

    args = parser.parse_args(argv)
    try:
        ok = args.func(args)
    except (config.ConfigError, runner.RunError) as e:
        print(f"gcf: error: {e}", file=sys.stderr)
        return 2
    return 0 if ok else 1
