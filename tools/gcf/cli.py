"""The `gcf` command: run pipeline stages from a config file or the command line.

  gcf run CONFIG.toml            run the stages a config describes
  gcf scan CONFIG.toml --set ... run it for several values of an option
  gcf status                     the runs in a directory and how far they are
  gcf generate [options]         PWIA events -> event file   (gcf_generate)
  gcf fsi [options]              event file -> + GENIE FSI   (gcf_fsi)
  gcf quicklook FILE             a page of control plots, with an optional selection
  gcf summary FILE               the numbers of an event file
  gcf gui                        the web interface
  gcf doctor                     check that this machine is ready
  gcf example [3N|2N|MF|fsi|dis] print an example config
  gcf env                        print the environment the stages run with
"""

import argparse
import json
import os
import shlex
import signal
import subprocess
import sys
import time
import tomllib
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
    "dis": """\
# Tagged DIS off an SRC pair on 12C at 6 GeV, with one of the suite's other generators.
[run]
name = "c12_DIS"
output_dir = "runs"
stages = ["program"]
jobs = 4

[program]
program = "genDIS"
events = 200000
target = [12, 6]
beam = 6.0
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
        if o.key == "scale":
            parser.add_argument("--scale", dest="scale", action="append", default=None, metavar="CHANNEL=FACTOR",
                                type=_scale, help="scale a cascade channel by a factor (repeatable)")
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


def _scale(text):
    name, eq, factor = text.partition("=")
    try:
        return name, float(factor) if eq else None
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected CHANNEL=FACTOR, got {text!r}") from None


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


def run_plan(cfg, run_dir, final_outputs=None, quiet=False, resume=False):
    cfg = runner.resolve_seed(cfg)
    lock = runner.RunLock(run_dir)     # held until this process ends  # noqa: F841
    if "chunks" not in cfg["run"] or cfg["run"]["chunks"] is None:
        # Part of what makes the run reproducible, so the saved config records it.
        stages = cfg["run"]["stages"]
        total = cfg["generate"]["events"] if "generate" in stages else cfg["program"]["events"] if "program" in stages else None
        if total is not None:
            cfg["run"]["chunks"] = runner.checkpoint_jobs(cfg["run"], total)
    plan = runner.plan_from_config(cfg, run_dir, final_outputs, resume=resume)
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    (Path(run_dir) / "config.toml").write_text(
        config.dumps(cfg, comment="Resolved configuration of this run (re-run with: gcf run config.toml)"))
    if not quiet:
        print(f"Run '{cfg['run']['name']}': {', '.join(cfg['run']['stages'])}, "
              f"{len(plan.jobs)} job(s), base seed {cfg['run']['seed']}, logs in {run_dir}/logs")
    return runner.execute(plan, quiet=quiet)


def detach(argv):
    """Start this command again, on its own: no terminal, its output in a log."""
    log = Path(argv["log"])
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as out:
        subprocess.Popen([sys.executable, str(runner.REPO / "gcf"), *argv["args"]], stdout=out, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True, cwd=str(log.parent))


def cmd_run(args):
    cfg = config.load(args.config)
    if args.jobs:
        cfg["run"]["jobs"] = args.jobs
    run_dir = Path(cfg["run"]["output_dir"]) / cfg["run"]["name"]
    if args.merge_finished:
        partial = runner.merge_finished(run_dir)
        for stage, p in partial.items():
            print(f"{stage}: merged {p['jobs']} of {p['jobs_total']} jobs, {p['events']} of {p['total']} events")
        return True
    if args.resume:
        if not cfg["run"]["seed"]:
            raise config.ConfigError(f"--resume needs the seed of the run: use {run_dir}/config.toml")
    elif (run_dir / "status.json").exists() and not args.force:
        raise config.ConfigError(f"{run_dir} already holds a run; pick another [run] name, pass --resume "
                                 f"with its config.toml to finish it, or --force to overwrite it")
    if args.detach:
        again = ["run", "-q", str(Path(args.config).resolve())]
        again += ["--resume"] if args.resume else []
        again += ["--force"] if args.force else []
        again += ["-j", str(args.jobs)] if args.jobs else []
        detach({"args": again, "log": run_dir / "gcf.log"})
        print(f"Started {cfg['run']['name']} on its own; it goes on after you log out.\n"
              f"Follow it with: ./gcf status --runs-dir {cfg['run']['output_dir']}")
        return True
    return run_plan(cfg, run_dir, quiet=args.quiet, resume=args.resume)


def cmd_scan(args):
    """Run the points of a scan one after the other, in this process."""
    from . import gui, scan
    with open(args.config, "rb") as f:
        raw = tomllib.load(f)
    axes = []
    for item in args.set:
        path, eq, values = item.partition("=")
        if not eq:
            raise config.ConfigError(f"--set expects OPTION=VALUES, got {item!r}")
        axes.append({"path": path, "values": values})
    base = Path(args.config).resolve().parent
    run = raw.setdefault("run", {})
    if args.name:
        run["name"] = args.name
    try:
        points = scan.expand(raw, axes)
        parsed = [(name, values, config.validate(point, base_dir=base)) for name, values, point in points]
    except ValueError as e:
        raise config.ConfigError(str(e)) from None
    runs_dir = Path(parsed[0][2]["run"]["output_dir"])
    runs_dir.mkdir(parents=True, exist_ok=True)
    group = raw["run"]["name"]
    record = runs_dir / gui.SCANS_DIR / f"{group}.json"
    record.parent.mkdir(exist_ok=True)
    record.write_text(json.dumps({"name": group, "created": time.time(),
                                  "axes": [{"path": a["path"], "label": scan.axis_label(a["path"])} for a in axes],
                                  "points": [{"name": n, "values": v} for n, v, _ in parsed]}, indent=1))
    if args.detach:
        again = ["scan", str(Path(args.config).resolve()), "-q", "--name", group]
        for item in args.set:
            again += ["--set", item]
        detach({"args": again, "log": runs_dir / gui.SCANS_DIR / f"{group}.log"})
        print(f"Started scan {group} with {len(parsed)} runs, on its own.\nFollow it with: ./gcf status --runs-dir {runs_dir}")
        return True
    ok = True
    for name, _, cfg in parsed:
        if (runs_dir / name / "status.json").exists():
            print(f"{name}: already there, skipping")
            continue
        ok = run_plan(cfg, runs_dir / name, quiet=args.quiet) and ok
    return ok


def cmd_status(args):
    from . import gui
    runs_dir = Path(args.runs_dir) if args.runs_dir else runner.default_runs_dir()
    app = gui.App(runs_dir, serve=False)
    runs = app.list_runs()
    if args.json:
        print(json.dumps(runs, indent=1))
        return True
    if not runs:
        print(f"No runs in {runs_dir}")
        return True
    print(f"Runs in {runs_dir}:")
    for r in runs:
        progress = ", ".join(f"{k} {100.0 * v['done'] / max(v['total'], 1):.0f}%" for k, v in r["stages"].items())
        print(f"  {r['name']:32s} {r['state']:10s} {progress:28s} {r['setup']}")
        if r["state"] in ("stopped", "failed", "cancelled") and r["resumable"]:
            print(f"      resume with: ./gcf run --resume --detach {runs_dir / r['name'] / 'config.toml'}")
    return True


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
    if "scale" in fsi:
        fsi["scale"] = dict(fsi["scale"])
    if args.target:
        fsi["target"] = args.target
    raw = {"run": {"name": Path(args.output).stem, "stages": ["fsi"], "jobs": args.jobs,
                   "seed": args.seed, "keep_chunks": args.keep_chunks},
           "fsi": fsi}
    cfg = config.validate(raw)
    return run_plan(cfg, work_dir_for(args.output), {"fsi": args.output})


def cmd_quicklook(args):
    from . import quicklook
    cuts = {k: getattr(args, k) for k in quicklook.CUT_LABELS if getattr(args, k) is not None}
    stats = {}
    out = quicklook.make(args.file, args.output, cuts=cuts, stats=stats)
    print(f"Wrote {out}")
    if cuts:
        print(f"{stats['description']}: {stats['events']} events, N_eff = {stats['n_eff']:.0f}, "
              f"{100 * stats['fraction']:.3g}% of sum(w)")
    return True


def cmd_summary(args):
    from . import results
    s = results.summarize(args.file)
    print(json.dumps(s, indent=1) if args.json else results.text(s))
    return True


def cmd_gui(args):
    from . import gui
    running = gui.read_server(args.port)
    token = gui.stable_token(args.new_token)
    if args.status or args.stop:
        if not running:
            print(f"No gcf gui of yours is running on port {args.port}.")
            return not args.stop
        if args.stop:
            os.kill(running["pid"], signal.SIGTERM)
            print(f"Stopped the gcf gui on port {args.port}. Runs already started keep going.")
            return True
        print(f"gcf gui is running on port {args.port}, with runs in {running['runs_dir']}:\n"
              f"    http://localhost:{args.port}/?token={token}")
        return True
    runs_dir = Path(args.runs_dir).expanduser().resolve() if args.runs_dir else runner.default_runs_dir()
    if args.background:
        if running:
            print(f"Already running on port {args.port}, with runs in {running['runs_dir']}:\n"
                  f"    http://localhost:{args.port}/?token={token}")
            return True
        runs_dir.mkdir(parents=True, exist_ok=True)
        again = ["gui", "--port", str(args.port), "--host", args.host, "--runs-dir", str(runs_dir), "--quiet"]
        again += ["--max-cpus", str(args.max_cpus)] if args.max_cpus else []
        log = gui.config_dir() / f"gui-{args.port}.log"
        detach({"args": again, "log": log})
        deadline = time.time() + 20
        while time.time() < deadline and not gui.read_server(args.port):
            time.sleep(0.1)
        if not gui.read_server(args.port):
            raise runner.RunError(f"the interface did not start; see {log}")
        print(gui.banner(gui.App(runs_dir, serve=False, token=token), args.host, args.port))
        print(f"\nRunning on its own: log out whenever you like. Stop it with: ./gcf gui --stop --port {args.port}\n"
              f"Find this URL again with: ./gcf gui --status --port {args.port}")
        return True
    if running:
        raise runner.RunError(f"a gcf gui of yours is already running on port {args.port}; see gcf gui --status, "
                              f"or use another --port")
    gui.serve(args.host, args.port, runs_dir, token=token, max_cpus=args.max_cpus, quiet=args.quiet)
    return True


def cmd_doctor(args):
    from . import doctor
    results = doctor.checks(args.runs_dir)
    print(doctor.text(results))
    return all(c["ok"] or not c["required"] for c in results)


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
    p.add_argument("--resume", action="store_true",
                   help="with a run's config.toml: finish it, keeping the jobs that finished and taking over "
                        "those still running")
    p.add_argument("--detach", action="store_true",
                   help="run on its own, without this terminal; follow it with gcf status")
    p.add_argument("--merge-finished", action="store_true",
                   help="with a stopped run's config.toml: make the final files from the jobs that finished")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("scan", help="run a config for several values of one or more options, one after the other")
    p.add_argument("config")
    p.add_argument("--set", action="append", required=True, metavar="OPTION=VALUES",
                   help="e.g. generate.sigma_cm=0.1,0.15,0.2 or fsi.scale:pion.absorption=0.5:1.5:5 (start:stop:count); repeat for a grid")
    p.add_argument("--name", help="the scan's name (default: the config's [run] name)")
    p.add_argument("--detach", action="store_true", help="run on its own, without this terminal")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_scan)

    p = sub.add_parser("status", help="the runs in a directory and how far they are")
    p.add_argument("--runs-dir", help="default: as for gcf gui")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_status)

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
    for flag, key, what in (("--q2", "q2", "Q^2 [GeV^2]"), ("--xb", "xb", "x_B"), ("--theta-e", "theta_e", "electron angle [deg]"),
                            ("--pmiss", "pmiss", "p_miss before FSI [GeV/c]"),
                            ("--protons", "protons", "number of protons (after FSI if the file has it)")):
        p.add_argument(flag, dest=key, metavar="MIN:MAX", type=_range, help=f"select on {what}")
    p.set_defaults(func=cmd_quicklook)

    p = sub.add_parser("summary", help="the numbers of an event file: N_eff, sum(w)/attempts, provenance")
    p.add_argument("file")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_summary)

    p = sub.add_parser("gui", help="web interface (open it through an SSH tunnel)")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--host", default="127.0.0.1",
                   help="interface to listen on (default: localhost only)")
    p.add_argument("--runs-dir", help="where GUI runs are written (default: $GCF_RUNS_DIR, else <repo>/runs "
                                      "if you can write there, else ~/gcf_runs)")
    p.add_argument("--background", action="store_true",
                   help="run on its own, so that it stays up after you log out; prints its URL")
    p.add_argument("--status", action="store_true", help="say whether it is running, and its URL")
    p.add_argument("--stop", action="store_true", help="stop one started with --background; runs keep going")
    p.add_argument("--max-cpus", type=int, help="CPUs its runs may use together; more runs wait in a queue "
                                                "(default: all of this machine's)")
    p.add_argument("--new-token", action="store_true", help="replace the token in the URL, so old links stop working")
    p.add_argument("--quiet", action="store_true", help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_gui)

    p = sub.add_parser("doctor", help="check that this machine is ready, and say how to fix what isn't")
    p.add_argument("--runs-dir")
    p.set_defaults(func=cmd_doctor)

    p = sub.add_parser("example", help="print an example config")
    p.add_argument("kind", nargs="?", default="3N", choices=sorted(EXAMPLES))
    p.set_defaults(func=cmd_example)

    p = sub.add_parser("env", help="print the environment the stages run with")
    p.set_defaults(func=cmd_env)

    args = parser.parse_args(argv)
    try:
        ok = args.func(args)
    except (config.ConfigError, runner.RunError, RuntimeError) as e:
        print(f"gcf: error: {e}", file=sys.stderr)
        return 2
    return 0 if ok else 1
