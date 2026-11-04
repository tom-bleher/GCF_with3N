"""Running the pipeline stages: environment, parallel jobs, merging, status.

A run is split into jobs. Each job runs its steps (gcf_generate and/or
gcf_fsi) one after the other as subprocesses; jobs run in parallel. Their
outputs are then merged with hadd. Progress is parsed from the stages'
"Progress: n/N" lines and written to <run_dir>/status.json, which the GUI
reads.

The eic stage is one job: tools/run_eic_sim.sh for each beam setting and Q^2
slice, one step per stage of that script, then the analysis. That script hides
failures inside the container from its exit code, so each step checks what it
was supposed to leave behind, and progress is counted from the chunks' files.
"""

import json
import os
import re
import secrets
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

from . import schema

REPO = Path(__file__).resolve().parents[2]
BIN_DIR = Path(os.environ.get("GCF_BIN_DIR", REPO / "src/programs/gcf_stages/build"))
GENIE_DIR = REPO / "src/programs/genQE_FSI/Generator-R-3_06_02"
FSI_CONFIG = REPO / "src/programs/genQE_FSI/config"


class RunError(Exception):
    pass


def proc_start(pid):
    """Start time (clock ticks since boot) of a process, from /proc; None if
    it does not exist. Together with the pid it identifies a process: a
    reused pid gets a different start time."""
    try:
        with open(f"/proc/{pid}/stat") as f:
            fields = f.read().rsplit(")", 1)[1].split()
        return int(fields[19])        # field 22 of /proc/<pid>/stat
    except (OSError, ValueError, IndexError):
        return None


def proc_is(pid, start):
    """True if `pid` is still the process that had start time `start`."""
    return pid is not None and start is not None and proc_start(pid) == start


def binary(name):
    path = BIN_DIR / name
    if not os.access(path, os.X_OK):
        raise RunError(f"{name} is not built ({path} missing). Build it with ./install.sh, "
                       f"or set GCF_BIN_DIR to where it is.")
    return str(path)


_env_cache = None


def environment():
    """Environment for the stage binaries: ROOT, GENIE and log4cpp.

    Taken from build/gcf_env.sh (written by install.sh) or, failing that, the
    tracked setup_env.sh, by sourcing it in bash.
    """
    global _env_cache
    if _env_cache is not None:
        return dict(_env_cache)
    env = dict(os.environ)
    for script in (REPO / "build/gcf_env.sh", REPO / "setup_env.sh"):
        if script.is_file():
            out = subprocess.run(["bash", "-c", f'source "{script}" >/dev/null 2>&1; env -0'],
                                 capture_output=True, check=False)
            if out.returncode == 0:
                env = dict(item.split("=", 1) for item in
                           out.stdout.decode(errors="replace").split("\0") if "=" in item)
                break
    env.setdefault("GENIE", str(GENIE_DIR))
    env.setdefault("GXMLPATH", f"{FSI_CONFIG}/genie_override:{env['GENIE']}/config/G18_10a:"
                               f"{env['GENIE']}/config")
    _env_cache = env
    return dict(env)


def hadd_command():
    hadd = shutil.which("hadd", path=environment().get("PATH"))
    if not hadd:
        raise RunError("hadd not found: ROOT is not set up (build/gcf_env.sh or setup_env.sh)")
    return hadd


# ---------------------------------------------------------------------------
# Command lines

def generate_args(gen, events, output, seed):
    args = [binary("gcf_generate"), "--model", gen["model"], "--beam", repr(gen["beam"]),
            "--events", str(events), "--output", str(output), "--seed", str(seed),
            "-A", str(gen["target"][0]), "-Z", str(gen["target"][1])]
    for o in schema.GENERATE:
        if o.key in ("model", "beam", "events", "target") or o.flag is None:
            continue
        v = gen.get(o.key)
        if v is None or gen["model"] not in o.models:
            continue
        if o.kind == "bool":
            if v != o.invert:   # pass the flag when it changes the binary's default
                args.append(o.flag)
        elif o.kind == "range":
            args += [o.flag, f"{v[0]}:{v[1]}"]
        else:
            args += [o.flag, str(v)]
    return args


def fsi_args(fsi, inp, output, seed, first=None, count=None, copy_objects=True):
    args = [binary("gcf_fsi"), "--input", str(inp), "--output", str(output),
            "--seed", str(seed), "--fsi-model", fsi["model"],
            "--vertex-power", fsi["vertex_power"]]
    if fsi["transport"] == "independent":
        args.append("--independent")
    if fsi.get("target"):
        args += ["--target", f"{fsi['target'][0]}:{fsi['target'][1]}"]
    if fsi.get("keep"):
        args += ["--keep", fsi["keep"]]
    for c in fsi.get("off") or []:
        args += ["--channel", f"{c}=0"]
    if fsi.get("nucleon_mfp_scale") is not None:
        args += ["--channel", f"nucleon.mfp={fsi['nucleon_mfp_scale']}"]
    if fsi.get("pion_mfp_scale") is not None:
        args += ["--channel", f"pion.mfp={fsi['pion_mfp_scale']}"]
    if first is not None:
        args += ["--first-event", str(first), "--max-events", str(count)]
    if not copy_objects:
        args.append("--no-copy-objects")
    return args


def check_setup(cfg):
    """Run the stages' own setup checks on a validated config, without
    generating anything: gcf_generate --check (nucleus tables, option
    combinations) or, for FSI-only runs, that the input file and tree can be
    read. Returns a list of problems (empty = fine)."""
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
    if "eic" in run["stages"]:
        return check_eic(cfg)
    problems = []
    env = environment()
    if "generate" in run["stages"]:
        args = generate_args(gen, 1, "/dev/null", 1)
        args = [a for a in args if a not in ("--output", "/dev/null")] + ["--check"]
        out = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120)
        if out.returncode != 0:
            msg = [l for l in (out.stderr + out.stdout).splitlines()
                   if l.strip() and not l.startswith(("eNCrossSection", "You are using"))]
            text = " ".join(msg[-3:]) or f"exit {out.returncode}"
            text = re.sub(r"\s*Add it to constants/nuclearMasses\.hh.*$", "", text)
            text = re.sub(r"^\w+::\w+: ", "", text)
            problems.append("[generate] target: the generator can't run this setup: " + text)
    else:
        try:
            n = count_entries(fsi["input"])
            if n == 0:
                problems.append(f"{fsi['input']} has no events")
        except RunError as e:
            problems.append(str(e))
    if "fsi" in run["stages"] and not (GENIE_DIR / "lib").is_dir():
        problems.append("GENIE is not built (no lib/ in the vendored tree); run ./install.sh")
    return problems


def count_entries(path):
    try:
        out = subprocess.run([binary("gcf_fsi"), "--input", str(path), "--count"],
                             capture_output=True, text=True, env=environment(), timeout=60)
    except subprocess.TimeoutExpired:
        raise RunError(f"reading {path} took more than 60 s") from None
    for line in out.stdout.splitlines():
        if line.startswith("Entries:"):
            return int(line.split()[1])
    raise RunError(f"cannot read {path}: {out.stderr.strip() or out.stdout.strip()}")


def split(total, jobs):
    """Split `total` into `jobs` near-equal positive parts."""
    jobs = max(1, min(jobs, total))
    base, rem = divmod(total, jobs)
    return [base + (1 if i < rem else 0) for i in range(jobs)]


# ---------------------------------------------------------------------------
# Plans: jobs of sequential steps, then merges

class Step:
    def __init__(self, stage, args, output, total, log, env=None, verify=None, poll=None, context=None):
        self.stage, self.args, self.output, self.total, self.log = stage, args, output, total, log
        self.env = env               # None: environment()
        self.verify = verify         # after exit 0: returns an error message, or None if all is well
        self.poll = poll             # while running: updates self.done / self.total from files
        self.context = context       # a log worth quoting when the step fails, besides its own
        self.done = 0
        self.state = "pending"       # pending, running, done, failed, cancelled
        self.returncode = None
        self.pid = None              # of the stage process while it runs ...
        self.pid_start = None        # ... and its start time, to identify it


class Plan:
    def __init__(self, run_dir, name):
        self.run_dir = Path(run_dir)
        self.name = name
        self.jobs = []               # list of lists of Step
        self.merges = []             # (stage, [inputs], output)
        self.outputs = {}            # stage -> final file
        self.cleanup = []            # files to delete after a successful run
        self.eic = False             # the ePIC chain: its steps poll for progress


def plan_from_config(cfg, run_dir, final_outputs=None, resume=False):
    """Jobs and merges for a validated config.

    Final files are <run_dir>/<name>_pwia.root and <name>_fsi.root unless
    final_outputs maps a stage to another path. With one job the stages
    write their final files directly; with several, per-job files go to
    <run_dir>/chunks and are merged with hadd.
    """
    if "eic" in cfg["run"]["stages"]:
        return eic_plan(cfg, run_dir, resume)
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
    name = run["name"]
    run_dir = Path(run_dir)
    chunks, logs = run_dir / "chunks", run_dir / "logs"
    plan = Plan(run_dir, name)
    seed = run["seed"]
    stages = run["stages"]
    final = {"generate": run_dir / f"{name}_pwia.root", "fsi": run_dir / f"{name}_fsi.root"}
    final.update({k: Path(v) for k, v in (final_outputs or {}).items()})
    for stage in stages:
        plan.outputs[stage] = final[stage]

    if "generate" in stages:
        parts = split(gen["events"], run["jobs"])
        single = len(parts) == 1
        outs = {"generate": [], "fsi": []}
        for i, n in enumerate(parts):
            g = final["generate"] if single else chunks / f"pwia_{i:03d}.root"
            steps = [Step("generate", generate_args(gen, n, g, seed + i), g, n,
                          logs / f"generate_{i:03d}.log")]
            outs["generate"].append(g)
            if "fsi" in stages:
                f = final["fsi"] if single else chunks / f"fsi_{i:03d}.root"
                steps.append(Step("fsi", fsi_args(fsi, g, f, seed + 1000003 + i), f, n,
                                  logs / f"fsi_{i:03d}.log"))
                outs["fsi"].append(f)
            plan.jobs.append(steps)
    else:
        total = count_entries(fsi["input"])
        if total == 0:
            raise RunError(f"{fsi['input']} has no events")
        parts = split(total, run["jobs"])
        single = len(parts) == 1
        outs = {"fsi": []}
        first = 0
        for i, n in enumerate(parts):
            f = final["fsi"] if single else chunks / f"fsi_{i:03d}.root"
            # Only the first chunk carries hAttempts etc., so the merge
            # doesn't add them up once per chunk.
            plan.jobs.append([Step("fsi", fsi_args(fsi, fsi["input"], f, seed + 1000003 + i,
                                                   first, n, copy_objects=(i == 0)),
                                   f, n, logs / f"fsi_{i:03d}.log")])
            outs["fsi"].append(f)
            first += n

    if not single:
        for stage, files in outs.items():
            if files:
                plan.merges.append((stage, files, final[stage]))
                if not run["keep_chunks"]:
                    plan.cleanup += files
    return plan


def resolve_seed(cfg):
    if cfg["run"]["seed"] == 0:
        top = schema.EIC_MAX_SEED if "eic" in cfg["run"]["stages"] else 2**31 - 2**21
        cfg["run"]["seed"] = 1 + secrets.randbelow(top)
    return cfg


# ---------------------------------------------------------------------------
# The ePIC chain

RUN_EIC_SIM = REPO / "tools/run_eic_sim.sh"
EIC_ANALYSIS = REPO / "tools/eic_src_analysis.py"
EIC_GENERATOR_BINARIES = {"genQE": "src/programs/genQE_FSI/build/genQE", "genQE_FSI": "src/programs/genQE_FSI/build/genQE_FSI",
                          "genQE_3N_FSI": "src/programs/genQE_3N/build/genQE_3N_FSI", "genDIS": "build/programs/genDIS/genDIS"}


def eic_slices(eic):
    """(label, beams, q2 min, q2 max, attempts) of every run_eic_sim.sh run of a config, in running order."""
    return [(f"{b}_q2_{lo:g}_{hi:g}", b, lo, hi, n) for b in eic["beams"] for lo, hi, n in eic["slices"]]


def eic_env(eic):
    env = dict(os.environ, EPIC_INSTALL=eic["epic_install"], EICRECON_INSTALL=eic["eicrecon_install"],
               EIC_SHELL=eic["eic_shell"], FIELDMAPS=eic["fieldmaps"])
    # these would override the per-beam defaults of run_eic_sim.sh without the saved config saying so
    for key in ("EPIC_CONFIG", "FIELD", "AB_PRESET"):
        env.pop(key, None)
    return env


def eic_sim_args(eic, beams, lo, hi, seed, out, stage, attempts=None):
    sample = schema.EIC_SAMPLES[eic["sample"]]
    args = [str(RUN_EIC_SIM), "--beams", beams, "--ion", eic["ion"], "--out", str(out), "--stage", stage,
            "--q2", f"{lo:g}:{hi:g}", "--seed", str(seed), "--generator", sample.generator,
            "--chunk-size", str(eic["chunk_size"]), "--jobs", str(eic["sim_jobs"]),
            "--sim-threads", str(eic["sim_threads"]), "--reco-jobs", str(eic["reco_jobs"]), "--reco-threads", "1",
            "--compress", str(eic["compress"]), "--genie-compound", eic["genie_compound"]]
    if attempts is not None:
        args += ["--attempts", str(attempts)]
    gen_args = " ".join(a for a in (sample.gen_args, eic.get("gen_args")) if a)
    if gen_args:
        args += ["--gen-args", gen_args]
    if not sample.cuts:
        args.append("--no-cuts")
    if eic["deexcite"] != "none":
        args += ["--deexcite", eic["deexcite"]]
    if not eic["beam_nucleon"]:
        args.append("--no-beam-nucleon")
    return args


def _read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _finished_chunks(out, step_name):
    """Chunks of a stage that have ended: their exit_* file is empty while the chunk runs."""
    n = 0
    try:
        with os.scandir(out) as entries:
            for e in entries:
                if e.name.startswith(f"exit_{step_name}_") and e.stat().st_size > 0:
                    n += 1
    except OSError:
        pass
    return n


def _chunk_progress(out, step_names, chunk_size):
    """poll callback of a step: events done, from finished chunks, and the true total once known."""
    def poll(step):
        written = (_read_json(Path(out) / "gen.hepmc3.json") or {}).get("events_written")
        if written:
            step.total = written * len(step_names)
        step.done = min(step.total, sum(_finished_chunks(out, s) for s in step_names) * chunk_size)
    return poll


def _tail(path, lines=12):
    try:
        return "\n".join(Path(path).read_text(errors="replace").splitlines()[-lines:])
    except OSError:
        return ""


def _verify_generated(out):
    def verify():
        info = _read_json(Path(out) / "gen.hepmc3.json")
        if info is None:
            return f"the generator stage left no gen.hepmc3.json in {out}:\n" + _tail(Path(out) / "run.log")
        if not info.get("events_written"):
            return f"no event of {info.get('attempts', 0):g} attempts passed the cuts in {out}: raise the attempts"
        return None
    return verify


def _verify_chunks(out, step_names):
    """run_eic_sim.sh exits 0 whatever happened in the container; run.json has the outcome."""
    def verify():
        jobs = (_read_json(Path(out) / "run.json") or {}).get("jobs")
        if jobs is None:
            return f"no run.json in {out}:\n" + _tail(Path(out) / "run.log")
        problems = []
        for s in step_names:
            if jobs.get(f"{s}_failed"):
                failed = sorted(p.name.rsplit("_", 1)[1] for p in Path(out).glob(f"exit_{s}_*")
                                if p.read_text().strip() != "0")
                log = "npsim" if s == "sim" else "eicrecon"
                problems.append(f"{s}: {jobs[f'{s}_failed']} of {jobs[f'{s}_ok'] + jobs[f'{s}_failed']} chunks failed "
                                f"(see {log}_{failed[0]}.log)" if failed else f"{s}: chunks failed")
            elif not jobs.get(f"{s}_ok"):
                problems.append(f"{s}: no chunk ran; the chain stopped inside the container")
        if not problems:
            return None
        return (f"{out}: " + "; ".join(problems) + ". Resume reruns only what did not succeed.\n"
                + _tail(Path(out) / "run.log") + "\n" + _tail(Path(out) / "container_stage.log", 6))
    return verify


def eic_slice_finished(out):
    """True if a run_eic_sim.sh run has all its chunks simulated and reconstructed."""
    jobs = (_read_json(Path(out) / "run.json") or {}).get("jobs")
    return bool(jobs and jobs["sim_ok"] and jobs["reco_ok"] and not jobs["sim_failed"] and not jobs["reco_failed"])


def container_command(eic_shell, command):
    """argv running a shell command inside eic-shell, with its exit code. eic-shell's `--` form splits on spaces."""
    return ["bash", "-c", 'set -o pipefail; printf "%s\\n" "$2" | "$1"', "bash", eic_shell, command]


def eic_plan(cfg, run_dir, resume=False):
    """One job: for each beam setting and Q^2 slice generate, simulate, reconstruct; then the figures.

    With resume, finished slices are skipped and a slice that has its chunks restarts with run_eic_sim.sh
    --stage resume, which runs only the chunks that have not succeeded; one stopped earlier starts over, and
    regenerates the same events since the seed is in the saved config.
    """
    run, eic = cfg["run"], cfg["eic"]
    run_dir = Path(run_dir)
    logs = run_dir / "logs"
    plan = Plan(run_dir, run["name"])
    plan.eic = True
    env = eic_env(eic)
    sample = schema.EIC_SAMPLES[eic["sample"]]
    chunk = eic["chunk_size"]
    steps, by_setting = [], {}
    for i, (label, beams, lo, hi, attempts) in enumerate(eic_slices(eic)):
        out = run_dir / label
        seed = run["seed"] + i
        by_setting.setdefault(beams, []).append(out)
        args = lambda stage, n=None: eic_sim_args(eic, beams, lo, hi, seed, out, stage, n)
        events = max(1, round(attempts * sample.pass_fraction))
        if resume and eic_slice_finished(out):
            continue
        if resume and (out / "geometry.xml").is_file() and (out / "prefetch/calibrations").is_dir() \
                and any(out.glob("chunk_*.hepmc")):
            steps.append(Step("resume", args("resume"), out / "run.json", 2 * events, logs / f"resume_{label}.log", env,
                              _verify_chunks(out, ("sim", "reco")), _chunk_progress(out, ("sim", "reco"), chunk)))
            continue
        steps.append(Step("generate", args("gen", attempts), out / "run.json", attempts, logs / f"generate_{label}.log", env,
                          _verify_generated(out), context=out / "generator.log"))
        steps.append(Step("simulate", args("sim"), out / "run.json", events, logs / f"simulate_{label}.log", env,
                          _verify_chunks(out, ("sim",)), _chunk_progress(out, ("sim",), chunk)))
        steps.append(Step("reconstruct", args("reco"), out / "run.json", events, logs / f"reconstruct_{label}.log", env,
                          _verify_chunks(out, ("sim", "reco")), _chunk_progress(out, ("reco",), chunk)))
    if eic["analysis"]:
        plots = run_dir / "plots"
        dirs = " ".join(shlex.quote(str(d)) for ds in by_setting.values() for d in ds)
        settings = " ".join("--setting " + shlex.quote(f"{b}={','.join(map(str, ds))}") for b, ds in by_setting.items())
        tool = f"python3 {shlex.quote(str(EIC_ANALYSIS))}"
        target = shlex.quote(str(plots))
        jobs = min(16, os.cpu_count() or 1)
        steps.append(Step("analyze", container_command(eic["eic_shell"], (
            f"set -e; {tool} extract {dirs} -j {jobs}; {tool} plots {settings} -o {target} --lumi-fb {eic['lumi_fb']!r}; "
            f'for f in {target}/*.pdf; do pdftoppm -png -r 130 -singlefile "$f" "${{f%.pdf}}"; done')),
            plots / "summary.json", 1, logs / "analyze.log", env,
            lambda: None if (plots / "summary.json").is_file() else f"no figures in {plots}, see {logs}/analyze.log"))
    plan.jobs.append(steps)
    return plan


def eic_cost(cfg):
    """(events, disk GB, core hours) of a config, roughly, from the production runs."""
    eic = cfg["eic"]
    sample = schema.EIC_SAMPLES[eic["sample"]]
    events = gb = core_h = 0.0
    for _, beams, _, _, attempts in eic_slices(eic):
        n = attempts * sample.pass_fraction
        events += n
        gb += n * schema.EIC_MB_PER_EVENT / 1000
        core_h += n * schema.EIC_CORE_S_PER_EVENT[beams] / 3600
    return events, gb, core_h


def _existing_parent(path):
    path = Path(path)
    while not path.exists() and path != path.parent:
        path = path.parent
    return path


def _mem_available_gb():
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 2**20
    except (OSError, ValueError, IndexError):
        pass
    return None


def _container_load_check(key, install, setup, libraries):
    """Shell test, for inside the container: the install's libraries find everything they link to."""
    q = shlex.quote(install)
    return (f'( . {q}/{setup} >/dev/null 2>&1; for f in {q}/{libraries}; do [ -f "$f" ] || continue; '
            f'm=$(ldd "$f" 2>/dev/null | grep -m1 "not found"); if [ -n "$m" ]; then '
            f'echo "[eic] {key}: {install} was built against another image of the container and does not load in '
            f'this one (${{f##*/}}:$m); rebuild it inside eic-shell, or use the one in the container"; break; fi; done )')


def check_eic(cfg):
    """Everything a chain would otherwise only trip over hours later, as "[eic] key: problem" lines."""
    run, eic = cfg["run"], cfg["eic"]
    sample = schema.EIC_SAMPLES[eic["sample"]]
    problems = []
    if os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER"):
        return ["[eic] eic_shell: start gcf on the host, not inside eic-shell: the generators run on the "
                "host and the chain enters the container itself"]
    if not os.access(eic["eic_shell"], os.X_OK):
        problems.append(f"[eic] eic_shell: {eic['eic_shell']} is not an executable eic-shell launcher")
    if not Path(eic["fieldmaps"]).is_dir():
        problems.append(f"[eic] fieldmaps: {eic['fieldmaps']} is not a directory")
    if not os.access(REPO / EIC_GENERATOR_BINARIES[sample.generator], os.X_OK):
        problems.append(f"[eic] sample: {sample.generator} is not built ({EIC_GENERATOR_BINARIES[sample.generator]}); "
                        f"run ./install.sh")
    if eic["deexcite"] == "geant4" and not os.access(REPO / "tools/deexcite/build/deexcite", os.X_OK):
        problems.append("[eic] deexcite: tools/deexcite is not built; inside eic-shell run cmake -S tools/deexcite "
                        "-B tools/deexcite/build && cmake --build tools/deexcite/build")
    if eic["deexcite"] == "beagle":
        home = Path.home()
        beagle = {"BEAGLESYS": (home / "src/BeAGLE", "BeAGLE"), "FLUPRO": (home / "src/beagle-deps/fluka", "data/nuclear.bin"),
                  "LHAPATH": (home / "src/beagle-deps/lhapdf5/share/lhapdf/PDFsets", "cteq6ll.LHpdf")}
        for var, (default, member) in beagle.items():
            if not (Path(os.environ.get(var, default)) / member).exists():
                problems.append(f"[eic] deexcite: BeAGLE needs {member} in {var}={os.environ.get(var, default)}")
    # resources
    where = _existing_parent(run["output_dir"])
    _, gb, _ = eic_cost(cfg)
    free = shutil.disk_usage(where).free / 1e9
    if free < gb:
        problems.append(f"[eic] slices: these slices write about {gb:.0f} GB and {where} has {free:.0f} GB free; "
                        f"fewer attempts, or start gcf gui with --runs-dir on a larger disk")
    need = max(eic["sim_jobs"] * schema.EIC_GB_PER_SIM_JOB[eic["sim_threads"] == 1],
               eic["reco_jobs"] * schema.EIC_GB_PER_RECO_JOB)
    have = _mem_available_gb()
    if have is not None and need > have:
        problems.append(f"[eic] sim_jobs: the npsim or eicrecon jobs need about {need:.0f} GB of memory and "
                        f"{have:.0f} GB is available")
    if problems or not os.access(eic["eic_shell"], os.X_OK):
        return problems
    # inside the container: paths only it can see, and whether it can see ours
    checks = [f'[ -f {shlex.quote(eic["epic_install"])}/bin/thisepic.sh ] || echo "[eic] epic_install: no bin/thisepic.sh in {eic["epic_install"]}"',
              f'[ -f {shlex.quote(eic["eicrecon_install"])}/bin/eicrecon-this.sh ] || echo "[eic] eicrecon_install: no bin/eicrecon-this.sh in {eic["eicrecon_install"]}"',
              f'[ -d {shlex.quote(eic["fieldmaps"])} ] || echo "[eic] fieldmaps: {eic["fieldmaps"]} is not visible in the container"',
              f'[ -w {shlex.quote(str(where))} ] || echo "[eic] slices: {where} is not writable in the container: see the bind paths in eic-shell"',
              f'[ -f {shlex.quote(str(EIC_ANALYSIS))} ] || echo "[eic] eic_shell: this repository is not visible in the container"']
    share = shlex.quote(eic["epic_install"]) + "/share/epic"
    for beams in eic["beams"]:
        config, field = schema.EIC_BEAMS[beams]
        checks.append(f'[ -f {share}/{config}.xml ] || echo "[eic] epic_install: {eic["epic_install"]} has no {config}.xml for {beams}"')
        if field:
            checks.append(f'[ -f {share}/compact/{field} ] || echo "[eic] epic_install: {eic["epic_install"]} has no compact/{field} for {beams}"')
    # a build against an older image of the moving nightly fails only when eicrecon starts, after the simulation
    checks += [_container_load_check("epic_install", eic["epic_install"], "bin/thisepic.sh", "lib/libepic.so"),
               _container_load_check("eicrecon_install", eic["eicrecon_install"], "bin/eicrecon-this.sh",
                                     "lib/EICrecon/plugins/*.so")]
    try:
        out = subprocess.run(container_command(eic["eic_shell"], "; ".join(checks) + '; echo "[eic] probe done"'),
                             capture_output=True, text=True, timeout=180, env=eic_env(eic))
    except subprocess.TimeoutExpired:
        return ["[eic] eic_shell: the container did not answer within 3 minutes"]
    lines = [l for l in out.stdout.splitlines() if l.startswith("[eic] ")]
    if "[eic] probe done" not in lines:
        tail = (out.stderr + out.stdout).strip().splitlines()[-3:]
        return [f"[eic] eic_shell: could not run a command in the container: {' '.join(tail)}"]
    return [l for l in lines if l != "[eic] probe done"]


# ---------------------------------------------------------------------------
# Execution

class Runner:
    def __init__(self, plan, quiet=False):
        self.plan = plan
        self.quiet = quiet
        self.lock = threading.RLock()   # the SIGTERM handler may re-enter
        self.procs = set()
        self.cancelled = False
        self.state = "running"
        self.error = None
        self.started = time.time()
        self.finished = None
        self.status_path = plan.run_dir / "status.json"

    # -- status ------------------------------------------------------------
    def status(self):
        with self.lock:
            stages = {}
            for steps in self.plan.jobs:
                for s in steps:
                    st = stages.setdefault(s.stage, {"done": 0, "total": 0, "jobs": 0,
                                                     "finished": 0})
                    st["done"] += s.done
                    st["total"] += s.total
                    st["jobs"] += 1
                    st["finished"] += s.state == "done"
            return {
                "name": self.plan.name, "state": self.state, "error": self.error,
                "pid": os.getpid(), "pid_start": proc_start(os.getpid()),
                "started": self.started, "finished": self.finished,
                "stages": stages,
                "outputs": {k: str(v) for k, v in self.plan.outputs.items()},
                "jobs": [[{"stage": s.stage, "state": s.state, "done": s.done, "total": s.total,
                           "log": str(s.log), "pid": s.pid, "pid_start": s.pid_start}
                          for s in steps] for steps in self.plan.jobs],
            }

    def write_status(self):
        st = self.status()
        tmp = self.status_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(st, indent=1))
        tmp.replace(self.status_path)
        return st

    def print_progress(self, st, final=False):
        if self.quiet:
            return
        parts = [f"{k} {100.0 * v['done'] / max(v['total'], 1):5.1f}%"
                 for k, v in st["stages"].items()]
        line = f"[{self.plan.name}] " + " | ".join(parts) + f" | {time.time() - self.started:6.0f} s"
        if sys.stdout.isatty():
            print("\r" + line, end="\n" if final else "", flush=True)
        elif final or int(time.time() - self.started) % 30 == 0:
            print(line, flush=True)

    # -- jobs --------------------------------------------------------------
    def run_step(self, step):
        step.log.parent.mkdir(parents=True, exist_ok=True)
        step.output.parent.mkdir(parents=True, exist_ok=True)
        with self.lock:
            if self.cancelled:
                step.state = "cancelled"
                return False
            step.state = "running"
        with open(step.log, "w") as log:
            log.write("+ " + " ".join(step.args) + "\n")
            log.flush()
            proc = subprocess.Popen(step.args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, errors="replace",
                                    env=environment() if step.env is None else step.env,
                                    start_new_session=True)
            with self.lock:
                self.procs.add(proc)
                step.pid, step.pid_start = proc.pid, proc_start(proc.pid)
            for line in proc.stdout:
                log.write(line)
                if line.startswith("Progress:"):
                    try:
                        step.done = int(line.split()[1].split("/")[0])
                    except (IndexError, ValueError):
                        pass
            proc.wait()
            with self.lock:
                self.procs.discard(proc)
                step.pid = step.pid_start = None
        step.returncode = proc.returncode
        with self.lock:
            if self.cancelled:
                step.state = "cancelled"
                return False
            step.state = "done" if proc.returncode == 0 else "failed"
            if proc.returncode != 0 and self.error is None:
                tail = step.log.read_text(errors="replace").splitlines()[-8:]
                self.error = (f"{step.stage} job failed (exit {proc.returncode}), "
                              f"log {step.log}:\n" + "\n".join(tail))
                if step.context:
                    self.error += f"\n{step.context}:\n" + _tail(step.context, 8)
        if proc.returncode == 0 and step.verify is not None:
            problem = step.verify()
            if problem:
                with self.lock:
                    step.state = "failed"
                    if self.error is None:
                        self.error = f"{step.stage}: {problem}"
                return False
        if proc.returncode == 0:
            if step.poll is not None:
                step.poll(step)        # the true total, once the generator has written it
            step.done = step.total
        return proc.returncode == 0

    def run_job(self, steps):
        for step in steps:
            if not self.run_step(step):
                self.cancel(fail=True)
                return

    def cancel(self, fail=False):
        with self.lock:
            if self.cancelled:
                return
            self.cancelled = True
            if not fail:
                self.state = "cancelled"
            for p in list(self.procs):
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    def merge(self):
        for stage, inputs, output in self.plan.merges:
            if len(inputs) == 1:
                shutil.copyfile(inputs[0], output)
                continue
            log = self.plan.run_dir / "logs" / f"merge_{stage}.log"
            with open(log, "w") as f:
                rc = subprocess.run([hadd_command(), "-f", str(output)] + [str(p) for p in inputs],
                                    stdout=f, stderr=subprocess.STDOUT, env=environment()).returncode
            if rc != 0:
                raise RunError(f"merging the {stage} outputs failed, see {log}")

    def poll_steps(self):
        for steps in self.plan.jobs:
            for s in steps:
                if s.state == "running" and s.poll is not None:
                    try:
                        s.poll(s)
                    except OSError:
                        pass

    def run(self):
        self.plan.run_dir.mkdir(parents=True, exist_ok=True)
        self.write_status()
        threads = [threading.Thread(target=self.run_job, args=(steps,), daemon=True)
                   for steps in self.plan.jobs]
        for t in threads:
            t.start()
        ticks = 0
        try:
            while any(t.is_alive() for t in threads):
                if ticks % 5 == 0:     # the run directory may be on a slow network disk
                    self.poll_steps()
                ticks += 1
                self.print_progress(self.write_status())
                time.sleep(1.0)
        except KeyboardInterrupt:
            self.cancel()
            for t in threads:
                t.join()
        if not self.cancelled:
            try:
                self.state = "merging"
                self.write_status()
                self.merge()
                for f in self.plan.cleanup:
                    Path(f).unlink(missing_ok=True)
                self.state = "done"
            except RunError as e:
                self.state, self.error = "failed", str(e)
        elif self.state != "cancelled":
            self.state = "failed"
        self.finished = time.time()
        st = self.write_status()
        self.print_progress(st, final=True)
        return self.state == "done"


def execute(plan, quiet=False):
    """Run a plan in this process; SIGTERM cancels it cleanly."""
    runner = Runner(plan, quiet=quiet)
    signal.signal(signal.SIGTERM, lambda *_: runner.cancel())
    ok = runner.run()
    if not quiet:
        if ok:
            for stage, path in plan.outputs.items():
                print(f"  {stage:8s} -> {path}")
        else:
            print(f"Run {runner.state}." + (f"\n{runner.error}" if runner.error else ""),
                  file=sys.stderr)
    return ok
