"""Running the pipeline stages: environment, parallel jobs, merging, status.

A run is split into jobs. Each job runs its steps (gcf_generate and/or
gcf_fsi) one after the other as subprocesses; `[run] jobs` of them run at a
time. Their outputs are then merged with hadd. Progress is parsed from the
stages' "Progress: n/N" lines and written to <run_dir>/status.json, which the
GUI reads.

A run is meant to outlive what started it. The process that supervises it and
every stage process are in sessions of their own, so closing a terminal or a
browser doesn't touch them. The stages write to their log files directly and
record their exit codes next to them, so they also outlive the supervisor: a new
one (`gcf run --resume`) takes over the stage processes that are still running
and reruns only the jobs that didn't finish. Long runs are split into more jobs
than run at a time, so that an interruption costs little (`checkpoint_jobs`).

The eic stage is one job: tools/run_eic_sim.sh for each beam setting and Q^2
slice, one step per stage of that script, then the analysis. That script hides
failures inside the container from its exit code, so each step checks what it
was supposed to leave behind, and progress is counted from the chunks' files.
"""

import fcntl
import json
import math
import os
import queue
import re
import secrets
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

from . import schema

REPO = Path(__file__).resolve().parents[2]
BIN_DIR = Path(os.environ.get("GCF_BIN_DIR", REPO / "src/programs/gcf_stages/build"))
# Runs a stage, recording its process and its exit code next to its log, so that
# neither depends on the supervisor still being there:
# sh -c WRAPPER <log without its suffix> <command...>
STAGE_WRAPPER = 'echo $$ > "$0.pid"; "$@"; echo $? > "$0.exit"'
POLL_SECONDS = 0.5
# Without [run] chunks, a run has at least one job per this many events.
EVENTS_PER_CHECKPOINT = 1000000
SUITE_DIR = Path(os.environ.get("GCF_SUITE_DIR", REPO / "build/programs"))
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
        if fields[0] == "Z":          # exited, not yet collected by its parent
            return None
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


FAILURE_HINTS = [
    (r"No space left|Disk quota", "The disk is full. Free some space, then resume the run: it keeps what finished."),
    (r"killed, for example|Killed|out of memory|bad_alloc", "It may have run out of memory. Resume the run with fewer "
     "parallel jobs in its config, or wait until the machine is less busy."),
    (r"error while loading shared libraries|cannot open shared object", "A program can't find its libraries. "
     "./gcf doctor shows what is missing, and ./install.sh usually fixes it."),
    (r"GENIE is not set|GXMLPATH", "The GENIE setup is missing. ./gcf doctor shows what is wrong."),
    (r"Permission denied", "A file or directory can't be written. Check the permissions of the runs directory."),
]


def failure_hint(error):
    """What someone can do about a failed run, from its error text."""
    for pattern, hint in FAILURE_HINTS:
        if error and re.search(pattern, error):
            return hint
    return None


def default_runs_dir():
    """Where runs go without --runs-dir: $GCF_RUNS_DIR, <repo>/runs, or, in a
    checkout you can't write to, ~/gcf_runs."""
    env = os.environ.get("GCF_RUNS_DIR")
    if env:
        return Path(env).expanduser()
    repo_runs = REPO / "runs"
    if os.access(repo_runs if repo_runs.exists() else REPO, os.W_OK):
        return repo_runs
    return Path.home() / "gcf_runs"


def hadd_command():
    hadd = shutil.which("hadd", path=environment().get("PATH"))
    if not hadd:
        raise RunError("hadd not found: ROOT is not set up (build/gcf_env.sh or setup_env.sh)")
    return hadd


# ---------------------------------------------------------------------------
# Command lines

def generate_args(gen, events, output, seed, base_seed=None):
    """gcf_generate's command line for one job. `base_seed` is the run's: the
    options that must be the same in every job of a run derive from it."""
    args = [binary("gcf_generate"), "--model", gen["model"], "--beam", repr(gen["beam"]),
            "--events", str(events), "--output", str(output), "--seed", str(seed),
            "-A", str(gen["target"][0]), "-Z", str(gen["target"][1])]
    for o in schema.GENERATE:
        if o.key in ("model", "beam", "events", "target") or o.flag is None:
            continue
        v = gen.get(o.key)
        if v is None or gen["model"] not in o.models:
            continue
        if o.key == "randomize":
            if v:
                args += [o.flag, str(base_seed or seed)]
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
    if fsi.get("tree"):
        args += ["--tree", fsi["tree"]]
    if fsi.get("overwrite"):
        args.append("--overwrite")
    for c in fsi.get("off") or []:
        args += ["--channel", f"{c}=0"]
    for c, factor in (fsi.get("scale") or {}).items():
        args += ["--channel", f"{c}={factor}"]
    if fsi.get("nucleon_mfp_scale") is not None:
        args += ["--channel", f"nucleon.mfp={fsi['nucleon_mfp_scale']}"]
    if fsi.get("pion_mfp_scale") is not None:
        args += ["--channel", f"pion.mfp={fsi['pion_mfp_scale']}"]
    if first is not None:
        args += ["--first-event", str(first), "--max-events", str(count)]
    if not copy_objects:
        args.append("--no-copy-objects")
    return args


_has_seed = {}


def program_binary(name):
    path = SUITE_DIR / ("sfSampler" if name == "sfSamplerLC" else name) / name
    if not os.access(path, os.X_OK):
        raise RunError(f"{name} is not built ({path} missing). Build it with ./install.sh, "
                       f"or set GCF_SUITE_DIR to where the generators are.")
    return str(path)


def program_has_seed(name):
    """Whether a generator takes -S yet (#137, PR #150): its usage says."""
    if name not in _has_seed:
        out = subprocess.run([program_binary(name), "-h"], capture_output=True, text=True, env=environment(), timeout=60)
        _has_seed[name] = "-S:" in (out.stdout + out.stderr)
    return _has_seed[name]


def program_args(prog, events, output, seed):
    """The command line of one of the suite's other generators, for one job."""
    from .schema import PROGRAMS, PROGRAM
    name = prog["program"]
    a, z = prog["target"]
    args = [program_binary(name), str(z), str(a - z)]
    if PROGRAMS[name].beam:
        args.append(repr(prog["beam"]))
    if PROGRAMS[name].pdf:
        args.append(prog.get("pdf") or str(REPO / "src/libraries/cross_sections/CTEQ/ct18nn.00.pds"))
    args += [str(output), str(events), "-v"]   # -v: what it prints is its progress
    for o in PROGRAM:
        v = prog.get(o.key)
        if o.flag is None or v is None or name not in o.models:
            continue
        if o.kind == "bool":
            if v:
                args.append(o.flag)
        else:
            args += [o.flag, str(v)]
    if program_has_seed(name):
        args += ["-S", str(seed)]
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
    if "program" in run["stages"]:
        with tempfile.TemporaryDirectory(prefix="gcf-check-") as tmp:
            try:
                args = program_args(cfg["program"], 200, Path(tmp) / "check.root", 1)
                out = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120, stdin=subprocess.DEVNULL)
            except RunError as e:
                return [str(e)]
            except subprocess.TimeoutExpired:
                return ["[program] program: a test of 200 attempts took more than 2 minutes"]
        if out.returncode != 0:
            text = [l for l in (out.stderr + out.stdout).splitlines() if l.strip() and not l.startswith(("eNCrossSection", "You are using", "Working on", "Random seed"))]
            problems.append("[program] target: the generator can't run this setup: " + (" ".join(text[-2:]) or f"exit {out.returncode}"))
        return problems
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
            # Point at the field that causes it.
            key = "nn" if "does not have contacts" in text else \
                  "target" if re.search(r"nucleus|nuclei|residual|mass", text) else "model"
            problems.append(f"[generate] {key}: the generator can't run this setup: " + text)
    else:
        problems += check_fsi_input(fsi, env)
    if "fsi" in run["stages"] and not (GENIE_DIR / "lib").is_dir():
        problems.append("GENIE is not built (no lib/ in the vendored tree); run ./install.sh")
    return problems


def check_fsi_input(fsi, env):
    """gcf_fsi --check on an FSI-only run's input: can it read the file and tree,
    has it the branches, has it been through FSI already. Problems as a list."""
    args = [binary("gcf_fsi"), "--input", str(fsi["input"]), "--check", "--fsi-model", fsi["model"]]
    if fsi.get("target"):
        args += ["--target", f"{fsi['target'][0]}:{fsi['target'][1]}"]
    if fsi.get("tree"):
        args += ["--tree", fsi["tree"]]
    if fsi.get("overwrite"):
        args.append("--overwrite")
    try:
        out = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120)
    except subprocess.TimeoutExpired:
        return [f"reading {fsi['input']} took more than 2 minutes"]
    if out.returncode == 0:
        return []
    errors = [l[len("Error: "):] for l in out.stderr.splitlines() if l.startswith("Error: ")]
    if any("already went through FSI" in e for e in errors):
        return ["[fsi] input: this file has already been through FSI. Choose the *_pwia.root file, "
                "or tick \"Replace earlier FSI in the input\" under Advanced."]
    return [f"[fsi] input: {errors[0] if errors else out.stderr.strip()[-300:] or f'exit {out.returncode}'}"]


def nuclear_params(model, target, nn=None):
    """The nuclear inputs of a setup, from the generator itself (the tables are C++)."""
    args = [binary("gcf_generate"), "--print-params", "--model", model, "--beam", "6",
            "-A", str(target[0]), "-Z", str(target[1])]
    if nn:
        args += ["--nn", nn]
    out = subprocess.run(args, capture_output=True, text=True, env=environment(), timeout=120)
    for line in out.stdout.splitlines():
        if line.startswith("PARAMS "):
            return json.loads(line[len("PARAMS "):])
    text = [l for l in (out.stderr + out.stdout).splitlines()
            if l.strip() and not l.startswith(("eNCrossSection", "You are using"))]
    raise RunError(" ".join(text[-2:]) or f"gcf_generate --print-params: exit {out.returncode}")


def count_entries(path, tree=None):
    try:
        out = subprocess.run([binary("gcf_fsi"), "--input", str(path), "--count"] + (["--tree", tree] if tree else []),
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
        self.base = Path(log).with_suffix("")
        self.exit_file = Path(log).with_suffix(".exit")
        self.pid_file = Path(log).with_suffix(".pid")
        self.env = env               # None: environment()
        self.verify = verify         # after exit 0: returns an error message, or None if all is well
        self.poll = poll             # while running: updates self.done / self.total from files
        self.context = context       # a log worth quoting when the step fails, besides its own
        self.done = 0
        self.state = "pending"       # pending, running, done, failed, cancelled
        self.returncode = None
        self.pid = None              # of the stage process while it runs ...
        self.pid_start = None        # ... and its start time, to identify it

    def live_pid(self):
        """The process running this step, whoever started it: its record says which
        one, and its command line shows it's still that one, not a reused number."""
        try:
            pid = int(self.pid_file.read_text().strip())
            with open(f"/proc/{pid}/cmdline", "rb") as f:
                arguments = f.read().split(b"\0")
        except (OSError, ValueError):
            return None
        if str(self.base).encode() in arguments and proc_start(pid) is not None:
            return pid
        return None

    def exit_code(self):
        try:
            return int(self.exit_file.read_text().strip())
        except (OSError, ValueError):
            return None


class Plan:
    def __init__(self, run_dir, name):
        self.run_dir = Path(run_dir)
        self.name = name
        self.jobs = []               # list of lists of Step
        self.workers = 1             # how many jobs run at a time
        self.merges = []             # (stage, [inputs], output)
        self.outputs = {}            # stage -> final file
        self.cleanup = []            # files to delete after a successful run
        self.eic = False             # the ePIC chain: its steps poll for progress


def checkpoint_jobs(run, events):
    """How many jobs a run is split into: [run] chunks, else one per
    EVENTS_PER_CHECKPOINT events and at least what runs at a time."""
    if run.get("chunks"):
        return run["chunks"]
    return max(run["jobs"], math.ceil(events / EVENTS_PER_CHECKPOINT))


def plan_from_config(cfg, run_dir, final_outputs=None, resume=False):
    """Jobs and merges for a validated config.

    Final files are <run_dir>/<name>_pwia.root and <name>_fsi.root unless
    final_outputs maps a stage to another path. With one job the stages
    write their final files directly; with several, per-job files go to
    <run_dir>/chunks and are merged with hadd.

    With `resume`, the steps that finished in the run's last status.json are
    kept, and those still running are taken over.
    """
    if "eic" in cfg["run"]["stages"]:
        return eic_plan(cfg, run_dir, resume)
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
    name = run["name"]
    run_dir = Path(run_dir)
    chunks, logs = run_dir / "chunks", run_dir / "logs"
    plan = Plan(run_dir, name)
    plan.workers = run["jobs"]
    seed = run["seed"]
    stages = run["stages"]
    final = {"generate": run_dir / f"{name}_pwia.root", "fsi": run_dir / f"{name}_fsi.root",
             "program": run_dir / f"{name}_gen.root"}
    final.update({k: Path(v) for k, v in (final_outputs or {}).items()})
    for stage in stages:
        plan.outputs[stage] = final[stage]

    if "program" in stages:
        prog = cfg["program"]
        parts = split(prog["events"], checkpoint_jobs(run, prog["events"]))
        single = len(parts) == 1
        outs = {"program": []}
        for i, n in enumerate(parts):
            g = final["program"] if single else chunks / f"gen_{i:03d}.root"
            plan.jobs.append([Step("program", program_args(prog, n, g, seed + i), g, n, logs / f"program_{i:03d}.log")])
            outs["program"].append(g)
    elif "generate" in stages:
        parts = split(gen["events"], checkpoint_jobs(run, gen["events"]))
        single = len(parts) == 1
        outs = {"generate": [], "fsi": []}
        for i, n in enumerate(parts):
            g = final["generate"] if single else chunks / f"pwia_{i:03d}.root"
            steps = [Step("generate", generate_args(gen, n, g, seed + i, seed), g, n,
                          logs / f"generate_{i:03d}.log")]
            outs["generate"].append(g)
            if "fsi" in stages:
                f = final["fsi"] if single else chunks / f"fsi_{i:03d}.root"
                steps.append(Step("fsi", fsi_args(fsi, g, f, seed + 1000003 + i), f, n,
                                  logs / f"fsi_{i:03d}.log"))
                outs["fsi"].append(f)
            plan.jobs.append(steps)
    else:
        total = count_entries(fsi["input"], fsi.get("tree"))
        if total == 0:
            raise RunError(f"{fsi['input']} has no events")
        parts = split(total, checkpoint_jobs(run, total))
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
    if resume:
        restore_steps(plan)
    return plan


def read_status(run_dir):
    try:
        return json.loads((Path(run_dir) / "status.json").read_text())
    except (OSError, ValueError):
        return {}


def restore_steps(plan):
    """Take up a plan where it stopped. Each step's own records decide, not the
    supervisor's last status, which can be a moment behind: a step still running is
    taken over, one that exited with 0 is kept, and any other runs again."""
    old = read_status(plan.run_dir).get("jobs", [])
    for i, steps in enumerate(plan.jobs):
        for j, step in enumerate(steps):
            has_output = step.output.is_file() and step.output.stat().st_size > 0
            pid = step.live_pid()
            if pid:
                step.pid, step.pid_start = pid, proc_start(pid)
                prev = old[i][j] if i < len(old) and j < len(old[i]) else {}
                step.done = prev.get("done", 0)
            elif step.exit_code() == 0 and has_output:
                step.state, step.done = "done", step.total
            elif not step.exit_file.exists() and not step.pid_file.exists() and has_output:
                # Finished under an older version, which kept no records of its own.
                prev = old[i][j] if i < len(old) and j < len(old[i]) else {}
                if prev.get("state") == "done" and prev.get("stage") == step.stage:
                    step.state, step.done = "done", step.total


def merge_finished(run_dir):
    """Make final files from the jobs that finished in a run that stopped. Their
    yield is right: each job's file carries its own attempts, and hadd adds them."""
    run_dir = Path(run_dir)
    st = read_status(run_dir)
    if not st.get("jobs"):
        raise RunError("this run has no jobs to merge")
    if "simulate" in st.get("stages", {}):
        raise RunError("an ePIC run has no files to merge: resume it, and its finished slices and chunks stay")
    if not {"generate", "program"} & set(st.get("stages", {})):
        raise RunError("an FSI-only run can't be merged in part: its first job carries the attempts "
                       "of the whole input file. Resume it instead")
    partial, outputs = {}, st.get("outputs", {})
    for stage, final in outputs.items():
        steps = [s for job in st["jobs"] for s in job if s["stage"] == stage]
        done = [s for s in steps if s["state"] == "done" and Path(s["output"]).is_file()]
        if not done:
            continue
        if done[0]["output"] != final:
            inputs = [s["output"] for s in done]
            log = run_dir / "logs" / f"merge_{stage}.log"
            log.parent.mkdir(exist_ok=True)
            with open(log, "w") as f:
                rc = subprocess.run([hadd_command(), "-f", final] + inputs, stdout=f, stderr=subprocess.STDOUT,
                                    env=environment()).returncode
            if rc != 0:
                raise RunError(f"merging the {stage} outputs failed, see {log}")
        partial[stage] = {"events": sum(s["total"] for s in done), "total": sum(s["total"] for s in steps),
                          "jobs": len(done), "jobs_total": len(steps)}
    if not partial:
        raise RunError("no job of this run has finished yet")
    st.update(state="partial", partial=partial, finished=time.time())
    tmp = run_dir / "status.json.tmp"
    tmp.write_text(json.dumps(st, indent=1))
    tmp.replace(run_dir / "status.json")
    return partial


class RunLock:
    """Held by a run's supervisor for as long as it lives, so that two can't work on
    a run. The kernel releases it when the process dies, however that happens."""

    def __init__(self, run_dir):
        Path(run_dir).mkdir(parents=True, exist_ok=True)
        self.file = open(Path(run_dir) / ".lock", "w")
        try:
            fcntl.flock(self.file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.file.close()
            raise RunError(f"{run_dir} is being run by another process") from None


def is_locked(run_dir):
    """True if a supervisor is working on this run."""
    try:
        with open(Path(run_dir) / ".lock", "w") as f:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(f, fcntl.LOCK_UN)
        return False
    except OSError:
        return True


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


def stop_leftover_steps(run_dir):
    """Stop what is still running of an earlier supervisor's steps. A resume of the eic
    stage plans other steps than the ones it continues, so it can't take them over, and
    two copies must not work on one set of chunks."""
    live = []
    for record in (Path(run_dir) / "logs").glob("*.pid"):
        step = Step("", [], Path(run_dir), 0, record.with_suffix(".log"))
        pid = step.live_pid()
        if pid:
            live.append(pid)
            try:
                os.killpg(pid, signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass
    deadline = time.time() + 20
    while time.time() < deadline and any(proc_start(p) is not None for p in live):
        time.sleep(0.2)


def eic_plan(cfg, run_dir, resume=False):
    """One job: for each beam setting and Q^2 slice generate, simulate, reconstruct; then the figures.

    With resume, finished slices are skipped and a slice that has its chunks restarts with run_eic_sim.sh
    --stage resume, which runs only the chunks that have not succeeded; one stopped earlier starts over, and
    regenerates the same events since the seed is in the saved config.
    """
    run, eic = cfg["run"], cfg["eic"]
    run_dir = Path(run_dir)
    if resume:
        stop_leftover_steps(run_dir)
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
        self.running = {}            # pid -> Step of the stage processes
        self.cancelled = False
        self.state = "running"
        self.error = None
        self.started = time.time()
        self.finished = None
        self.status_path = plan.run_dir / "status.json"
        self.queue = queue.SimpleQueue()
        previous = read_status(plan.run_dir) if any(s.state == "done" or s.pid for j in plan.jobs for s in j) else {}
        self.started = previous.get("started") or self.started

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
                "workers": self.plan.workers,
                "jobs": [[{"stage": s.stage, "state": s.state, "done": s.done, "total": s.total,
                           "log": str(s.log), "output": str(s.output), "pid": s.pid, "pid_start": s.pid_start}
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
        try:
            if sys.stdout.isatty():
                print("\r" + line, end="\n" if final else "", flush=True)
            elif final or int(time.time() - self.started) % 30 == 0:
                print(line, flush=True)
        except OSError:
            self.quiet = True      # the terminal is gone; the run carries on

    # -- jobs --------------------------------------------------------------
    @staticmethod
    def read_progress(step, offset, partial):
        """Take the new "Progress: n/N" lines of a step's log; returns where to go on from."""
        try:
            with open(step.log, "rb") as f:
                f.seek(offset)
                data = f.read()
        except OSError:
            return offset, partial
        lines = (partial + data.decode(errors="replace")).split("\n")
        for line in lines[:-1]:
            try:
                if line.startswith("Progress:"):
                    step.done = int(line.split()[1].split("/")[0])
                elif line.startswith("Working on event "):      # the suite's other generators
                    step.done = int(line.split()[3])
            except (IndexError, ValueError):
                pass
        return offset + len(data), lines[-1]

    def run_step(self, step):
        if step.state == "done":
            return True
        step.log.parent.mkdir(parents=True, exist_ok=True)
        step.output.parent.mkdir(parents=True, exist_ok=True)
        proc = None
        with self.lock:
            if self.cancelled:
                step.state = "cancelled"
                return False
            step.state = "running"
            if step.pid is None:     # otherwise a previous supervisor started it
                leftover = step.live_pid()
                if leftover:         # of a run being overwritten: two can't write one file
                    try:
                        os.killpg(leftover, signal.SIGTERM)
                    except (ProcessLookupError, PermissionError):
                        pass
                step.exit_file.unlink(missing_ok=True)
                step.pid_file.unlink(missing_ok=True)
                with open(step.log, "w") as log:
                    log.write("+ " + " ".join(step.args) + "\n")
                    log.flush()
                    proc = subprocess.Popen(["/bin/sh", "-c", STAGE_WRAPPER, str(step.base), *step.args],
                                            stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                            env=environment() if step.env is None else step.env,
                                            start_new_session=True)
                step.pid, step.pid_start = proc.pid, proc_start(proc.pid)
            self.running[step.pid] = step
        offset, partial = 0, ""
        while True:
            alive = proc.poll() is None if proc else proc_is(step.pid, step.pid_start)
            offset, partial = self.read_progress(step, offset, partial)
            if not alive:
                break
            time.sleep(POLL_SECONDS)
        code = step.exit_code()      # None: killed before it could record one
        with self.lock:
            self.running.pop(step.pid, None)
            step.pid = step.pid_start = None
            step.returncode = code
            if self.cancelled:
                step.state = "cancelled"
                return False
            step.state = "done" if code == 0 else "failed"
            if code != 0 and self.error is None:
                tail = step.log.read_text(errors="replace").splitlines()[-8:]
                why = f"exit {code}" if code is not None else "killed, for example because the machine ran out of memory"
                self.error = f"{step.stage} job failed ({why}), log {step.log}:\n" + "\n".join(tail)
                if step.context:
                    self.error += f"\n{step.context}:\n" + _tail(step.context, 8)
        if code == 0 and step.verify is not None:
            problem = step.verify()
            if problem:
                with self.lock:
                    step.state = "failed"
                    if self.error is None:
                        self.error = f"{step.stage}: {problem}"
                return False
        if code == 0:
            if step.poll is not None:
                step.poll(step)        # the true total, once the generator has written it
            step.done = step.total
        return code == 0

    def run_job(self, steps):
        for step in steps:
            if not self.run_step(step):
                return

    def worker(self):
        while True:
            try:
                steps = self.queue.get_nowait()
            except queue.Empty:
                return
            self.run_job(steps)

    def cancel(self, fail=False):
        with self.lock:
            if self.cancelled:
                return
            self.cancelled = True
            if not fail:
                self.state = "cancelled"
            for pid in list(self.running):
                try:
                    os.killpg(pid, signal.SIGTERM)
                except (ProcessLookupError, PermissionError):
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
        # Stage processes taken over from a previous supervisor come first, and each
        # needs a thread to follow it, whatever the number of jobs at a time.
        pending = sorted((steps for steps in self.plan.jobs if any(s.state != "done" for s in steps)),
                         key=lambda steps: not any(s.pid for s in steps))
        adopted = sum(1 for steps in pending if any(s.pid for s in steps))
        for steps in pending:
            self.queue.put(steps)
        threads = [threading.Thread(target=self.worker, daemon=True)
                   for _ in range(max(self.plan.workers, adopted) if pending else 0)]
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
        failed = [s for steps in self.plan.jobs for s in steps if s.state == "failed"]
        if self.cancelled:
            if self.state != "cancelled":
                self.state = "failed"
        elif failed:
            self.state = "failed"
            finished = sum(all(s.state == "done" for s in steps) for steps in self.plan.jobs)
            self.error = (f"{len(failed)} of {len(self.plan.jobs)} jobs failed; {finished} finished and are kept, "
                          f"so a resume reruns only the rest. The first failure:\n{self.error}")
        else:
            try:
                self.state = "merging"
                self.write_status()
                self.merge()
                for f in self.plan.cleanup:
                    Path(f).unlink(missing_ok=True)
                self.state = "done"
            except RunError as e:
                self.state, self.error = "failed", str(e)
        self.finished = time.time()
        st = self.write_status()
        self.print_progress(st, final=True)
        return self.state == "done"


def execute(plan, quiet=False):
    """Run a plan in this process; SIGTERM cancels it cleanly."""
    runner = Runner(plan, quiet=quiet)
    signal.signal(signal.SIGTERM, lambda *_: runner.cancel())
    # A closed terminal is not a reason to stop.
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    ok = runner.run()
    if not quiet:
        if ok:
            for stage, path in plan.outputs.items():
                print(f"  {stage:8s} -> {path}")
        else:
            try:
                print(f"Run {runner.state}." + (f"\n{runner.error}" if runner.error else ""),
                      file=sys.stderr)
            except OSError:
                pass
    return ok
