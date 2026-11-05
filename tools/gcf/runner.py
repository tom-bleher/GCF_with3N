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
"""

import fcntl
import json
import math
import os
import queue
import re
import secrets
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
# Runs a stage, recording its process and its exit code next to its log, so that
# neither depends on the supervisor still being there:
# sh -c WRAPPER <log without its suffix> <command...>
STAGE_WRAPPER = 'echo $$ > "$0.pid"; "$@"; echo $? > "$0.exit"'
POLL_SECONDS = 0.5
# Without [run] chunks, a run has at least one job per this many events.
EVENTS_PER_CHECKPOINT = 1000000
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


def check_setup(cfg):
    """Run the stages' own setup checks on a validated config, without
    generating anything: gcf_generate --check (nucleus tables, option
    combinations) or, for FSI-only runs, that the input file and tree can be
    read. Returns a list of problems (empty = fine)."""
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
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
    def __init__(self, stage, args, output, total, log):
        self.stage, self.args, self.output, self.total, self.log = stage, args, output, total, log
        self.base = Path(log).with_suffix("")
        self.exit_file = Path(log).with_suffix(".exit")
        self.pid_file = Path(log).with_suffix(".pid")
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
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
    name = run["name"]
    run_dir = Path(run_dir)
    chunks, logs = run_dir / "chunks", run_dir / "logs"
    plan = Plan(run_dir, name)
    plan.workers = run["jobs"]
    seed = run["seed"]
    stages = run["stages"]
    final = {"generate": run_dir / f"{name}_pwia.root", "fsi": run_dir / f"{name}_fsi.root"}
    final.update({k: Path(v) for k, v in (final_outputs or {}).items()})
    for stage in stages:
        plan.outputs[stage] = final[stage]

    if "generate" in stages:
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
    if "generate" not in st.get("stages", {}):
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
        cfg["run"]["seed"] = 1 + secrets.randbelow(2**31 - 2**21)
    return cfg


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
            if line.startswith("Progress:"):
                try:
                    step.done = int(line.split()[1].split("/")[0])
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
                                            env=environment(), start_new_session=True)
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
        if code == 0:
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
        try:
            while any(t.is_alive() for t in threads):
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
