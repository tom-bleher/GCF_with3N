"""Running the pipeline stages: environment, parallel jobs, merging, status.

A run is split into jobs. Each job runs its steps (gcf_generate and/or
gcf_fsi) one after the other as subprocesses; jobs run in parallel. Their
outputs are then merged with hadd. Progress is parsed from the stages'
"Progress: n/N" lines and written to <run_dir>/status.json, which the GUI
reads.
"""

import json
import os
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
GENIE_DIR = REPO / "src/programs/genQE_FSI/Generator-R-3_06_02"
FSI_CONFIG = REPO / "src/programs/genQE_FSI/config"


class RunError(Exception):
    pass


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
    if first is not None:
        args += ["--first-event", str(first), "--max-events", str(count)]
    if not copy_objects:
        args.append("--no-copy-objects")
    return args


def count_entries(path):
    out = subprocess.run([binary("gcf_fsi"), "--input", str(path), "--count"],
                         capture_output=True, text=True, env=environment())
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
        self.done = 0
        self.state = "pending"       # pending, running, done, failed, cancelled
        self.returncode = None


class Plan:
    def __init__(self, run_dir, name):
        self.run_dir = Path(run_dir)
        self.name = name
        self.jobs = []               # list of lists of Step
        self.merges = []             # (stage, [inputs], output)
        self.outputs = {}            # stage -> final file
        self.cleanup = []            # files to delete after a successful run


def plan_from_config(cfg, run_dir, final_outputs=None):
    """Jobs and merges for a validated config.

    Final files are <run_dir>/<name>_pwia.root and <name>_fsi.root unless
    final_outputs maps a stage to another path. With one job the stages
    write their final files directly; with several, per-job files go to
    <run_dir>/chunks and are merged with hadd.
    """
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
        cfg["run"]["seed"] = 1 + secrets.randbelow(2**31 - 2**21)
    return cfg


# ---------------------------------------------------------------------------
# Execution

class Runner:
    def __init__(self, plan, quiet=False):
        self.plan = plan
        self.quiet = quiet
        self.lock = threading.Lock()
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
                "pid": os.getpid(), "started": self.started, "finished": self.finished,
                "stages": stages,
                "outputs": {k: str(v) for k, v in self.plan.outputs.items()},
                "jobs": [[{"stage": s.stage, "state": s.state, "done": s.done, "total": s.total,
                           "log": str(s.log)} for s in steps] for steps in self.plan.jobs],
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
                                    text=True, errors="replace", env=environment(),
                                    start_new_session=True)
            with self.lock:
                self.procs.add(proc)
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
        if proc.returncode == 0:
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

    def run(self):
        self.plan.run_dir.mkdir(parents=True, exist_ok=True)
        self.write_status()
        threads = [threading.Thread(target=self.run_job, args=(steps,), daemon=True)
                   for steps in self.plan.jobs]
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
