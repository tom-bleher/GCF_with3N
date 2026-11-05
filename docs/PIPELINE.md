# GCF_with3N pipeline: generation and FSI as separate stages

This page is about running things. What the generators compute (GCF
assumptions, the expressions in the weights, how FSI is done) is in
[PHYSICS.md](PHYSICS.md).

The generators used to do everything in one run: `genQE_3N_FSI` / `genQE_FSI`
generate GCF quasi-elastic events and immediately transport the outgoing
nucleons through GENIE. The pipeline splits this into independent stages that
talk through one event-file format, so you can run only the part you need:

| Stage | Program | In | Out |
|---|---|---|---|
| 1. Generate (PWIA) | `gcf_generate` | settings | event file |
| 2. FSI | `gcf_fsi` | **any** event file in the format below | same file + the post-FSI final state |

Stage 2 does not care who wrote its input: events from another generator get
the same GENIE treatment, as long as the file follows the format. Later stages
(radiative corrections, acceptance maps) are meant to follow the same pattern:
event file in, event file out ([Adding a stage](#adding-a-stage)).

The existing combined programs are unchanged and still work.

## Contents
- [Quick start](#quick-start)
- [The `gcf` command](#the-gcf-command)
- [The web interface](#the-web-interface)
- [Event-file format](#event-file-format)
- [FSI output](#fsi-output)
- [Weights and normalization](#weights-and-normalization)
- [Seeds and provenance](#seeds-and-provenance)
- [Differences from the combined programs](#differences-from-the-combined-programs)
- [Validation](#validation)
- [Adding a stage](#adding-a-stage)

## Quick start

```bash
./install.sh                         # once: builds GENIE, the generators and the stages
./gcf example 3N > my_run.toml       # an annotated example config (3N, 2N, MF or fsi)
./gcf run my_run.toml                # -> runs/<name>/<name>_pwia.root and <name>_fsi.root
./gcf quicklook runs/<name>/<name>_fsi.root    # a page of control plots (PNG)
```

`./gcf` sets up ROOT/GENIE by itself (from `build/gcf_env.sh`, written by
`install.sh`, or `setup_env.sh`), so there is nothing to source first.

## The `gcf` command

```text
gcf run CONFIG.toml        run the stages a config describes (--detach, --resume)
gcf scan CONFIG.toml --set generate.sigma_cm=0.1,0.15   the config for several values
gcf status                 the runs in a directory and how far they are
gcf doctor                 check that this machine is ready, and say how to fix what isn't
gcf generate [options]     PWIA events -> event file           (gcf_generate)
gcf fsi [options]          event file -> + GENIE FSI           (gcf_fsi)
gcf quicklook FILE         control plots of an event file, with an optional selection
gcf summary FILE           the numbers of an event file: N_eff, sum(w)/attempts, provenance
gcf gui                    the web interface
gcf example [3N|2N|MF|fsi] print an example config
gcf env                    print the environment the stages run with
```

Every subcommand has `--help`. All runs can be split into parallel jobs
(`[run] jobs` or `-j N`): each job gets its own seed, and the outputs are
merged with `hadd`.

### Config files

A config has up to three tables. Unset options take the defaults shown by
`gcf <command> --help`; options that don't apply to the chosen model are an
error rather than silently ignored, both here and in `gcf_generate` itself
(`src/programs/gcf_stages/test_model_options.py`, part of the `install.sh`
smoke test, checks that the two agree).

```toml
[run]
name = "c12_3N_6GeV"          # outputs: <output_dir>/<name>/<name>_pwia.root, _fsi.root
output_dir = "runs"           # relative to this file
stages = ["generate", "fsi"]  # or ["generate"], or ["fsi"] with [fsi] input
jobs = 16
seed = 0                      # 0: draw one; the run's saved config records it

[generate]
model = "3N"                  # 3N, 2N or MF
beam = 6.0                    # GeV
events = 1000000              # accepted events, over all jobs
target = [12, 6]
ff = "kelly"                  # nucleon form factors: kelly, ye or dipole
wave_function = "AV18"        # 3N: AV18, AV18+UIX, AV18+TM, CDBonn, CDBonn+TM, or ...-pD
theta_e = [5.0, 40.0]         # deg
# 2N: nn = "AV18" | "N2LO" | "N2LO12" | "N3LO" | "AV4" | "NV"; MF (12C only): fermi_gas = "global" | "local"
# 2N nuclear inputs, default from the tables (gcf_generate --print-params shows them):
# cpp0, cnn0, cpn0, cpn1 (contacts); sigma_cm; estar; sigma_e; randomize = true; light_cone = true

[fsi]
model = "hN"                  # hN or hA
transport = "shared"          # or "independent"
off = ["pion.absorption"]     # cascade channels to switch off (see below)
scale = { "pion.cex" = 0.8 }  # or to scale by a factor
```

Every run directory holds the resolved `config.toml` (including the seed),
per-job logs and a `status.json`; `gcf run runs/<name>/config.toml` repeats it.

### Interrupted runs

A run has a supervisor process (`gcf run`) and stage processes, each in a session
of its own, so a closed terminal doesn't reach them. The stages write straight to
their logs and record their exit codes next to them, so they also outlive the
supervisor. `gcf run --resume <run>/config.toml` starts a new one: it keeps the
jobs that finished, takes over stage processes that are still running, and runs
the others again with the same seeds. The result is the same events as a run that
was never interrupted. A file lock stops two supervisors working on one run.

`[run] jobs` is how many jobs run at a time. `[run] chunks` is how many jobs the
run is split into; by default one per million events, and at least `jobs`. A
failed job doesn't stop the others, and a resume reruns only it. `gcf run
--merge-finished <run>/config.toml` makes output files from the jobs that finished;
each job's file carries its own `hAttempts`, so Σw/attempts stays right. It isn't
available for FSI-only runs, whose first job carries the whole input's attempts.

### Interaction channels of the cascade

`gcf_fsi --channel NAME=FACTOR` scales the rate of one interaction channel
of the cascade by FACTOR; 0 switches it off and the remaining fates are
renormalized. In a config, `[fsi] off = [...]` lists the channels to switch
off and `[fsi] scale = { "NAME" = FACTOR }` scales them; `gcf fsi --off NAME`
and `gcf fsi --scale NAME=FACTOR` do the same on the command line. Which channels
exist depends on the model (`gcf_fsi --list-channels`):

| channel | hN | hA | |
|---|---|---|---|
| `nucleon.elastic` | yes | – | elastic is not a separate hA fate |
| `nucleon.inelastic` | yes | yes | |
| `nucleon.cex` | – | yes | |
| `nucleon.absorption` | – | yes | |
| `nucleon.piprod` | – | yes | |
| `nucleon.compound` | yes | yes | compound-nucleus (pre-equilibrium) fate |
| `pion.elastic` | yes | – | |
| `pion.inelastic` | yes | yes | |
| `pion.cex` | yes | yes | |
| `pion.absorption` | yes | yes | |
| `pion.piprod` | – | yes | |
| `nucleon.mfp`, `pion.mfp` | yes | yes | mean-free-path scales (`nucleon_mfp_scale`, `pion_mfp_scale` in a config), not switches |

These map onto GENIE's own `FSI-*-Frac*Scale` / `INUKE-*Fac` configuration
parameters, written into a private override of `HN/HAIntranuke2018.xml` for
the run (recorded in the file's provenance). Three of the hN factors
(nucleon elastic and inelastic, pion inelastic) are a small local GENIE
patch, since stock hN only has them for the pion cex/absorption/elastic
channels. In hA at least one of the four nucleon channels must stay on.

A nucleon left with no allowed fate (e.g. a slow nucleon with elastic
scattering off) leaves the nucleus without interacting; GENIE logs that at
ERROR level for every step, so with a channel off `gcf_fsi` quiets the two
cascade log streams and says so.

### Scans

`gcf scan CONFIG.toml --set OPTION=VALUES` runs a config for several values of an
option: `generate.sigma_cm=0.1,0.15,0.2`, `fsi.scale:pion.absorption=0.5:1.5:5`
(start:stop:count) or `run.seed=1:20:20`. Several `--set` make a grid. Each point
is an ordinary run named `<name>_<option>-<value>`, one after the other. All points
share a base seed, so for an option that only changes the weights, such as a
contact, they hold the same events, and differences between them are smaller than
their statistical errors. With `randomize = true`, a scan of `run.seed` makes
replicas of the nuclear inputs. The web interface's **Scan** tab does the same
through its queue.

### The suite's other generators

`stages = ["program"]` with a `[program]` table runs `genDIS`, `genCCQE`, `genPhoto`,
`genElectro`, `genIon`, `sfSampler` or `sfSamplerLC` the same way as a stage: split
into jobs, merged, resumable (`./gcf example dis`). They write their own `genT` tree,
not the event-file format, so `gcf quicklook`, `gcf_fsi` and the comparison don't
apply. `[program] events` counts attempts. Only `genDIS` takes a seed so far (#137,
PR #150), so the others can't be repeated exactly; the runner passes `-S` as soon as
a generator's usage lists it.

### Without a config

```bash
# 2N events, 8 jobs, with hA FSI into a second file
./gcf generate --model 2N --beam 4.5 -A 12 -Z 6 --events 200000 -j 8 \
               -o c12_2N.root --fsi-output c12_2N_hA.root --fsi-model hA

# FSI on events from another generator that has no target branches
./gcf fsi --input other.root --target 40:20 --fsi-model hN -j 16 -o other_fsi.root
```

The stage programs can also be called directly
(`src/programs/gcf_stages/build/gcf_generate --help`, `gcf_fsi --help`) after
`source build/gcf_env.sh`.

## The web interface

[GUI.md](GUI.md) is the guide to using it. In brief:

```bash
./gcf gui --background       # prints an ssh -L command and a URL with a token
```

The interface runs on the analysis machine. Open an SSH tunnel from your laptop
(`ssh -N -L 8765:localhost:8765 you@machine`; VS Code Remote-SSH forwards the
port automatically) and open the printed URL. From the page you can:

- choose the stages (generate + FSI, generate only, FSI on an existing file);
- fill in the essential settings (model, beam, events, target with presets of
  supported nuclei, wave function, cascade model, interaction channels); the
  rest is under "Advanced", which says how many of its options you changed;
- set the physics model on the "Physics inputs" tab, which shows the table value
  and uncertainty of each nuclear input for the chosen target and interaction;
- see mistakes as you type (the offending field is highlighted), and run
  "Check" to have the generator itself validate the settings (nucleus tables
  and all) before anything starts;
- preview or download the config, and start the run;
- follow the runs: progress with rate and time left, live logs, output files
  with their paths, cancel;
- read the results of a finished run: Σw/attempts with its error, N_eff, the
  heaviest event's share of the weight, warnings about a low N_eff, and the
  provenance of each stage (`gcf summary` prints the same);
- look at quick-look plots, drawn when you open a finished run, and apply a
  selection on Q², x_B, θ_e, p_miss and the proton count (`gcf quicklook --q2 ...`
  does the same);
- compare the outputs of two runs (`tools/validation/compare_samples.py
  --weighted`, see [Validation](#validation));
- reuse a finished run as a template for a new one, or delete it;
- resume a run that stopped, or use the jobs that finished; download outputs and
  full logs; add FSI to a run's PWIA events; load a saved config or an example;
- leave: runs don't depend on the page or a terminal. Runs the page starts share a
  CPU limit (`--max-cpus`) and wait in a queue that survives restarts.

Everything the page does goes through `gcf run`, and each run keeps its
config, so a run set up in the browser can be repeated or scaled up from a
terminal. The server listens on localhost only and every request needs the
token from the URL, because other users of a shared machine can reach localhost
ports. The token is kept in `~/.config/gcf/gui-token`, readable only by you, so the
URL survives restarts; `--new-token` replaces it. `--host` changes the interface it
listens on, over plain HTTP, so anyone on that network can read the token. Stopping
the page does not stop runs already started. When started again it continues runs
that were interrupted, such as by a machine restart, up to 3 times per run.

## Event-file format

One tree, named `events` by default, one entry per event. Units are GeV, GeV/c
and fm; the frame is the target rest frame with the beam along +z.

| Branch | Type | Meaning |
|---|---|---|
| `weight` | double | event weight ([normalization](#weights-and-normalization)) |
| `target_A`, `target_Z` | int | target nucleus |
| `beam_E` | double | beam energy |
| `e_px e_py e_pz e_E` | double | scattered electron |
| `nuc_pdg` | `vector<int>` | particles leaving the primary vertex: the struck nucleon and its correlated partners |
| `nuc_role` | `vector<int>` | 0 = struck nucleon, 1, 2 = partners |
| `nuc_px nuc_py nuc_pz nuc_E` | `vector<double>` | their 4-momenta |
| `vtx_valid` | int | 1 if `vtx_*` hold the interaction vertex |
| `vtx_x vtx_y vtx_z` | double | vertex inside the nucleus (fm) |
| `res_px res_py res_pz res_E` | double | residual nucleus before FSI |

Per file, outside the tree:

- `hAttempts`: a 1-bin TH1D holding the number of generation attempts,
  including rejected ones that are not in the tree. `hadd` sums it.
- `gcf_generate`, `gcf_fsi`, ...: one TNamed per stage with its provenance
  ([Seeds and provenance](#seeds-and-provenance)).

`gcf_generate` writes all of these. The **vertex** matters for FSI: the 2N
generator draws it during generation (it fixes k_F(r) for the local onset and
the local Fermi gas), so stage 1 writes it and stage 2 uses it.

### What `gcf_fsi` needs from a file

Only `nuc_pdg` and `nuc_px nuc_py nuc_pz nuc_E`, plus the target, from
`target_A`/`target_Z` or `--target A:Z`. Everything else is optional:

- no `nuc_role`: the roles are the list positions;
- no vertex (or `vtx_valid == 0`): one is drawn from ρⁿ(r) r² of the target,
  with n = 1, 2, 3 for 1, 2, ≥ 3 nucleons (`--vertex-power` overrides);
- no `weight`, no `hAttempts`: fine, they are only passed through.

Protons and neutrons in the list are transported; anything else (a pion, say)
is passed to the final state unchanged. All input branches are copied to the
output (`--keep REGEX` selects a subset). The input may be a TTree with
`std::vector` or C-array branches of any numeric type, or an RNTuple; for
example, from Python:

```python
import awkward as ak, numpy as np, uproot
with uproot.recreate("my_events.root") as f:
    f.mktree("events", {"nuc_pdg": "var * int32", "nuc_px": "var * float64",
                        "nuc_py": "var * float64", "nuc_pz": "var * float64",
                        "nuc_E": "var * float64", "weight": "float64"})
    f["events"].extend({"nuc_pdg": ak.Array(pdg), "nuc_px": ak.Array(px), ...})
# (f["events"] = {...} writes an RNTuple, which gcf_fsi reads as well)
```

## FSI output

`gcf_fsi` adds, per event:

| Branch | Meaning |
|---|---|
| `fsi_ran` | 1 if the GENIE cascade ran; 0 if it was skipped: residual A < 3 (e.g. 3N on ⁴He) or nothing to transport |
| `fsi_vtx_x y z`, `fsi_vtx_sampled` | vertex used, and whether it was sampled here |
| `fsi_res_A`, `fsi_res_Z` | transport medium: the target minus the transported nucleons |
| `part_n`, `part_pdg`, `part_px py pz E` | the final state: GENIE's stable particles verbatim, or the input unchanged when `fsi_ran == 0` |
| `part_parentRole` | role of the input nucleon a particle descends from |
| `part_rescatterCode` | GENIE fate of that input nucleon: **1 = left without interacting** (the particle then has exactly its input momentum); other values are cascade products; -1 = the cascade did not run |
| `part_generation` | mother-chain depth. Intranuke re-emits even non-interacting nucleons as daughters, so this is ≥ 1 for every cascade output; use `part_rescatterCode`, not this, to find untouched nucleons |

No particle is labelled "lead" or "recoil": that is a detection-level choice
made in the analysis from types and momenta. Pauli blocking is GENIE's own; in
this GENIE version nucleon-nucleon elastic scattering is not Pauli-blocked
(#10).

## Weights and normalization

Events are weighted. `gcf_generate` keeps generating until it has the
requested number of events with weight > 0, and counts every attempt in
`hAttempts`. A weighted yield normalized per attempt, Σw / attempts, is the
quantity to compare between samples. It is printed at the end of each
`gcf_generate` run and in the title of `gcf quicklook`.

Weights in the 3N generator vary over orders of magnitude, so a sample's
effective size N_eff = (Σw)² / Σw² can be far below its event count (also
shown by `quicklook`). FSI does not change weights.

## Seeds and provenance

Each stage takes `--seed`; 0 (the default) draws a seed and prints it. `gcf`
gives job i of a run the seed `seed + i` for generation and `seed + 1000003 + i`
for FSI, so a run with a fixed `[run] seed` is reproducible. Each output file
has a TNamed per stage recording the command line, seed, date and code version:
`git describe` of the tree when the stage was built, with `-dirty` if the code
under `src/` (the vendored GENIE aside) differed from that commit.

## Differences from the combined programs

The split pipeline runs the same physics code (`QEGenerator_3N`,
`QEGeneratorFSI`, the same GENIE calls) but is not bit-identical to the
combined programs: generation and FSI now draw from separate random-number
streams instead of one interleaved stream. The results agree statistically
([Validation](#validation)). Smaller differences:

- **3N vertex**: drawn in stage 1 (ρ³-weighted, as before) instead of right
  before the cascade; the same distribution.
- **2N Coulomb correction** (`--coulomb`): applied at generation, so the
  cascade transports the corrected nucleon. The combined program corrects the
  generated lead nucleon after FSI and leaves the FSI particle list
  uncorrected.
- **Event counts**: `gcf_generate --events N` means N accepted events for
  every model; the combined 2N program counted attempts.
- **Randomized nuclear properties** (`--randomize S`, `genQE_FSI -r`): one draw
  per run, from a random stream of its own, so that all jobs of a run use the
  same inputs. The combined program draws from its event stream.

## Validation

`tools/validation/compare_samples.py A.root B.root` compares two samples: for
each observable an unweighted two-sample χ² test (FSI does not depend on the
weight, so this uses the full statistics), plus the weighted cross section per
attempt. It reads both this format and the combined programs' `genT` trees,
and exits non-zero on disagreement.

An unweighted test can't see a change that only moves weights, such as another
contact, because the events are the same. `--weighted` decides on a weighted
test instead, as the web interface does: the χ² of the normalized weighted
histograms, with a p-value from randomly reassigning events between the samples.
The weights are too heavy-tailed for that statistic to follow a χ² distribution
(on identical settings, 3N gave p < 0.001 in 5-8% of the tests with the
χ² distribution, and in none of 224 with the permutations). It also covers the
observables that respond to the nuclear inputs: p_miss, the partner's momentum,
the angle between them and the residual's momentum. `--max-events` limits the
events used for the shapes; Σw/attempts always uses the whole file.

Split vs combined, same GENIE build, 12C, hN shared transport:

| | 3N, 6 GeV, u = 1 | 2N, 4.5 GeV |
|---|---|---|
| combined | `genQE_3N_FSI`, 400k events | `genQE_FSI`, 640k attempts (405k events) |
| split | `gcf run`, 400k events, 32 jobs | `gcf run`, 430k events, 32 jobs |
| Σw / attempts, pull | +1.71 σ | −0.98 σ |
| smallest p-value of 9 observables (unweighted χ²) | 0.105 | 0.119 |

The observables cover the generated kinematics (Q², x_B, struck-nucleon
momentum) and the post-FSI final state (proton, neutron and charged-pion
multiplicities, fastest-proton momentum and angle, fraction of events without
any FSI interaction). As a negative control, the same 3N events run through hA
instead of hN fail the test decisively (χ²/ndf = 6852/10 for the proton
multiplicity), so the test is sensitive to a real change in the cascade.

## Adding a stage

A new stage (radiative corrections, an acceptance map, ...) should:

1. read the `events` tree, needing only the branches it uses;
2. write all input branches back, plus its own, with an obvious prefix;
3. copy `hAttempts` and the other per-file objects, and add a TNamed with its
   provenance (`src/libraries/io/Provenance.hh`);
4. take `--seed` if it is random, and print `Progress: n/N` lines so `gcf` and
   the GUI can follow it;
5. get a table in `tools/gcf/schema.py`, which makes it available in config
   files, on the command line and in the web form.

`gcf_fsi` (`src/programs/gcf_stages/gcf_fsi.cc`) is the template.
