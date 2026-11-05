# GCF_Generator_Suite
A suite of event generators utilizing the Generalized Contact Formalism to describe interactions with nuclear SRC pairs.

This repository will contain the code and necessary libraries for generating events involving nuclear SRCs interacting in the PWIA framework.

## Quick start: the pipeline

Generation (PWIA) and final-state interactions (GENIE) can be run as separate
stages, from config files, the command line or a web interface:

```
./install.sh                      # builds GENIE, the generators and the stages
./gcf example 3N > my_run.toml    # then edit it
./gcf run my_run.toml
./gcf gui                         # web interface, see "Using the web interface" below
```

The FSI stage also accepts events from other generators, and lets you switch
individual cascade channels (absorption, charge exchange, ...) off. What the
generators compute, assumptions and expressions: [docs/PHYSICS.md](docs/PHYSICS.md). Event format, options
and validation: [docs/PIPELINE.md](docs/PIPELINE.md).

## Using the web interface

`./gcf gui` starts a web page for setting up quasi-elastic runs, starting them
and looking at the results. It runs everything through `./gcf run` and saves each
run's config, so anything you start in the browser can be repeated from a
terminal.

The page covers the `gcf_generate` and `gcf_fsi` stages only. Run the other
generators (`genDIS`, `genPhoto`, ...) and the EIC chain
(`tools/run_eic_sim.sh`) from the terminal.

### Start it

You need the programs built with `./install.sh` and Python 3.11 or newer.

The plots and summaries also need `uproot`, `awkward`, `numpy` and `matplotlib`,
and comparing runs needs `scipy`. `install.sh` does not install them. The rest of
the page works without them. A virtual environment is the simplest way to add
them; activate it in the shell where you start the page:

```
python3 -m venv ~/gcf-venv        # on Ubuntu this needs the python3-venv package
source ~/gcf-venv/bin/activate
pip install uproot awkward numpy matplotlib scipy
```

If you use `uv`, `uv venv ~/gcf-venv` and `uv pip install` do the same job.

Then start the page on the machine that will run the jobs:

```
./gcf gui                                # runs are written to <repo>/runs
./gcf gui --runs-dir /path/with/space    # choose where runs are written
./gcf gui --port 8800                    # if the default port 8765 is taken
```

It prints an SSH tunnel command and a URL:

1. On your laptop, run the `ssh -N -L ...` command it printed and leave it
   running. VS Code Remote-SSH forwards the port for you, and you don't need a
   tunnel if the browser runs on the same machine.
2. Open the URL in your browser, including the `?token=...` part. The token is
   the password, and it changes every time the page is started.

Ctrl-C or closing the terminal stops the page but not the runs. Start it again
with the same `--runs-dir` to see them, using the new URL. To keep the page up
after you log out, start it inside `tmux`.

### Set up a run

1. Choose what to run at the top of the form: **Generate + FSI**, **Generate
   only**, or **FSI on a file**.
2. Fill in the form. The fields change with the generator model (3N, 2N or MF),
   and the **?** next to a field explains it. **Events** is the number of
   accepted events over all jobs. Less common settings, such as the eN cross
   section, `sigma_CM`, the Q² window and E\*, are under **Advanced**.
3. For FSI, choose the cascade model (hN or hA). Untick an interaction channel
   to switch it off, or type a factor next to it to scale its rate.
4. Fix anything marked in red; the form checks your values as you type.
   **Check & show config** also asks the generator to validate the settings and
   shows the resulting config. **Download config** saves it as a TOML file.
5. Press **Start run**. The run appears under **Runs**.

The **Physics inputs** tab holds the settings of the physics model: the 2N
contacts, `sigma_CM`, E\* and its width, the cross section and form factors,
radiation, the cutoffs. Each field shows the value from the generator's tables
for your target and interaction, with its uncertainty, and an empty field uses
that value. None of this needs a rebuild. **Vary the inputs within their errors**
makes one random draw per run, for systematics studies. Contacts only change event
weights, so for a contact study you can also reweight existing events (`band` in
`tools/eic_src_analysis.py` does this).

For **FSI on a file**, use **Browse…** to pick an event file that has not been
through FSI yet: a `*_pwia.root` file, or another generator's file in the
[event-file format](docs/PIPELINE.md#event-file-format). A `*_fsi.root` file
fails as soon as the run starts.

### Follow a run and use the results

Click a run to open its details:

- **Progress and logs:** while it runs you see the rate, the time left and the
  last 40 lines of each job's log. **Cancel run** stops it.
- **Output files:** these are `<name>_pwia.root` and `<name>_fsi.root` in
  `<runs-dir>/<name>/`, next to `config.toml` and the full logs in `logs/`. The
  page can't download event files, so use **Copy path** and read or `scp` them
  from the machine.
- **Numbers:** each output lists its events, attempts, Σw / attempts (the
  [quantity to compare between samples](docs/PIPELINE.md#weights-and-normalization)),
  the effective number of events N_eff and what made it. Read the warnings before
  trusting a shape: 4,000 3N events gave N_eff = 20 in one test. `./gcf summary
  FILE` prints the same numbers in a terminal.
- **Plots:** control plots are drawn when you open a finished run. Under
  **Selection**, set limits on Q², x_B, θ_e, p_miss or the proton count to redraw
  them and see what share of the yield passes.
- **Compare…:** compares the outputs of two runs: Σw / attempts with the pull
  between them, and whether the shapes of the generated and final-state
  distributions differ. The test uses the event weights, so it also sees changes
  that only move weights. "No difference" means none larger than the statistics
  can show, so check the N_eff it quotes. In a terminal:
  `python3 tools/validation/compare_samples.py --weighted A.root B.root`.
- **Use as template:** fills the form from that run, with a new seed.
- **Copy run command:** copies the command that repeats the run from the
  repository root. Add `--force` to overwrite the run, or change `name` in a
  copy of its `config.toml` to keep both. With the saved seed and the same
  number of jobs, you get the same events.
- **Delete:** removes the run's whole directory, including its event files.

The **docs** link at the top of the page opens the pipeline and physics
documentation.

### Things to watch for

- **A passed check doesn't guarantee a good run.** It checks the options and
  that the generator can produce events, not everything about the input. For
  example, FSI on a `*_fsi.root` file passes and then fails when the run starts.
- **Switching off the compound-nucleus channel does nothing in hA** (issue #177).
- **Cancelled runs can't be resumed.** Finished job files stay in the run's
  `chunks/` directory, unmerged. Start a new run instead.
- **Runs don't share a CPU limit.** **Parallel jobs** is only checked against
  the machine's CPU count for each run, so check what is already running before
  you start another large one.

More detail: [docs/PIPELINE.md](docs/PIPELINE.md#the-web-interface).

## Compile

Navigate to a build directory of your choice and execute the following:

```
ccmake [path/to]/GCF_Generator_Suite/src
```
Press [c] and then [g]
```
make
```
The executables will then be located in subdirectories of the ```[build]/programs``` directory. For instance, the Quasielastic Scattering executable is located at:
```
[build]/programs/genQE/genQE
```

## Authors

In order of first contribution (from the git history):

- Jackson Pybus ([@JacksonPybus](https://github.com/JacksonPybus)), 2020–2023
- Phoebe Sharp ([@psharp15](https://github.com/psharp15)), 2020–2021
- Alex Kiral ([@a-kiral](https://github.com/a-kiral)), 2020
- Axel Schmidt ([@schmidta87](https://github.com/schmidta87)), 2020–2023
- Andrew Denniston ([@awild7](https://github.com/awild7)), 2025 (3N generator)
- Idey ([@ilayws](https://github.com/ilayws)), 2026
- Tom Bleher ([@tom-bleher](https://github.com/tom-bleher)), 2026

`git shortlog -sn` gives per-author commit counts; `.mailmap` merges each author's name and email variants.

## Meta

This repository is maintained by Jackson Pybus (jrpybus@mit.edu).

[https://github.com/JacksonPybus](https://github.com/JacksonPybus)