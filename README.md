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
./gcf gui --background            # web interface, see "Using the web interface" below
```

The FSI stage also accepts events from other generators, and lets you switch
individual cascade channels (absorption, charge exchange, ...) off. The `eic`
stage (`./gcf example eic`, or "EIC / ePIC" in the web interface) runs the events
through the ePIC simulation and reconstruction. What the
generators compute, assumptions and expressions: [docs/PHYSICS.md](docs/PHYSICS.md). Event format, options
and validation: [docs/PIPELINE.md](docs/PIPELINE.md).

## Using the web interface

`./gcf gui` is a web page for setting up event generation, running it and reading the
results. **[docs/GUI.md](docs/GUI.md) is the full guide**, and the page's **docs** link
opens it too. In short:

```
./install.sh                 # once: builds the programs and sets up the Python packages
./gcf doctor                 # checks that everything is ready, and says how to fix what isn't
./gcf gui --background       # starts the page; prints an SSH command and an address
```

On your laptop, run the `ssh -N -L ...` command it printed, then open the address in a
browser, including the `?token=...` part. The token is the password, so don't share the
address. `./gcf gui --status` shows it again, and `./gcf gui --stop` stops the page.

Then:

1. On a new page, press **Fill in a small test run**, then **Start run**. It takes under
   a minute.
2. Click the run to see its numbers, warnings, plots and output files. **Σw / attempts**
   is the quantity to compare between samples, and **N_eff** tells you how far to trust
   them.
3. Set up your own run in the form: **Run setup**, **Physics inputs** (contacts,
   `sigma_CM`, E\* and more, with their table values) and **Scan** (several values of a
   parameter). **Check & show config** tests the settings and estimates the time and
   disk space.

**You can leave.** Runs don't need `tmux` or `screen`. They are started on their own, so
you can close the page or log out, and they carry on. Come back to the same address, or
use `./gcf status`. A run that was interrupted can always continue from what it finished,
with the same events. See [Leaving and coming back](docs/GUI.md#leaving-and-coming-back).

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