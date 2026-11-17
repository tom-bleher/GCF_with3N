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
opens it too.

Once, on the analysis machine:

```
./install.sh                 # builds the programs and sets up the Python packages
./gcf doctor                 # checks that everything is ready, and says how to fix what isn't
```

### From a Mac over SSH (Terminal or iTerm2)

1. Logged in to the analysis machine, in this checkout, run `./gcf gui`. It starts the
   page in the background (it keeps running after you log out) and prints one line,
   with your user name, the machine's address and this checkout filled in, like:

   ```
   ssh you@machine cat src/GCF_with3N/tools/gcf/gcf-open | sh -s you@machine
   ```

2. Press **Cmd-T** for a new tab on your Mac, paste that line there and press Enter.
   It starts the page if it isn't running, opens a tunnel to it and opens your browser
   on it.
3. It asks `Add the shortcut gcf-gui to ~/.zshrc ...? [Y/n]`. Press Enter.
4. From then on, type **`gcf-gui`** in any terminal on your Mac and the page opens.
   No need to log in first. After the Mac has slept or changed network, just run it
   again: it reconnects, at the same address, so a bookmark keeps working.

With SSH keys you're never asked for a password; without them, once per connection.
The same works from a Linux laptop, WSL or Git Bash. A link printed on the remote
machine can't open the page by itself: only your laptop's `ssh` can open the tunnel,
which is what that line does.

### Elsewhere

In a VS Code or Cursor remote terminal, or on a desktop you're sitting at, `./gcf gui`
opens the page by itself.

The token in the address (`?token=...`) is the password, so don't share it.
`./gcf gui --status` shows how to open the page again, and `./gcf gui --stop` stops it
(runs already started keep going).

### Your first run

Once the page is open:

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