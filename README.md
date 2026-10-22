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
./gcf gui                         # web interface (open it through an SSH tunnel)
```

The FSI stage also accepts events from other generators, and lets you switch
individual cascade channels (absorption, charge exchange, ...) off. What the
generators compute, assumptions and expressions: [docs/PHYSICS.md](docs/PHYSICS.md). Event format, options
and validation: [docs/PIPELINE.md](docs/PIPELINE.md).

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