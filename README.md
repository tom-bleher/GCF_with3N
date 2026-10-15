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
individual cascade channels (absorption, charge exchange, ...) off. Event format, options
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

## Meta

This repository is maintained by Jackson Pybus (jrpybus@mit.edu).

[https://github.com/JacksonPybus](https://github.com/JacksonPybus)