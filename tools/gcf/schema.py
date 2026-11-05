"""Options of the pipeline stages, in one place.

Every option the front end knows about is described here once: its type,
default, the gcf_generate / gcf_fsi flag it maps to, which generator models it
applies to, and a help text. The TOML config validation, the command lines of
the stage binaries and the GUI form are all built from these tables.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional

MODELS = ["3N", "2N", "MF"]

# 3N wave functions (Faddeev 3He densities): name -> table number of gcf_generate.
WAVE_FUNCTIONS_3N = ["AV18", "AV18+UIX", "AV18+TM", "CDBonn", "CDBonn+TM",
                     "AV18-pD", "AV18+UIX-pD", "AV18+TM-pD", "CDBonn-pD", "CDBonn+TM-pD"]
NN_INTERACTIONS = ["AV18", "N2LO", "N2LO12", "N3LO", "AV4", "NV"]

# Targets whose masses (and, for 3N, the residual's) are in the generators' tables.
TARGETS = {
    "3N": [(3, 2, "3He"), (4, 2, "4He"), (6, 3, "6Li"), (12, 6, "12C"), (14, 7, "14N")],
    "2N": [(2, 1, "2H"), (3, 1, "3H"), (3, 2, "3He"), (4, 2, "4He"), (6, 3, "6Li"), (12, 6, "12C"),
           (27, 13, "27Al"), (40, 18, "40Ar"), (40, 20, "40Ca"), (56, 26, "56Fe"), (208, 82, "208Pb")],
}
# MF needs the A-1 residual masses, which the nucleus tables have only for 12C.
TARGETS["MF"] = [(12, 6, "12C")]

# Interaction channels of the GENIE cascade that gcf_fsi --channel can switch
# off, and whether each is a separate fate in the hN / hA model.
CHANNELS = [
    ("nucleon.elastic",    "elastic scattering",   True,  False),
    ("nucleon.inelastic",  "inelastic scattering", True,  True),
    ("nucleon.cex",        "charge exchange",      False, True),
    ("nucleon.absorption", "absorption",           False, True),
    ("nucleon.piprod",     "pion production",      False, True),
    ("nucleon.compound",   "compound nucleus",     True,  False),   # hA reads the flag but never uses it (#177)
    ("pion.elastic",       "elastic scattering",   True,  False),
    ("pion.inelastic",     "inelastic scattering", True,  True),
    ("pion.cex",           "charge exchange",      True,  True),
    ("pion.absorption",    "absorption",           True,  True),
    ("pion.piprod",        "pion production",      False, True),
]
# Channels gcf_fsi can only switch, not scale (a yes/no parameter of GENIE's).
SWITCH_CHANNELS = {"nucleon.compound"}


@dataclass
class Opt:
    key: str                     # name in the TOML config / GUI
    kind: str                    # float, int, str, bool, range, target
    flag: Optional[str]          # binary flag (None: handled by the runner)
    default: Any = None          # None = not passed (binary default applies)
    help: str = ""
    label: str = ""
    choices: Optional[List[Any]] = None
    # Models the option applies to; gcf_generate rejects it for the others with
    # the same table (kModelOptions in gcf_generate.cc), and
    # src/programs/gcf_stages/test_model_options.py checks the two agree.
    models: List[str] = field(default_factory=lambda: list(MODELS))
    invert: bool = False         # bool flag that is passed when the value is False
    unit: str = ""
    advanced: bool = False       # GUI: shown under "Advanced", folded by default
    tab: str = ""                # GUI: shown on this tab of the form ("" = the main one)
    ge: Optional[float] = None   # allowed values (for a range, of both ends): >= this,
    gt: Optional[float] = None   # > this,
    le: Optional[float] = None   # and <= this
    choice_labels: Optional[dict] = None   # GUI: text of each choice, when not the value itself


GENERATE = [
    Opt("model", "str", "--model", "3N", choices=MODELS, label="Generator model",
        help="3N: struck nucleon from an SRC triplet. 2N: from an SRC pair. "
             "MF: single mean-field nucleon from a Fermi gas."),
    Opt("beam", "float", "--beam", 6.0, label="Beam energy", unit="GeV", gt=0, le=1e5),
    Opt("events", "int", "--events", 10000, label="Events", ge=1,
        help="Number of accepted (weight > 0) events in total, over all jobs."),
    Opt("target", "target", None, [12, 6], label="Target (A, Z)"),
    # Kelly: PRC 70, 068202 (2004), doi:10.1103/PhysRevC.70.068202
    # Ye et al.: PLB 777, 8 (2018), doi:10.1016/j.physletb.2017.11.023
    Opt("ff", "str", "--ff", "kelly", choices=["kelly", "ye", "dipole"],
        label="Nucleon form factors",
        help="Kelly 2004 (default), Ye et al. 2018 z expansion, or dipole with Galster G_En.",
        tab="physics"),
    # De Forest: Nucl. Phys. A 392, 232 (1983), doi:10.1016/0375-9474(83)90124-0
    Opt("cs", "str", "--cs", "cc1", choices=["cc1", "cc2", "onshell"],
        label="eN cross section", help="De Forest cc1 (default), cc2, or on-shell.", tab="physics"),
    Opt("sigma_cm", "float", "--sigma-cm", None, models=["3N", "2N"], label="sigma_CM", unit="GeV/c", ge=0, le=1,
        help="CM momentum width of the SRC cluster per component. Default: 0.15 (3N), "
             "the nucleus table value (2N).", tab="physics"),
    Opt("q2", "range", "--q2", None, label="Q^2 window", unit="GeV^2", ge=0,
        help="Generation-time Q^2 window MIN:MAX. Default: open (3N), 1:10 (2N/MF).", advanced=True),
    # 3N
    Opt("wave_function", "str", "--wave-function", "AV18", choices=WAVE_FUNCTIONS_3N,
        models=["3N"], label="3N wave function",
        help="Nucleon interaction behind the three-body density (AV18, AV18+UIX, AV18+TM, "
             "CDBonn, CDBonn+TM). The 3He Faddeev density is used for the SRC triplet in "
             "every nucleus. The -pD variants keep only the deuteron (pD) channel."),
    Opt("theta_e", "range", "--theta-e", [5.0, 40.0], models=["3N"], ge=0, le=180,
        label="Electron angle window", unit="deg", advanced=True),
    Opt("kF", "float", "--kF", None, models=["3N"], label="k_F cutoff", unit="GeV/c", gt=0, le=0.5,
        help="Internal momenta k > k_F (default 0.25); also sets sigma_CM = "
             "sqrt(3/5) k_F unless sigma_CM is given.", tab="physics"),
    Opt("cm", "bool", "--no-cm", True, models=["3N"], invert=True,
        label="CM motion of the triplet", tab="physics"),
    # 2N
    Opt("nn", "str", "--nn", None, models=["2N"], label="NN interaction",
        choices=NN_INTERACTIONS,
        help="NN interaction behind the pair wave function and contacts (default AV18)."),
    Opt("estar", "float", "--estar", None, models=["3N", "2N"], label="E*", unit="GeV", ge=0, le=1,
        help="Excitation energy of the residual nucleus; default 0. Not for a target with no residual.", tab="physics"),
    Opt("sigma_e", "float", "--sigma-e", None, models=["2N"], label="E* width", unit="GeV", gt=0, le=1,
        help="Draw E* per event from a Gaussian of this width around E*, keeping E* > 0. "
             "Empty: the same E* for every event.", tab="physics"),
    Opt("cpp0", "float", "--cpp0", None, models=["2N"], label="C_pp (S=0)", ge=0, tab="physics",
        help="pp contact, in the units of the nucleus table. A contact only changes the event weights."),
    Opt("cnn0", "float", "--cnn0", None, models=["2N"], label="C_nn (S=0)", ge=0, tab="physics",
        help="nn contact, in the units of the nucleus table. A contact only changes the event weights."),
    Opt("cpn0", "float", "--cpn0", None, models=["2N"], label="C_pn (S=0)", ge=0, tab="physics",
        help="Spin-singlet pn contact, in the units of the nucleus table. A contact only changes the event weights."),
    Opt("cpn1", "float", "--cpn1", None, models=["2N"], label="C_pn (S=1)", ge=0, tab="physics",
        help="Spin-triplet (deuteron-channel) pn contact, in the units of the nucleus table. "
             "A contact only changes the event weights."),
    Opt("randomize", "bool", "--randomize", False, models=["2N"], label="Vary the inputs within their errors",
        help="One draw per run: sigma_CM and the contacts from Gaussians of their uncertainties, E* uniform up to the "
             "table's maximum, using the run's seed. Values you set here override the draw. For systematics studies.",
        tab="physics"),
    Opt("light_cone", "bool", "--light-cone", False, models=["2N"], label="Light-cone cross section",
        help="Instead of the instant-form one. Excludes importance sampling of |p_rel| and the local onset.", tab="physics"),
    Opt("kcut", "float", "--kcut", None, models=["2N"], label="k_rel cutoff", unit="GeV/c", gt=0, le=1.05,
        help="Relative-momentum floor of the SRC pair (default 0.25).", tab="physics"),
    Opt("onset", "str", "--onset", "fixed", choices=["fixed", "local"], models=["2N"],
        label="SRC onset", help="fixed: k_rel > cutoff. local: k_rel > k_F(r) at the vertex.", tab="physics"),
    Opt("prel_is", "bool", "--prel-is", False, models=["2N"],
        label="Importance-sample |p_rel|", help="Lower weight variance; same expectation values.", advanced=True),
    # MF
    Opt("fermi_gas", "str", "--fermi-gas", "global", choices=["global", "local"],
        models=["MF"], label="Fermi gas",
        help="global: one Fermi momentum for the whole nucleus. local: k_F(r) from the density."),
    Opt("kF_fg", "float", "--kF-fg", None, models=["MF"], label="Global Fermi momentum", gt=0, le=1,
        unit="GeV/c", help="Default 0.25.", tab="physics"),
    # 2N and MF
    Opt("rad", "bool", "--rad", False, models=["2N", "MF"], label="Radiation (peaking approx.)", tab="physics"),
    Opt("coulomb", "bool", "--coulomb", False, models=["2N", "MF"], label="Coulomb correction",
        help="Applied at generation, i.e. before FSI.", tab="physics"),
    Opt("phase_space", "str", "--phase-space", None, models=["2N", "MF"],
        label="Phase-space file", help="gcfGenerator phase-space limits file.", advanced=True),
]

FSI = [
    Opt("input", "str", "--input", None, label="Input event file",
        help="Only for FSI-only runs: any file in the GCF event-file format."),
    Opt("model", "str", "--fsi-model", "hN", choices=["hN", "hA"], label="Cascade model",
        help="hN: full hadron-nucleon cascade. hA: effective hadron-nucleus model (faster)."),
    Opt("off", "list", None, [], choices=[c[0] for c in CHANNELS],
        label="Interaction channels switched off",
        help="Fates of the cascade to switch off (their rate is scaled to 0, the rest "
             "is renormalized). Which ones exist depends on the model: see --list-channels."),
    Opt("scale", "scales", None, {}, label="Interaction channel scale factors",
        help="Scale the rate of a fate of the cascade by a factor; the rest is renormalized. "
             "Which ones exist depends on the model: see --list-channels."),
    Opt("nucleon_mfp_scale", "float", None, None, label="Nucleon mean-free-path scale", gt=0,
        help="Multiplies the nucleon mean free path (> 1: fewer interactions).", advanced=True),
    Opt("pion_mfp_scale", "float", None, None, label="Pion mean-free-path scale", gt=0,
        help="Multiplies the pion mean free path.", advanced=True),
    Opt("transport", "str", None, "shared", choices=["shared", "independent"],
        label="Transport",
        help="shared: all nucleons in one GENIE record (absorption depletes the "
             "remnant for the others). independent: one record each.", advanced=True),
    Opt("vertex_power", "str", "--vertex-power", "auto", choices=["auto", "1", "2", "3"],
        label="Vertex weighting (if the file has none)",
        help="rho^n(r) r^2; auto picks n = number of nucleons (max 3).", advanced=True),
    Opt("target", "target", None, None, label="Target (A, Z) for files without one", advanced=True),
    Opt("keep", "str", "--keep", None, label="Input branches to keep (regex)",
        help="Default: all.", advanced=True),
    Opt("tree", "str", "--tree", None, label="Input tree", help="Default: events.", advanced=True),
    Opt("overwrite", "bool", "--overwrite", False, label="Replace earlier FSI in the input",
        help="Without this, a file that has already been through FSI is refused.", advanced=True),
]

RUN = [
    Opt("name", "str", None, "run", label="Run name",
        help="Output files are <output_dir>/<name>/<name>_pwia.root and _fsi.root."),
    Opt("output_dir", "str", None, "runs", label="Output directory",
        help="Relative paths are relative to the config file's directory."),
    Opt("stages", "list", None, ["generate", "fsi"], label="Stages",
        choices=["generate", "fsi", "program", "eic"]),
    Opt("jobs", "int", None, 1, label="Parallel jobs",
        help="Independent jobs run at the same time; their outputs are merged. "
             "Each needs one CPU core."),
    Opt("seed", "int", None, 0, label="Base seed",
        help="0 = draw one and record it in the saved config. Each job's seed is derived from it (ePIC: slice i uses seed + i).", advanced=True),
    Opt("note", "str", None, None, label="Note", advanced=True,
        help="What this run is for. Shown in the list of runs and searched by its filter."),
    Opt("chunks", "int", None, None, label="Split into jobs", advanced=True,
        help="More jobs than run at a time means an interruption loses less: a resume reruns only the jobs "
             "that didn't finish. Default: one per million events, and at least Parallel jobs. "
             "Each job has its own seed, so this is part of what makes a run reproducible."),
    Opt("keep_chunks", "bool", None, False, label="Keep per-job files", advanced=True),
    Opt("scan", "str", None, None, label="Scan", help="The parameter scan this run is a point of."),
]

# ---------------------------------------------------------------------------
# The suite's other generators. They write their own tree (genT) rather than the
# event-file format, and most have no seed option yet (#137, PR #150).

@dataclass
class Program:
    label: str
    beam: str = ""               # what the beam-energy argument means ("" = no such argument)
    pdf: bool = False            # takes a PDF table
    note: str = ""               # what a user must know before trusting the output


PROGRAMS = {
    "genDIS": Program("DIS off an SRC pair", beam="Electron beam energy", pdf=True,
                      note="Known issues: W² comes out too high (#63, PR #88), and radiation has no effect (#75, PR #95)."),
    "genCCQE": Program("Neutrino CCQE off an SRC pair",
                       note="Known issue: its cross section is a rescaled electromagnetic one, not a true CCQE cross section (#66)."),
    "genPhoto": Program("Photoproduction off an SRC pair",
                        note="Known issue: the default pi- p channel can put most of the summed weight in one event (#132). "
                             "Check N_eff."),
    "genElectro": Program("Electroproduction of J/psi off an SRC pair", beam="Electron beam energy",
                          note="Known issues: J/psi events get the wrong outgoing baryon (#67, PR #96), and the form factor "
                               "squares a mass that is already squared (#76, PR #100)."),
    "genIon": Program("Proton-nucleus scattering off an SRC pair", beam="Beam energy per nucleon",
                      note="Known issues: the pp cross section is zero above s = 8.3 GeV², which excludes the BM@N energy "
                           "(#138, PR #174), and it is evaluated at a t that can be unphysical (#175)."),
    "sfSampler": Program("Spectral function samples"),
    "sfSamplerLC": Program("Light-cone spectral function samples"),
}
_ALL_PROGRAMS = list(PROGRAMS)
_PAIR_INPUTS = ["genPhoto", "sfSampler", "sfSamplerLC"]
_SAMPLERS = ["sfSampler", "sfSamplerLC"]

PROGRAM = [
    Opt("program", "str", None, "genDIS", choices=_ALL_PROGRAMS, models=_ALL_PROGRAMS, label="Generator",
        help="The suite's other generators. They write their own tree, genT, so plots and FSI are not available for them."),
    Opt("events", "int", None, 100000, models=_ALL_PROGRAMS, label="Attempts", ge=1,
        help="Number of attempts in total, over all jobs. Only those with a positive weight are stored."),
    Opt("target", "target", None, [12, 6], models=_ALL_PROGRAMS, label="Target (A, Z)"),
    Opt("beam", "float", None, 6.0, models=[k for k, v in PROGRAMS.items() if v.beam], label="Beam energy", unit="GeV", gt=0, le=1e5),
    Opt("pdf", "str", None, None, models=["genDIS"], label="PDF table",
        help="A CTEQ .pds file. Default: CT18NN, from the repository."),
    Opt("reaction", "str", "-R", None, models=["genPhoto"], label="Reaction",
        choices=["pim", "rho0", "rhom", "omega", "phi", "phin", "deltapp", "deltappim"], help="Default: pi- p."),
    Opt("fixed_beam", "float", "-B", None, models=["genPhoto"], label="Fixed beam energy", unit="GeV", gt=0,
        help="Default: the GlueX spectrum."),
    Opt("amorphous", "bool", "-A", False, models=["genPhoto"], label="Amorphous radiator spectrum"),
    Opt("light_cone", "bool", "-l", False, models=["genPhoto"], label="Light-cone cross section"),
    Opt("inverse", "bool", "-I", False, models=["genIon"], label="Inverse kinematics",
        help="Boost to inverse kinematics; the beam energy is then the kinetic energy per nucleon."),
    Opt("nn", "str", "-u", None, models=_PAIR_INPUTS, label="NN interaction", choices=NN_INTERACTIONS),
    Opt("kcut", "float", "-k", None, models=_PAIR_INPUTS, label="k_rel cutoff", unit="GeV/c", gt=0, le=1.05),
    Opt("sigma_cm", "float", "-s", None, models=_SAMPLERS, label="sigma_CM", unit="GeV/c", ge=0, le=1),
    Opt("estar", "float", "-E", None, models=["genDIS"] + _SAMPLERS, label="E*", unit="GeV", ge=0, le=1),
    Opt("estar_12c", "bool", "-M", False, models=["genDIS", "genCCQE"], label="12C E* distribution",
        help="E* = 17.32 MeV with a 9.57 MeV spread, measured on 12C; it is applied whatever the target (#153)."),
    Opt("randomize", "bool", "-r", False, models=_SAMPLERS, label="Vary the nuclear inputs",
        help="Varied separately in each job, so use one job if they must all have the same inputs."),
    Opt("rad", "bool", "-O", False, models=["genDIS"], label="Radiation (peaking approx.)",
        help="Has no effect until PR #95 is merged (#75)."),
    Opt("phase_space", "str", "-P", None, models=_ALL_PROGRAMS, label="Phase-space file",
        help="gcfGenerator phase-space limits file."),
]

# ---------------------------------------------------------------------------
# The ePIC chain: tools/run_eic_sim.sh, one run of it per beam setting and Q^2 slice.

# The light-ion settings of run_eic_sim.sh: geometry config and the beamline field swapped into it.
# Defaults of that script, repeated here only to check the install before a run starts.
EIC_BEAMS = {"5x41": ("epic_craterlake_5x41", "fields/beamline_5x41_He4.xml"),
             "10x110": ("epic_craterlake_10x110_H2", None)}
# run_eic_sim.sh's targets: Z/A = 1/2, the rigidity of those lattices.
EIC_IONS = {"C12": "12C", "d": "2H", "He4": "4He"}
# seed * 100000 + chunk seeds npsim; this keeps it within 32 bits
EIC_MAX_SEED = 21000


@dataclass
class EicSample:
    label: str
    generator: str               # run_eic_sim.sh --generator
    gen_args: str = ""           # its --gen-args
    cuts: bool = True            # Q2 > 3, xB > 1.2 of PRC 105, 034001 (False: --no-cuts)
    deexcite: tuple = ("none", "geant4", "beagle")   # allowed --deexcite; the first but "none" is "auto"
    pass_fraction: float = 1.0   # events written per --attempts in the production runs
    production: dict = field(default_factory=dict)   # ion -> slices [q2 min, q2 max, attempts]
    ions: tuple = ("C12", "He4")  # the deuteron has no A-2 system and no cascade: SRC pairs only


_SRC_D = [[3, 5, 811000], [5, 12, 338000], [12, 25, 205000]]
_DE_EXCITE_RESIDUAL = ("none", "geant4")   # genQE_FSI: BeAGLE's GCF interface needs the 2N event layout
EIC_SAMPLES = {
    "src": EicSample("SRC pairs, no FSI (genQE)", "genQE", pass_fraction=0.125,
                     ions=("C12", "d", "He4"),
                     production={"C12": [[3, 5, 960000], [5, 12, 377000], [12, 25, 226000]], "d": _SRC_D}),
    "fsi": EicSample("SRC pairs + GENIE hN FSI (genQE_FSI)", "genQE_FSI", deexcite=_DE_EXCITE_RESIDUAL,
                     pass_fraction=0.125, production={"C12": [[3, 5, 496000], [5, 12, 192000], [12, 25, 115000]]}),
    "mf": EicSample("Mean field, no FSI (genQE_FSI -n -F 0.25)", "genQE_FSI", "-n -F 0.25", cuts=False,
                    deexcite=_DE_EXCITE_RESIDUAL,
                    production={"C12": [[3, 5, 25000], [5, 12, 15000], [12, 25, 10000]]}),
    "3n": EicSample("3N SRC, no FSI (genQE_3N_FSI -n)", "genQE_3N_FSI", "-n", deexcite=("none",), pass_fraction=0.08,
                    production={"C12": [[3, 5, 375000], [5, 12, 150000], [12, 25, 100000]]}),
    "3n_fsi": EicSample("3N SRC + GENIE hN FSI (genQE_3N_FSI)", "genQE_3N_FSI", deexcite=("none",), pass_fraction=0.08,
                        production={"C12": [[3, 5, 375000], [5, 12, 150000], [12, 25, 100000]]}),
    "dis": EicSample("Tagged DIS off an SRC pair (genDIS)", "genDIS", cuts=False, pass_fraction=0.4,
                     production={"C12": [[3, 6, 56000], [6, 12, 22500], [12, 25, 15000]]}),
}


def eic_production(sample, ion):
    """Production slices of a sample; targets it was not run on borrow the 12C statistics."""
    table = EIC_SAMPLES[sample].production
    return table.get(ion, table["C12"])


def eic_test_slices(sample, ion):
    """1% of the production statistics: enough to see every step work in minutes."""
    return [[a, b, max(100, n // 100)] for a, b, n in eic_production(sample, ion)]


# Rough cost per event through npsim + eicrecon, from the 2026-10 productions on 12C at 5x41 and 10x110.
EIC_MB_PER_EVENT = 1.05
EIC_CORE_S_PER_EVENT = {"5x41": 3.0, "10x110": 4.5}
EIC_GB_PER_SIM_JOB = {True: 2.5, False: 2.9}     # single-threaded npsim job: True; multithreaded: False
EIC_GB_PER_RECO_JOB = 3.1

EIC = [
    Opt("sample", "str", None, "src", choices=list(EIC_SAMPLES), label="Sample",
        choice_labels={k: v.label for k, v in EIC_SAMPLES.items()},
        help="What goes through the detector. SRC samples keep the Q^2 > 3, x_B > 1.2 cuts of PRC 105, 034001; "
             "mean field and DIS convert every event."),
    Opt("beams", "list", None, ["5x41"], choices=list(EIC_BEAMS), label="Beam settings",
        help="Electron GeV x ion GeV/c per nucleon. Each setting runs after the other; the figures overlay them."),
    Opt("ion", "str", None, "C12", choices=list(EIC_IONS), label="Target", choice_labels=EIC_IONS,
        help="Z/A = 1/2 only: the rigidity of the light-ion lattices in epic."),
    Opt("slices", "slices", None, eic_test_slices("src", "C12"), label="Q^2 slices",
        help="Each slice is its own run with its own seed and attempts, normalized by its own attempts, so disjoint "
             "slices add. Attempts are genQE attempts (for the 3N samples, accepted generator events)."),
    Opt("deexcite", "str", None, "auto", choices=["auto", "none", "geant4", "beagle"], label="Remnant de-excitation",
        choice_labels={"auto": "auto", "none": "none", "geant4": "Geant4", "beagle": "BeAGLE + FLUKA"},
        help="Evaporation nucleons, fragments and photons from the A-2 remnant, before the afterburner: Geant4 "
             "(tools/deexcite) or BeAGLE with FLUKA, as PRC 105, 034001. auto: Geant4 if the sample has a remnant. "
             "genQE_FSI samples: Geant4 only. 3N samples have no remnant."),
    Opt("sim_jobs", "int", None, 8, label="npsim jobs",
        help="Parallel npsim jobs; about 2.5 GB each single-threaded, 2.9 GB with threads."),
    Opt("sim_threads", "int", None, 1, label="Threads per npsim job",
        help="32 jobs x 8 threads ran 1.3x the events/s of 70 single-threaded jobs at half the memory; "
             "one job with very many threads is far less efficient."),
    Opt("reco_jobs", "int", None, 8, label="eicrecon jobs",
        help="Parallel eicrecon jobs, about 3.1 GB each. Always single-threaded: multithreaded eicrecon "
             "segfaults in the dRICH ring finder."),
    Opt("analysis", "bool", None, True, label="Figures when done",
        help="tools/eic_src_analysis.py extract and plots over all slices: PRC 105, 034001 figure set and summary.json."),
    # Software; the placeholder of each field says what an empty one resolves to, and the saved config records the path.
    Opt("epic_install", "str", None, None, label="epic install", advanced=True,
        help="Install prefix of the geometry. Empty: $EPIC_INSTALL, else the container's nightly. Reproducing the "
             "productions needs the pinned tree with the realistic B0 (see docs/PIPELINE.md)."),
    Opt("eicrecon_install", "str", None, None, label="EICrecon install", advanced=True,
        help="Install prefix of EICrecon. Empty: $EICRECON_INSTALL, else the one in the container."),
    Opt("eic_shell", "str", None, None, label="eic-shell", advanced=True,
        help="The eic-shell launcher. Empty: $EIC_SHELL, else ~/eic/eic-shell, else eic-shell on PATH."),
    Opt("fieldmaps", "str", None, None, label="Field maps directory", advanced=True,
        help="Linked as ./fieldmaps in every job. Empty: $FIELDMAPS, else ~/eic/fieldmaps."),
    Opt("chunk_size", "int", None, 500, label="Events per chunk", advanced=True,
        help="Events per npsim/eicrecon job; the productions used 1000."),
    Opt("compress", "int", None, 505, label="ROOT compression", advanced=True,
        help="Of sim and reco files: 505 is ZSTD level 5, about a third smaller than the tools' ZLIB 1; "
             "0 keeps their output as is."),
    Opt("genie_compound", "str", None, "off", choices=["off", "on"], label="GENIE compound-nucleus fate", advanced=True,
        help="hN compound-nucleus / pre-equilibrium fate of the FSI samples (off: GENIE R-3_08_00 withdrew it as "
             "experimental). Set in a copy of the override inside the run, recorded in run.json."),
    Opt("beam_nucleon", "bool", None, True, label="Per-nucleon beam proton", advanced=True,
        help="Added after the afterburner; EICrecon's InclusiveKinematics* and HadronicFinalState need it. Roman Pot / "
             "OMD matrix reconstruction is not valid for Z/A = 1/2 beams either way."),
    Opt("gen_args", "str", None, None, label="Extra generator flags", advanced=True,
        help="Appended to the generator command line of the sample, e.g. -F 0.22."),
    Opt("lumi_fb", "float", None, 10.0, label="Luminosity per nucleus", unit="fb^-1", advanced=True,
        help="Of the yields in the figures. 10 fb^-1 per nucleus is A x 10 fb^-1 per nucleon, the ePIC convention for eA."),
]

SECTIONS = {"run": RUN, "generate": GENERATE, "fsi": FSI, "program": PROGRAM, "eic": EIC}


def to_json():
    """The schema as plain data, for the GUI."""
    out = {"_channels": [dict(name=n, label=l, hN=h, hA=a, scalable=n not in SWITCH_CHANNELS) for n, l, h, a in CHANNELS],
           "_targets": {m: [dict(A=a, Z=z, name=n) for a, z, n in t] for m, t in TARGETS.items()},
           "_programs": {k: dict(label=v.label, beam=v.beam, note=v.note) for k, v in PROGRAMS.items()},
           "_eic": {"samples": {k: dict(pass_fraction=v.pass_fraction, deexcite=list(v.deexcite),
                                        production={i: eic_production(k, i) for i in EIC_IONS},
                                        test={i: eic_test_slices(k, i) for i in EIC_IONS})
                                for k, v in EIC_SAMPLES.items()},
                    "mb_per_event": EIC_MB_PER_EVENT, "core_s_per_event": EIC_CORE_S_PER_EVENT}}
    for name, opts in SECTIONS.items():
        out[name] = [dict(key=o.key, kind=o.kind, default=o.default, help=o.help,
                          label=o.label or o.key, choices=o.choices, models=o.models,
                          unit=o.unit, advanced=o.advanced, tab=o.tab, choice_labels=o.choice_labels,
                          min=o.ge if o.ge is not None else o.gt, max=o.le) for o in opts]
    return out
