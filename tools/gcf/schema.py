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
        choices=["generate", "fsi"]),
    Opt("jobs", "int", None, 1, label="Parallel jobs",
        help="Independent jobs run at the same time; their outputs are merged. "
             "Each needs one CPU core."),
    Opt("seed", "int", None, 0, label="Base seed",
        help="0 = draw one and record it in the saved config. Job i uses seed + i.", advanced=True),
    Opt("note", "str", None, None, label="Note", advanced=True,
        help="What this run is for. Shown in the list of runs and searched by its filter."),
    Opt("chunks", "int", None, None, label="Split into jobs", advanced=True,
        help="More jobs than run at a time means an interruption loses less: a resume reruns only the jobs "
             "that didn't finish. Default: one per million events, and at least Parallel jobs. "
             "Job i uses seed + i, so this is part of what makes a run reproducible."),
    Opt("keep_chunks", "bool", None, False, label="Keep per-job files", advanced=True),
]

SECTIONS = {"run": RUN, "generate": GENERATE, "fsi": FSI}


def to_json():
    """The schema as plain data, for the GUI."""
    out = {"_channels": [dict(name=n, label=l, hN=h, hA=a, scalable=n not in SWITCH_CHANNELS) for n, l, h, a in CHANNELS],
           "_targets": {m: [dict(A=a, Z=z, name=n) for a, z, n in t] for m, t in TARGETS.items()}}
    for name, opts in SECTIONS.items():
        out[name] = [dict(key=o.key, kind=o.kind, default=o.default, help=o.help,
                          label=o.label or o.key, choices=o.choices, models=o.models,
                          unit=o.unit, advanced=o.advanced, tab=o.tab,
                          min=o.ge if o.ge is not None else o.gt, max=o.le) for o in opts]
    return out
