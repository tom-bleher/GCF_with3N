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
    ("nucleon.compound",   "compound nucleus",     True,  True),
    ("pion.elastic",       "elastic scattering",   True,  False),
    ("pion.inelastic",     "inelastic scattering", True,  True),
    ("pion.cex",           "charge exchange",      True,  True),
    ("pion.absorption",    "absorption",           True,  True),
    ("pion.piprod",        "pion production",      False, True),
]


@dataclass
class Opt:
    key: str                     # name in the TOML config / GUI
    kind: str                    # float, int, str, bool, range, target
    flag: Optional[str]          # binary flag (None: handled by the runner)
    default: Any = None          # None = not passed (binary default applies)
    help: str = ""
    label: str = ""
    choices: Optional[List[Any]] = None
    models: List[str] = field(default_factory=lambda: list(MODELS))
    invert: bool = False         # bool flag that is passed when the value is False
    unit: str = ""
    advanced: bool = False       # GUI: shown under "Advanced", folded by default


GENERATE = [
    Opt("model", "str", "--model", "3N", choices=MODELS, label="Generator model",
        help="3N: struck nucleon from an SRC triplet. 2N: from an SRC pair. "
             "MF: single mean-field nucleon from a Fermi gas."),
    Opt("beam", "float", "--beam", 6.0, label="Beam energy", unit="GeV"),
    Opt("events", "int", "--events", 10000, label="Events",
        help="Number of accepted (weight > 0) events in total, over all jobs."),
    Opt("target", "target", None, [12, 6], label="Target (A, Z)"),
    Opt("cs", "str", "--cs", "cc1", choices=["cc1", "cc2", "onshell"],
        label="eN cross section", help="De Forest cc1 (default), cc2, or on-shell.", advanced=True),
    Opt("sigma_cm", "float", "--sigma-cm", None, models=["3N", "2N"], label="sigma_CM", unit="GeV/c",
        help="CM momentum width of the SRC cluster per component. Default: 0.15 (3N), "
             "the nucleus table value (2N).", advanced=True),
    Opt("q2", "range", "--q2", None, label="Q^2 window", unit="GeV^2",
        help="Generation-time Q^2 window MIN:MAX. Default: open (3N), 1:10 (2N/MF).", advanced=True),
    # 3N
    Opt("wave_function", "str", "--wave-function", "AV18", choices=WAVE_FUNCTIONS_3N,
        models=["3N"], label="3N wave function",
        help="Nucleon interaction behind the three-body density (AV18, AV18+UIX, AV18+TM, "
             "CDBonn, CDBonn+TM). The 3He Faddeev density is used for the SRC triplet in "
             "every nucleus. The -pD variants keep only the deuteron (pD) channel."),
    Opt("theta_e", "range", "--theta-e", [5.0, 40.0], models=["3N"],
        label="Electron angle window", unit="deg", advanced=True),
    Opt("kF", "float", "--kF", None, models=["3N"], label="k_F cutoff", unit="GeV/c",
        help="Internal momenta k > k_F (default 0.25); also sets sigma_CM = "
             "sqrt(3/5) k_F unless sigma_CM is given.", advanced=True),
    Opt("cm", "bool", "--no-cm", True, models=["3N"], invert=True,
        label="CM motion of the triplet", advanced=True),
    # 2N / MF
    Opt("nn", "str", "--nn", None, models=["2N"], label="NN interaction",
        choices=NN_INTERACTIONS,
        help="NN interaction behind the pair wave function and contacts (default AV18)."),
    Opt("estar", "float", "--estar", None, models=["2N"], label="E*", unit="GeV",
        help="Excitation energy of the residual nucleus.", advanced=True),
    Opt("kcut", "float", "--kcut", None, models=["2N"], label="k_rel cutoff", unit="GeV/c",
        help="Relative-momentum floor of the SRC pair (default 0.25).", advanced=True),
    Opt("onset", "str", "--onset", "fixed", choices=["fixed", "local"], models=["2N"],
        label="SRC onset", help="fixed: k_rel > cutoff. local: k_rel > k_F(r) at the vertex.", advanced=True),
    Opt("prel_is", "bool", "--prel-is", False, models=["2N"],
        label="Importance-sample |p_rel|", help="Lower weight variance; same expectation values.", advanced=True),
    Opt("fermi_gas", "str", "--fermi-gas", "global", choices=["global", "local"],
        models=["MF"], label="Fermi gas",
        help="global: one Fermi momentum for the whole nucleus. local: k_F(r) from the density."),
    Opt("kF_fg", "float", "--kF-fg", None, models=["MF"], label="Global Fermi momentum",
        unit="GeV/c", help="Default 0.25.", advanced=True),
    Opt("rad", "bool", "--rad", False, models=["2N", "MF"], label="Radiation (peaking approx.)", advanced=True),
    Opt("coulomb", "bool", "--coulomb", False, models=["2N", "MF"], label="Coulomb correction",
        help="Applied at generation, i.e. before FSI.", advanced=True),
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
    Opt("nucleon_mfp_scale", "float", None, None, label="Nucleon mean-free-path scale",
        help="Multiplies the nucleon mean free path (> 1: fewer interactions).", advanced=True),
    Opt("pion_mfp_scale", "float", None, None, label="Pion mean-free-path scale",
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
    Opt("keep_chunks", "bool", None, False, label="Keep per-job files", advanced=True),
]

SECTIONS = {"run": RUN, "generate": GENERATE, "fsi": FSI}


def to_json():
    """The schema as plain data, for the GUI."""
    out = {"_channels": [dict(name=n, label=l, hN=h, hA=a) for n, l, h, a in CHANNELS],
           "_targets": {m: [dict(A=a, Z=z, name=n) for a, z, n in t] for m, t in TARGETS.items()}}
    for name, opts in SECTIONS.items():
        out[name] = [dict(key=o.key, kind=o.kind, default=o.default, help=o.help,
                          label=o.label or o.key, choices=o.choices, models=o.models,
                          unit=o.unit, advanced=o.advanced) for o in opts]
    return out
