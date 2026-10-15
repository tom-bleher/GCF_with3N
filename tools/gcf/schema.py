"""Options of the pipeline stages, in one place.

Every option the front end knows about is described here once: its type,
default, the gcf_generate / gcf_fsi flag it maps to, which generator models it
applies to, and a help text. The TOML config validation, the command lines of
the stage binaries and the GUI form are all built from these tables.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional

MODELS = ["3N", "2N", "MF"]


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


GENERATE = [
    Opt("model", "str", "--model", "3N", choices=MODELS, label="Generator model",
        help="3N: struck nucleon from an SRC triplet. 2N: from an SRC pair. "
             "MF: single mean-field nucleon from a Fermi gas."),
    Opt("beam", "float", "--beam", 6.0, label="Beam energy", unit="GeV"),
    Opt("events", "int", "--events", 10000, label="Events",
        help="Number of accepted (weight > 0) events in total, over all jobs."),
    Opt("target", "target", None, [12, 6], label="Target (A, Z)"),
    Opt("cs", "str", "--cs", "cc1", choices=["cc1", "cc2", "onshell"],
        label="eN cross section", help="De Forest cc1 (default), cc2, or on-shell."),
    Opt("sigma_cm", "float", "--sigma-cm", None, label="sigma_CM", unit="GeV/c",
        help="CM momentum width per component. Default: 0.15 (3N), the "
             "nucleus value (2N)."),
    Opt("q2", "range", "--q2", None, label="Q^2 window", unit="GeV^2",
        help="Generation-time Q^2 window MIN:MAX. Default: open (3N), 1:10 (2N/MF)."),
    # 3N
    Opt("interaction", "int", "--interaction", 1, choices=list(range(1, 11)),
        models=["3N"], label="3N density table",
        help="1-5: AV18, AV18+UIX, AV18+TM, CD-Bonn, CD-Bonn+TM (3He Faddeev); "
             "6-10: the same, deuteron (pD) channel only."),
    Opt("theta_e", "range", "--theta-e", [5.0, 40.0], models=["3N"],
        label="Electron angle window", unit="deg"),
    Opt("kF", "float", "--kF", None, models=["3N"], label="k_F cutoff", unit="GeV/c",
        help="Internal momenta k > k_F (default 0.25); also sets sigma_CM = "
             "sqrt(3/5) k_F unless sigma_CM is given."),
    Opt("cm", "bool", "--no-cm", True, models=["3N"], invert=True,
        label="CM motion of the triplet"),
    # 2N / MF
    Opt("nn", "str", "--nn", None, models=["2N", "MF"], label="NN interaction",
        help="NN interaction for the 2N contacts (default AV18)."),
    Opt("estar", "float", "--estar", None, models=["2N", "MF"], label="E*", unit="GeV",
        help="Excitation energy of the residual nucleus."),
    Opt("kcut", "float", "--kcut", None, models=["2N"], label="k_rel cutoff", unit="GeV/c",
        help="Relative-momentum floor of the SRC pair (default 0.25)."),
    Opt("onset", "str", "--onset", "fixed", choices=["fixed", "local"], models=["2N"],
        label="SRC onset", help="fixed: k_rel > cutoff. local: k_rel > k_F(r) at the vertex."),
    Opt("prel_is", "bool", "--prel-is", False, models=["2N"],
        label="Importance-sample |p_rel|", help="Lower weight variance; same expectation values."),
    Opt("fermi_gas", "str", "--fermi-gas", "global", choices=["global", "local"],
        models=["MF"], label="Fermi gas"),
    Opt("kF_fg", "float", "--kF-fg", None, models=["MF"], label="Global Fermi momentum",
        unit="GeV/c", help="Default 0.25."),
    Opt("rad", "bool", "--rad", False, models=["2N", "MF"], label="Radiation (peaking approx.)"),
    Opt("coulomb", "bool", "--coulomb", False, models=["2N", "MF"], label="Coulomb correction",
        help="Applied at generation, i.e. before FSI."),
    Opt("phase_space", "str", "--phase-space", None, models=["2N", "MF"],
        label="Phase-space file", help="gcfGenerator phase-space limits file."),
]

FSI = [
    Opt("input", "str", "--input", None, label="Input event file",
        help="Only for FSI-only runs: any file in the GCF event-file format."),
    Opt("model", "str", "--fsi-model", "hN", choices=["hN", "hA"], label="Cascade model",
        help="hN: full hadron-nucleon cascade. hA: effective hadron-nucleus model (faster)."),
    Opt("transport", "str", None, "shared", choices=["shared", "independent"],
        label="Transport",
        help="shared: all nucleons in one GENIE record (absorption depletes the "
             "remnant for the others). independent: one record each."),
    Opt("vertex_power", "str", "--vertex-power", "auto", choices=["auto", "1", "2", "3"],
        label="Vertex weighting (if the file has none)",
        help="rho^n(r) r^2; auto picks n = number of nucleons (max 3)."),
    Opt("target", "target", None, None, label="Target (A, Z) for files without one"),
    Opt("keep", "str", "--keep", None, label="Input branches to keep (regex)",
        help="Default: all."),
]

RUN = [
    Opt("name", "str", None, "run", label="Run name",
        help="Output files are <output_dir>/<name>/<name>_pwia.root and _fsi.root."),
    Opt("output_dir", "str", None, "runs", label="Output directory",
        help="Relative paths are relative to the config file's directory."),
    Opt("stages", "list", None, ["generate", "fsi"], label="Stages",
        choices=["generate", "fsi"]),
    Opt("jobs", "int", None, 1, label="Parallel jobs"),
    Opt("seed", "int", None, 0, label="Base seed",
        help="0 = draw one and record it in the saved config. Job i uses seed + i."),
    Opt("keep_chunks", "bool", None, False, label="Keep per-job files"),
]

SECTIONS = {"run": RUN, "generate": GENERATE, "fsi": FSI}


def to_json():
    """The schema as plain data, for the GUI."""
    out = {}
    for name, opts in SECTIONS.items():
        out[name] = [dict(key=o.key, kind=o.kind, default=o.default, help=o.help,
                          label=o.label or o.key, choices=o.choices, models=o.models,
                          unit=o.unit) for o in opts]
    return out
