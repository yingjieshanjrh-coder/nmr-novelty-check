"""Residual solvent, water and common laboratory solvent signals (1H, ppm).

Residual solvent shifts and solvent aliases come from the parent project
(nmrcheck.solvents). Water and the trace-solvent table are transcribed from
Fulmer et al., Organometallics 2010, 29, 2176 (Table 1).
"""
from .parent import ensure_importable

ensure_importable()
from nmrcheck import solvents as _parent  # noqa: E402

WATER = {
    "CDCl3": 1.56, "DMSO-d6": 3.33, "CD3OD": 4.87, "D2O": 4.79,
    "CD3COCD3": 2.84, "CD3CN": 2.13, "C6D6": 0.40, "CD2Cl2": 1.52,
}

RESIDUAL_NAME = {
    "CDCl3": "CHCl3", "DMSO-d6": "DMSO-d5", "CD3OD": "CHD2OD", "D2O": "HDO",
    "CD3COCD3": "acetone-d5", "CD3CN": "CHD2CN", "C6D6": "C6D5H",
    "CD2Cl2": "CHDCl2",
}

# name -> [(ppm, multiplicity, protons)]; protons None = no fixed ratio
TRACE = {
    "CDCl3": {
        "acetone": [(2.17, "s", 6)],
        "acetonitrile": [(2.10, "s", 3)],
        "acetic acid": [(2.10, "s", 3)],
        "benzene": [(7.36, "s", 6)],
        "cyclohexane": [(1.43, "s", 12)],
        "dichloromethane": [(5.30, "s", 2)],
        "diethyl ether": [(3.48, "q", 4), (1.21, "t", 6)],
        "DMF": [(8.02, "s", 1), (2.96, "s", 3), (2.88, "s", 3)],
        "DMSO": [(2.62, "s", 6)],
        "1,4-dioxane": [(3.71, "s", 8)],
        "ethanol": [(3.72, "q", 2), (1.25, "t", 3)],
        "ethyl acetate": [(4.12, "q", 2), (2.05, "s", 3), (1.26, "t", 3)],
        "n-hexane": [(1.26, "m", 8), (0.88, "t", 6)],
        "methanol": [(3.49, "s", 3)],
        "2-propanol": [(4.04, "sep", 1), (1.22, "d", 6)],
        "pyridine": [(8.62, "m", 2), (7.68, "m", 1), (7.29, "m", 2)],
        "THF": [(3.76, "m", 4), (1.85, "m", 4)],
        "toluene": [(7.25, "m", 2), (7.17, "m", 3), (2.36, "s", 3)],
        "triethylamine": [(2.53, "q", 6), (1.03, "t", 9)],
        "grease": [(1.26, "br", None), (0.86, "m", None)],
        "silicone grease": [(0.07, "s", None)],
        "TMS": [(0.00, "s", None)],
    },
    "DMSO-d6": {
        "acetone": [(2.09, "s", 6)],
        "acetonitrile": [(2.07, "s", 3)],
        "acetic acid": [(1.91, "s", 3)],
        "benzene": [(7.37, "s", 6)],
        "chloroform": [(8.32, "s", 1)],
        "cyclohexane": [(1.40, "s", 12)],
        "dichloromethane": [(5.76, "s", 2)],
        "diethyl ether": [(3.38, "q", 4), (1.09, "t", 6)],
        "DMF": [(7.95, "s", 1), (2.89, "s", 3), (2.73, "s", 3)],
        "1,4-dioxane": [(3.57, "s", 8)],
        "ethanol": [(3.44, "q", 2), (1.06, "t", 3)],
        "ethyl acetate": [(4.03, "q", 2), (1.99, "s", 3), (1.17, "t", 3)],
        "n-hexane": [(1.25, "m", 8), (0.86, "t", 6)],
        "methanol": [(3.16, "s", 3)],
        "2-propanol": [(3.78, "sep", 1), (1.04, "d", 6)],
        "pyridine": [(8.58, "m", 2), (7.79, "m", 1), (7.39, "m", 2)],
        "THF": [(3.60, "m", 4), (1.76, "m", 4)],
        "toluene": [(7.25, "m", 2), (7.18, "m", 3), (2.30, "s", 3)],
        "triethylamine": [(2.43, "q", 6), (0.93, "t", 9)],
        "TMS": [(0.00, "s", None)],
    },
}

TRACE_ALIASES = {
    "丙酮": "acetone", "乙腈": "acetonitrile", "mecn": "acetonitrile",
    "乙酸": "acetic acid", "醋酸": "acetic acid", "acoh": "acetic acid",
    "苯": "benzene", "环己烷": "cyclohexane",
    "氯仿": "chloroform", "chcl3": "chloroform",
    "dcm": "dichloromethane", "二氯甲烷": "dichloromethane",
    "ch2cl2": "dichloromethane",
    "乙醚": "diethyl ether", "ether": "diethyl ether", "et2o": "diethyl ether",
    "二甲亚砜": "DMSO", "dmso": "DMSO", "dmf": "DMF",
    "二氧六环": "1,4-dioxane", "dioxane": "1,4-dioxane",
    "乙醇": "ethanol", "etoh": "ethanol",
    "乙酸乙酯": "ethyl acetate", "ea": "ethyl acetate",
    "etoac": "ethyl acetate",
    "正己烷": "n-hexane", "己烷": "n-hexane", "hexane": "n-hexane",
    "hexanes": "n-hexane", "石油醚": "n-hexane", "pe": "n-hexane",
    "甲醇": "methanol", "meoh": "methanol",
    "异丙醇": "2-propanol", "isopropanol": "2-propanol", "ipa": "2-propanol",
    "吡啶": "pyridine", "四氢呋喃": "THF", "thf": "THF",
    "甲苯": "toluene", "三乙胺": "triethylamine", "et3n": "triethylamine",
    "tea": "triethylamine", "油脂": "grease", "真空脂": "grease",
    "硅脂": "silicone grease", "tms": "TMS",
}

RESIDUAL_WINDOW = 0.03
WATER_WINDOW = 0.15
TRACE_WINDOW = 0.04


def resolve_solvent(name):
    """Canonical name of a deuterated solvent, or None if it is not known."""
    return _parent.resolve(name)


def known_solvents():
    return sorted(_parent.RESIDUAL)


def solvent_lines(solvent):
    """Signals that are always looked for: residual solvent and water."""
    lines = []
    for ppm in _parent.RESIDUAL[solvent]["1H"]:
        lines.append({"name": f"{RESIDUAL_NAME[solvent]}（{solvent} 残留峰）",
                      "kind": "residual", "ppm": ppm,
                      "window": RESIDUAL_WINDOW})
    water = WATER.get(solvent)
    if water is not None and all(abs(water - line["ppm"]) > 0.01
                                 for line in lines):
        lines.append({"name": f"H2O（{solvent} 中的水峰）", "kind": "water",
                      "ppm": water, "window": WATER_WINDOW})
    return lines


def resolve_trace(name, solvent):
    """Canonical trace-solvent name if it is tabulated for `solvent`."""
    table = TRACE.get(solvent, {})
    key = (name or "").strip()
    by_lower = {k.lower(): k for k in table}
    key = TRACE_ALIASES.get(key.lower(), key)
    return by_lower.get(key.lower())


def trace_lines(name, solvent):
    return list(TRACE[solvent][name])


def nearby_trace(solvent, ppm, window=TRACE_WINDOW):
    """Tabulated trace solvents with a signal next to `ppm`; used for hints."""
    hits = []
    for name, lines in TRACE.get(solvent, {}).items():
        for shift, mult, _ in lines:
            if abs(shift - ppm) <= window:
                hits.append((name, shift, mult, len(lines)))
    return hits
