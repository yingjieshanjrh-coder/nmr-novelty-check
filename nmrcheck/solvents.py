"""Residual solvent signals (ppm), after Fulmer et al., Organometallics 2010, 29, 2176."""

RESIDUAL = {
    "CDCl3": {"1H": [7.26], "13C": [77.16]},
    "DMSO-d6": {"1H": [2.50], "13C": [39.52]},
    "CD3OD": {"1H": [3.31], "13C": [49.00]},
    "D2O": {"1H": [4.79], "13C": []},
    "CD3COCD3": {"1H": [2.05], "13C": [29.84, 206.26]},
    "CD3CN": {"1H": [1.94], "13C": [1.32, 118.26]},
    "C6D6": {"1H": [7.16], "13C": [128.06]},
    "CD2Cl2": {"1H": [5.32], "13C": [53.84]},
}

ALIASES = {
    "cdcl3": "CDCl3", "chloroform": "CDCl3", "chloroform-d": "CDCl3",
    "dmso": "DMSO-d6", "dmso-d6": "DMSO-d6",
    "cd3od": "CD3OD", "methanol-d4": "CD3OD", "meod": "CD3OD",
    "d2o": "D2O",
    "acetone-d6": "CD3COCD3", "cd3cocd3": "CD3COCD3", "acetone": "CD3COCD3",
    "cd3cn": "CD3CN", "acetonitrile-d3": "CD3CN",
    "c6d6": "C6D6", "benzene-d6": "C6D6",
    "cd2cl2": "CD2Cl2", "dcm-d2": "CD2Cl2",
}


def resolve(name):
    return ALIASES.get((name or "").strip().lower())
