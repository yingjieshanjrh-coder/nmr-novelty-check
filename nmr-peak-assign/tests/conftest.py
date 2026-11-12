import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from peakassign.assign import Group, Species  # noqa: E402
from peakassign.inputs import Peak  # noqa: E402

EXAMPLE = ROOT / "examples" / "ethyl_benzoate"

# hand-written predictions so the assignment logic can be tested without
# RDKit or the shift index: smiles -> [(label, ppm, protons, labile)]
BOOK = {
    "CCOC(=O)c1ccccc1": [("ArH o", 8.04, 2, False), ("ArH p", 7.54, 1, False),
                         ("ArH m", 7.43, 2, False), ("OCH2", 4.37, 2, False),
                         ("CH3", 1.38, 3, False)],
    "OC(=O)c1ccccc1": [("COOH", 12.0, 1, True), ("ArH o", 8.12, 2, False),
                       ("ArH p", 7.62, 1, False), ("ArH m", 7.48, 2, False)],
    "CCO": [("OCH2", 3.70, 2, False), ("OH", 2.60, 1, True),
            ("CH3", 1.22, 3, False)],
    "CC(C)=O": [("CH3", 2.17, 6, False)],
    "COc1ccccc1": [("ArH", 7.29, 2, False), ("ArH", 6.93, 3, False),
                   ("OCH3", 3.80, 3, False)],
}


def book_predict(smiles, key, role, name="", index=None, tol=0.08):
    groups = [Group(label=label, shift=ppm, protons=n, window=0.20,
                    tight=0.10, labile=labile)
              for label, ppm, n, labile in BOOK[smiles]]
    return Species(key=key, role=role, name=name or smiles, groups=groups,
                   smiles=smiles)


def species(smiles, key="product", role="product", name=""):
    return book_predict(smiles, key, role, name)


def peaks(*rows):
    return [Peak(ppm, integral, mult[0] if mult else "")
            for ppm, integral, *mult in rows]


def by_ppm(result):
    return {round(p.ppm, 2): p for p in result.peaks}


@pytest.fixture
def predict():
    return book_predict
