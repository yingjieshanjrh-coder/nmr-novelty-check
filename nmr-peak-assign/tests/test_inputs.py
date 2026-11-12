import pytest

from conftest import EXAMPLE
from peakassign import impurities
from peakassign.inputs import InputError, load_peaks, load_targets


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_example_files_load():
    targets = load_targets(EXAMPLE / "targets.yaml")
    assert targets.product.smiles == "CCOC(=O)c1ccccc1"
    assert targets.product.name == "苯甲酸乙酯"
    assert [r.smiles for r in targets.reactants] == ["OC(=O)c1ccccc1", "CCO"]
    assert targets.solvent == "CDCl3"
    peaks = load_peaks(EXAMPLE / "measured_peaks.csv")
    assert len(peaks) == 10
    assert peaks[0].ppm == 8.05 and peaks[0].multiplicity == "d"


def test_plain_smiles_and_solvent_alias(tmp_path):
    path = write(tmp_path, "t.yaml",
                 "product: CCO\nreactants: [CC=O]\nsolvent: dmso\n"
                 "extra_solvents: [EtOAc, 丙酮]\n"
                 "reference: {ppm: 3.44, protons: 2}\n")
    targets = load_targets(path)
    assert targets.solvent == "DMSO-d6"
    assert targets.reactants[0].smiles == "CC=O"
    assert targets.extra_solvents == ["ethyl acetate", "acetone"]
    assert targets.reference == {"ppm": 3.44, "protons": 2.0}


@pytest.mark.parametrize("text, word", [
    ("reactants: []\nsolvent: CDCl3\n", "product"),
    ("product: CCO\nsolvent: CDCl3\n", "reactants"),
    ("product: CCO\nreactants: []\nsolvent: water-d9\n", "solvent"),
    ("product: CCO\nreactants: []\nsolvent: CDCl3\nextra_solvents: [xenon]\n",
     "xenon"),
    ("product: CCO\nreactants: [{name: x}]\nsolvent: CDCl3\n", "reactants[0]"),
])
def test_bad_targets(tmp_path, text, word):
    with pytest.raises(InputError, match=word.replace("[", r"\[")):
        load_targets(write(tmp_path, "t.yaml", text))


def test_peak_list_details(tmp_path):
    path = write(tmp_path, "p.csv",
                 "﻿Shift,Area\n1.20,3\n\n7.30, 5.5\n")
    peaks = load_peaks(path)
    assert [(p.ppm, p.integral, p.multiplicity) for p in peaks] == \
        [(7.3, 5.5, ""), (1.2, 3.0, "")]


@pytest.mark.parametrize("text", [
    "", "ppm,height\n1.0,2.0\n", "ppm,integral\n1.0,abc\n",
    "ppm,integral\n1.0,-2\n", "ppm,integral\n",
])
def test_bad_peak_lists(tmp_path, text):
    with pytest.raises(InputError):
        load_peaks(write(tmp_path, "p.csv", text))


def test_solvent_lines():
    lines = {l["kind"]: l["ppm"] for l in impurities.solvent_lines("CDCl3")}
    assert lines == {"residual": 7.26, "water": 1.56}
    lines = {l["kind"]: l["ppm"] for l in impurities.solvent_lines("DMSO-d6")}
    assert lines == {"residual": 2.50, "water": 3.33}
    # in D2O the residual HDO signal is the water signal
    assert len(impurities.solvent_lines("D2O")) == 1


def test_trace_lookup():
    assert impurities.resolve_trace("乙酸乙酯", "CDCl3") == "ethyl acetate"
    assert impurities.resolve_trace("DCM", "DMSO-d6") == "dichloromethane"
    assert impurities.resolve_trace("acetone", "CD3OD") is None
    names = [hit[0] for hit in impurities.nearby_trace("CDCl3", 3.50)]
    assert "methanol" in names
