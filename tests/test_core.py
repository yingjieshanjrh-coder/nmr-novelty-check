import pytest
from rdkit import Chem

from nmrcheck.chem import identify_smiles, parse_reaction, symmetry_classes
from nmrcheck.match import agreement, align, overlaps
from nmrcheck.predict import predict


def test_identity_is_canonical():
    a = identify_smiles("C1=CC=CC=C1")
    b = identify_smiles("c1ccccc1")
    assert a == b
    assert a[0] == "UHOVQNZJYSORNB-UHFFFAOYSA-N"
    assert a[2] == "C6H6"
    assert identify_smiles("not-a-smiles") is None
    assert identify_smiles("") is None


def test_symmetry_classes():
    assert len(symmetry_classes(Chem.MolFromSmiles("c1ccccc1"))) == 1
    assert len(symmetry_classes(Chem.MolFromSmiles("Cc1ccccc1"))) == 5


def test_reaction_parsing():
    assert parse_reaction("CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC") == \
        [["CC(=O)O", "OCC"], ["OS(=O)(=O)O"], ["CC(=O)OCC"]]
    with pytest.raises(ValueError):
        parse_reaction("CC>>CC>O")


def test_align_shares_peaks_when_signals_coincide():
    pairs, unexplained = align([129.0, 128.9, 21.0], [128.9, 21.3])
    assert dict(pairs) == {0: 0, 1: 0, 2: 1}
    assert unexplained == []


def test_align_reports_unexplained_peaks():
    pairs, unexplained = align([128.5], [128.4, 55.0, 20.1])
    assert pairs == [(0, 0)]
    assert unexplained == [1, 2]


def test_agreement_scores():
    good = agreement("13C", [128.5], [128.4])
    bad = agreement("1H", [7.36], [0.82, 0.97])
    assert good["score"] > 0.9 and good["unexplained"] == 0
    assert bad["score"] < 0.05
    assert agreement("13C", [], [1.0])["score"] == 0.0


def test_overlaps_window():
    assert overlaps("13C", [128.5], [128.9, 137.8]) == [(128.5, 128.9)]
    assert overlaps("1H", [7.36], [7.26]) == []


def test_prediction_uses_trained_environments(index):
    benzene = Chem.MolFromSmiles("c1ccccc1")
    peaks = predict(benzene, "13C", index.lookup_env)
    assert len(peaks) == 1
    assert peaks[0]["shift"] == 128.5
    assert peaks[0]["count"] == 6
    assert peaks[0]["radius"] >= 2
    proton = predict(benzene, "1H", index.lookup_env)
    assert proton[0]["shift"] == 7.36 and proton[0]["count"] == 6


def test_prediction_backs_off_to_smaller_environments(index):
    # propylbenzene is not in the fixture; its ring carbons are still covered
    mol = Chem.MolFromSmiles("CCCc1ccccc1")
    peaks = predict(mol, "13C", index.lookup_env)
    assert len(peaks) == 7
    aromatic = [p for p in peaks if p["shift"] > 100]
    assert len(aromatic) == 4
    assert all(p["radius"] >= 0 for p in aromatic)


def test_prediction_fallback_without_table():
    mol = Chem.MolFromSmiles("CC#N")
    peaks = predict(mol, "13C", lambda *a: None)
    assert {p["radius"] for p in peaks} == {-1}
    assert sorted(p["shift"] for p in peaks) == [30.0, 80.0]
