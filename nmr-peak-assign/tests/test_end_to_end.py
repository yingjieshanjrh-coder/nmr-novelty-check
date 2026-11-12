"""Runs the example through the real predictor.

Needs RDKit and the index the parent project builds in ../data; skipped
otherwise.
"""
import pytest

pytest.importorskip("rdkit")

from conftest import EXAMPLE, by_ppm  # noqa: E402
from peakassign.cli import main  # noqa: E402
from peakassign.inputs import load_peaks, load_targets  # noqa: E402
from peakassign.parent import open_index  # noqa: E402
from peakassign.pipeline import run  # noqa: E402
from peakassign.predictor import build_species  # noqa: E402
from peakassign.report import ALERT  # noqa: E402


@pytest.fixture(scope="module")
def index():
    idx = open_index()
    if idx is None:
        pytest.skip("parent index not built (run `nmrcheck build` first)")
    yield idx
    idx.close()


def test_predicted_signals_of_ethyl_benzoate(index):
    species = build_species("CCOC(=O)c1ccccc1", "product", "product",
                            index=index)
    assert sum(g.protons for g in species.groups) == 10
    assert sorted(g.protons for g in species.groups) == [1, 2, 2, 2, 3]
    assert not any(g.labile for g in species.groups)
    assert species.atom_map.count(":") == 11


def test_exchangeable_protons_are_marked(index):
    species = build_species("CCO", "reactant1", "reactant", index=index)
    labile = [g for g in species.groups if g.labile]
    assert len(labile) == 1 and labile[0].protons == 1
    assert labile[0].label.startswith("OH")


def test_example_with_real_predictions(index):
    targets = load_targets(EXAMPLE / "targets.yaml")
    peaks = load_peaks(EXAMPLE / "measured_peaks.csv")
    result = run(targets, peaks, index=index)
    cats = {ppm: p.category for ppm, p in by_ppm(result).items()}
    assert cats == {
        8.05: "product", 7.55: "product", 7.44: "product", 4.38: "product",
        1.39: "product", 3.72: "reactant", 1.24: "reactant",
        7.26: "solvent", 1.56: "solvent", 3.5: "unknown",
    }
    assert [p.ppm for p in result.category("unknown")] == [3.50]
    ethanol = result.reactants[1]
    assert ethanol.status == "detected"
    assert ethanol.ratio == pytest.approx(0.3, abs=0.03)


def test_command_line(index, tmp_path, capsys):
    out = tmp_path / "report.md"
    code = main([str(EXAMPLE / "targets.yaml"),
                 str(EXAMPLE / "measured_peaks.csv"),
                 "-o", str(out), "--json", str(tmp_path / "report.json")])
    assert code == 0
    assert ALERT in out.read_text(encoding="utf-8")
    assert ALERT in capsys.readouterr().out
    assert (tmp_path / "report.json").exists()
