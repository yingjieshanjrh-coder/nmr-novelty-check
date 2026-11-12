import json

import pytest

from conftest import EXAMPLE, by_ppm
from peakassign.inputs import Compound, Peak, Targets, load_peaks, load_targets
from peakassign.pipeline import run
from peakassign.report import ALERT, to_json, to_markdown, to_text

EXPECTED = {
    8.05: "product", 7.55: "product", 7.44: "product", 4.38: "product",
    1.39: "product", 3.72: "reactant", 1.24: "reactant",
    7.26: "solvent", 1.56: "solvent", 3.5: "unknown",
}


@pytest.fixture
def example(predict):
    targets = load_targets(EXAMPLE / "targets.yaml")
    peaks = load_peaks(EXAMPLE / "measured_peaks.csv")
    return run(targets, peaks, index=object(), predict=predict)


def test_example_is_sorted_into_the_right_bins(example):
    assert {ppm: p.category for ppm, p in by_ppm(example).items()} == EXPECTED
    assert example.unit == pytest.approx(1.01, abs=0.02)
    acid, ethanol = example.reactants
    assert acid.status == "masked"
    assert ethanol.status == "detected" and ethanol.alert
    assert ethanol.ratio == pytest.approx(0.3, abs=0.02)
    assert not example.warnings


def test_example_warns_about_the_stray_peak(example):
    assert example.unknown_alert
    unknown = example.category("unknown")
    assert len(unknown) == 1
    assert unknown[0].ppm == 3.50
    assert unknown[0].h == pytest.approx(2.0, abs=0.05)


def test_markdown_report(example):
    text = to_markdown(example)
    for heading in ("## 1. 产物特征峰", "## 2. 未反应完的原料",
                    "## 3. 溶剂与已知杂质", "## 4. ⚠️ 未知副产物预警",
                    "## 5. 冲突与重叠"):
        assert heading in text
    assert text.count(ALERT) == 2
    assert "**3.50**" in text
    assert "原料残留较多" in text
    assert "methanol" in text          # hint only, the peak stays unknown
    assert "diethyl ether" not in text  # a singlet cannot be the ether quartet
    assert "30 mol%" in text


def test_json_and_text(example):
    doc = json.loads(to_json(example))
    assert doc["unknown_alert"] is True
    assert doc["unknown"] == [{"ppm": 3.5, "h": pytest.approx(1.98, abs=0.03)}]
    assert len(doc["peaks"]) == 10
    assert doc["reactants"][1]["molar_ratio"] == pytest.approx(0.3, abs=0.02)
    assert ALERT in to_text(example)


def test_report_without_unknowns(predict):
    targets = load_targets(EXAMPLE / "targets.yaml")
    peaks = [p for p in load_peaks(EXAMPLE / "measured_peaks.csv")
             if p.ppm != 3.50]
    result = run(targets, peaks, index=object(), predict=predict)
    text = to_markdown(result)
    assert ALERT not in text
    assert "没有发现积分超过" in text


ETOAC = [Peak(4.12, 0.4, "q"), Peak(2.05, 0.6, "s"), Peak(1.26, 0.6, "t")]


def acetone_targets(**kwargs):
    return Targets(product=Compound("CC(C)=O", "acetone"), reactants=[],
                   solvent="CDCl3", **kwargs)


def test_undeclared_lab_solvent_is_flagged_with_a_hint(predict):
    peaks = [Peak(2.17, 6.0, "s")] + ETOAC
    result = run(acetone_targets(), peaks, index=object(), predict=predict)
    cats = {ppm: p.category for ppm, p in by_ppm(result).items()}
    assert cats[2.05] == "unknown"
    assert "ethyl acetate" in to_markdown(result)


def test_declared_lab_solvent_is_recognised(predict):
    peaks = [Peak(2.17, 6.0, "s")] + ETOAC
    result = run(acetone_targets(extra_solvents=["ethyl acetate"]), peaks,
                 index=object(), predict=predict)
    cats = {ppm: p.category for ppm, p in by_ppm(result).items()}
    assert cats == {2.17: "product", 4.12: "solvent", 2.05: "solvent",
                    1.26: "solvent"}
    assert not result.unknown_alert


def test_declared_solvent_needs_all_its_signals(predict):
    # only the 2.05 singlet: one of three ethyl acetate signals is not enough
    peaks = [Peak(2.17, 6.0, "s"), Peak(2.05, 0.6, "s")]
    result = run(acetone_targets(extra_solvents=["ethyl acetate"]), peaks,
                 index=object(), predict=predict)
    assert by_ppm(result)[2.05].category == "unknown"
    assert result.impurities[0].status == "doubtful"


def test_user_supplied_impurity_line(predict):
    peaks = [Peak(2.17, 6.0, "s"), Peak(0.07, 0.9, "s"), Peak(5.01, 0.8, "s")]
    targets = acetone_targets(extra_solvents=["silicone grease"],
                              extra_peaks=[{"name": "BHT OH", "ppm": 5.01}])
    result = run(targets, peaks, index=object(), predict=predict)
    cats = {ppm: p.category for ppm, p in by_ppm(result).items()}
    assert cats[0.07] == cats[5.01] == "solvent"


def test_missing_index_is_mentioned(predict):
    targets = load_targets(EXAMPLE / "targets.yaml")
    peaks = load_peaks(EXAMPLE / "measured_peaks.csv")
    result = run(targets, peaks, index=None, predict=predict)
    assert any("index.sqlite" in w for w in result.warnings)
