import pytest

from conftest import by_ppm, peaks, species
from peakassign.assign import Group, Settings, Species, analyse, fit_species

WATER = {"name": "H2O", "kind": "water", "ppm": 1.56, "window": 0.15}
CHCL3 = {"name": "CHCl3", "kind": "residual", "ppm": 7.26, "window": 0.03}

ESTER = "CCOC(=O)c1ccccc1"
CLEAN = [(8.05, 2.0), (7.55, 1.0), (7.44, 2.0), (4.38, 2.0), (1.39, 3.0)]


def ester():
    return species(ESTER)


def ethanol():
    return species("CCO", "reactant1", "reactant", "ethanol")


def categories(result):
    return {ppm: p.category for ppm, p in by_ppm(result).items()}


def test_clean_product_is_fully_assigned():
    result = analyse(peaks(*CLEAN), ester(), [])
    assert set(categories(result).values()) == {"product"}
    assert result.unit == pytest.approx(1.0)
    assert not result.product.missing
    assert not result.unknown_alert
    assert all(p.parts[0].confidence == "高" for p in result.peaks)


def test_scale_of_raw_integrals_does_not_matter():
    rows = CLEAN + [(3.50, 2.0), (1.56, 0.4)]
    small = analyse(peaks(*rows), ester(), [], lines=[WATER])
    big = analyse(peaks(*[(ppm, area * 137.5) for ppm, area in rows]),
                  ester(), [], lines=[WATER])
    assert categories(small) == categories(big)
    assert big.unit == pytest.approx(137.5)
    assert by_ppm(big)[3.5].h == pytest.approx(2.0)


def test_large_unexplained_peak_raises_the_alert():
    result = analyse(peaks(*CLEAN, (3.50, 2.0, "s")), ester(), [ethanol()])
    unknown = result.category("unknown")
    assert [p.ppm for p in unknown] == [3.50]
    assert unknown[0].left_h == pytest.approx(2.0)
    assert result.unknown_alert


def test_small_unexplained_peak_is_listed_quietly():
    result = analyse(peaks(*CLEAN, (5.90, 0.2)), ester(), [])
    assert categories(result)[5.9] == "minor"
    assert not result.unknown_alert


def test_integral_too_small_for_the_signal_is_a_conflict():
    # 1.0H where the ester CH3 (3H) is predicted
    rows = [r for r in CLEAN if r[0] != 1.39] + [(1.39, 1.0)]
    result = analyse(peaks(*rows), ester(), [])
    peak = by_ppm(result)[1.39]
    assert peak.category == "ambiguous"
    assert "3H" in peak.conflict and not peak.parts
    assert [g.label for g in result.product.missing] == ["CH3"]


def test_predicted_signal_is_used_once():
    result = analyse(peaks(*CLEAN, (4.30, 2.0)), ester(), [])
    cats = categories(result)
    assert cats[4.38] == "product"
    assert cats[4.3] == "unknown"


def test_reactant_ratio_comes_from_its_own_signals():
    rows = CLEAN + [(3.72, 0.6, "q"), (1.24, 0.9, "t")]
    result = analyse(peaks(*rows), ester(), [ethanol()])
    cats = categories(result)
    assert cats[3.72] == cats[1.24] == "reactant"
    report = result.reactants[0]
    assert report.status == "detected"
    assert report.ratio == pytest.approx(0.3, abs=0.01)
    assert report.alert


def test_trace_of_reactant_is_not_called_a_lot():
    rows = CLEAN + [(3.72, 0.1), (1.24, 0.15)]
    report = analyse(peaks(*rows), ester(), [ethanol()]).reactants[0]
    assert report.status == "detected" and not report.alert


def test_peak_next_to_reactant_signal_with_wrong_integral_stays_unknown():
    # 3.50/2.0H sits inside the ethanol OCH2 window, but ethanol is present
    # at 0.3 equivalents and its OCH2 is already accounted for by 3.72
    rows = CLEAN + [(3.72, 0.6), (3.50, 2.0), (1.24, 0.9)]
    result = analyse(peaks(*rows), ester(), [ethanol()])
    cats = categories(result)
    assert cats[3.72] == "reactant"
    assert cats[3.5] == "unknown"


def test_lone_peak_in_a_reactant_window_is_not_enough():
    anisole = species("COc1ccccc1", "reactant1", "reactant", "anisole")
    result = analyse(peaks(*CLEAN, (3.80, 1.5)), ester(), [anisole])
    report = result.reactants[0]
    assert report.status == "doubtful"
    peak = by_ppm(result)[3.8]
    assert peak.category == "unknown"
    assert "anisole" in peak.hints[0]


def test_product_and_reactant_under_one_peak_is_ambiguous():
    # ester CH3 (3H) and ethanol CH3 (0.3 eq -> 0.9H) merged at 1.32
    rows = [r for r in CLEAN if r[0] != 1.39] + [(3.72, 0.6), (1.32, 3.9)]
    result = analyse(peaks(*rows), ester(), [ethanol()])
    peak = by_ppm(result)[1.32]
    assert peak.category == "ambiguous"
    assert {p.role for p in peak.parts} == {"product", "reactant"}
    taken = {p.role: p.h for p in peak.parts}
    assert taken["product"] == pytest.approx(3.0, abs=0.05)
    assert taken["reactant"] == pytest.approx(0.9, abs=0.05)
    assert not result.unknown_alert


def test_reactant_signal_hidden_under_a_product_peak_is_flagged():
    # the reactant shows at 2.60; its other signal sits under the ester CH3,
    # whose integral is still within tolerance
    guest = Species("reactant1", "reactant", "guest", [
        Group("CH2", 2.60, 2, 0.2), Group("CH3", 1.35, 3, 0.2)])
    result = analyse(peaks(*CLEAN, (2.60, 0.2)), ester(), [guest])
    assert result.reactants[0].status == "detected"
    peak = by_ppm(result)[1.39]
    assert peak.category == "ambiguous"
    assert peak.overlaps[0]["name"] == "guest"
    assert peak.overlaps[0]["h"] == pytest.approx(0.3, abs=0.02)


def test_absent_reactant_only_leaves_a_note():
    acid = species("OC(=O)c1ccccc1", "reactant1", "reactant", "acid")
    result = analyse(peaks(*CLEAN), ester(), [acid])
    assert result.reactants[0].status == "masked"
    peak = by_ppm(result)[8.05]
    assert peak.category == "product"
    assert "acid" in peak.also[0]


def test_strict_overlap_moves_such_peaks_to_ambiguous():
    acid = species("OC(=O)c1ccccc1", "reactant1", "reactant", "acid")
    result = analyse(peaks(*CLEAN), ester(), [acid],
                     cfg=Settings(strict_overlap=True))
    cats = categories(result)
    assert cats[8.05] == cats[7.55] == cats[7.44] == "ambiguous"
    assert cats[4.38] == "product"


def test_one_multiplet_can_hold_two_product_signals():
    rows = [(8.05, 2.0), (7.50, 3.0, "m"), (4.38, 2.0), (1.39, 3.0)]
    result = analyse(peaks(*rows), ester(), [])
    peak = by_ppm(result)[7.5]
    assert peak.category == "product"
    assert sorted(peak.parts[0].labels) == ["ArH m", "ArH p"]
    assert not result.product.missing


def test_signal_reported_as_two_peaks_is_put_back_together():
    rows = [r for r in CLEAN if r[0] != 1.39] + [(1.41, 1.5), (1.37, 1.5)]
    result = analyse(peaks(*rows), ester(), [])
    for ppm in (1.41, 1.37):
        peak = by_ppm(result)[ppm]
        assert peak.category == "product"
        assert peak.parts[0].status == "split"
    assert not result.product.missing


def test_extra_area_under_a_product_peak_is_reported():
    rows = [r for r in CLEAN if r[0] != 4.38] + [(4.38, 3.2)]
    result = analyse(peaks(*rows), ester(), [])
    peak = by_ppm(result)[4.38]
    assert peak.category == "product"
    assert peak.parts[0].status == "excess"
    assert peak.parts[0].confidence == "中"
    assert result.hidden == [peak]
    assert peak.left_h == pytest.approx(1.2, abs=0.05)
    assert result.unknown_alert


def test_solvent_and_water_lines():
    rows = CLEAN + [(7.26, 0.2, "s"), (1.56, 0.4, "s")]
    result = analyse(peaks(*rows), ester(), [], lines=[CHCL3, WATER])
    cats = categories(result)
    assert cats[7.26] == cats[1.56] == "solvent"
    assert cats[1.39] == "product"


def test_water_must_look_like_a_singlet():
    rows = CLEAN + [(1.56, 0.8, "q")]
    result = analyse(peaks(*rows), ester(), [], lines=[WATER])
    assert categories(result)[1.56] == "unknown"


def test_exchangeable_proton_is_matched_loosely():
    alcohol = Species("product", "product", "alcohol", [
        Group("CH2", 3.70, 2, 0.2), Group("OH", 2.60, 1, 0.2, labile=True),
        Group("CH3", 1.22, 3, 0.2)])
    rows = [(3.72, 2.0), (1.95, 0.8, "br"), (1.24, 3.0)]
    result = analyse(peaks(*rows), alcohol, [])
    peak = by_ppm(result)[1.95]
    assert peak.category == "product"
    assert peak.parts[0].status == "labile"
    assert peak.parts[0].confidence == "中"


def test_reference_peak_fixes_the_scale():
    rows = [(ppm, area * 5) for ppm, area in CLEAN]
    result = analyse(peaks(*rows), ester(), [],
                     reference={"ppm": 4.38, "protons": 2})
    assert result.scale_source == "reference"
    assert result.unit == pytest.approx(5.0)
    with pytest.raises(ValueError):
        analyse(peaks(*rows), ester(), [],
                reference={"ppm": 6.0, "protons": 1})


def test_nothing_matching_the_product_is_reported():
    result = analyse(peaks((5.5, 1.0), (2.9, 3.0)), ester(), [])
    assert result.scale_source == "total"
    assert result.product.status == "absent"
    assert result.warnings
    assert len(result.product.missing) == 5


def test_fit_rejects_a_peak_that_is_too_small():
    group = Group("CH3", 2.0, 3, 0.2)
    fit = fit_species([(0, 2.0, 1.0)], [group], Settings(), noise=0.1,
                      ratio=1.0)
    assert fit is None
