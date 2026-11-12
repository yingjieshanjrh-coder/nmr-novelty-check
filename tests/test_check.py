from nmrcheck.check import Log, run
from nmrcheck.report import to_html, to_json, to_markdown


def codes(report):
    return {e["code"] for e in report["degradations"]}


def test_index_contents(index):
    meta = index.meta
    assert meta["sources"]["nmrshiftdb2"]["spectra"] == 10
    assert meta["sources"]["nmrexp"]["spectra"] == 5
    assert meta["missing_sources"] == ["nmrbank"]
    assert len(index.compounds(inchikey="UHOVQNZJYSORNB-UHFFFAOYSA-N")) == 1


def test_benzene_is_reported(index):
    report = run(index, "c1ccccc1")
    v = report["verdict"]
    assert v["status"] == "reported"
    assert v["confidence"] == 0.95
    target = report["target"]
    assert target["record_counts"]["13C"]["exact_experimental"] == 2
    assert target["record_counts"]["1H"]["exact_experimental"] == 2
    sources = {e["source"] for e in target["evidence"]["13C"]}
    assert sources == {"nmrshiftdb2", "nmrexp"}
    row = target["assignment"]["13C"][0]
    assert row["atoms"] == "C1,2,3,4,5,6" and row["observed"] == 128.5
    assert "source_not_indexed" in codes(report)
    assert "prediction_not_independent" in codes(report)


def test_mispaired_record_scores_zero(index):
    report = run(index, "c1ccccc1")
    wrong = [e for e in report["target"]["evidence"]["1H"]
             if e["doi"] == "10.1000/demo.0003"]
    assert wrong and wrong[0]["agreement"]["score"] < 0.05


def test_stereo_variant(index):
    report = run(index, "C[C@@H](N)C(=O)O")
    assert report["verdict"]["status"] == "reported_as_variant"
    assert report["verdict"]["confidence"] == 0.6


def test_calculated_only(index):
    report = run(index, "CC(C)=O")
    assert report["verdict"]["status"] == "listed_without_experiment"


def test_unknown_compound(index):
    report = run(index, "CCCc1ccccc1")
    v = report["verdict"]
    assert v["status"] == "not_found"
    assert v["confidence"] == 0.5
    assert "nearest_known_isomers" in report["target"]
    assert "assignment_predicted_only" in codes(report)
    assert all(r["observed"] is None for r in report["target"]["assignment"]["13C"])


def test_nearest_isomer(index):
    report = run(index, "Cc1ccccc1C")       # o-xylene, isomer of ethylbenzene
    near = report["target"]["nearest_known_isomers"]
    assert [n["smiles"] for n in near] == ["CCc1ccccc1"]


def test_invalid_target(index):
    report = run(index, "this is not smiles")
    assert report["verdict"]["status"] == "invalid_input"
    assert "invalid_structure" in codes(report)
    to_markdown(report)


def test_reaction_input(index):
    report = run(index, "c1ccccc1.CCl>[Al](Cl)(Cl)Cl>Cc1ccccc1",
                 impurities=["CCc1ccccc1", "bad smiles"], solvent="CDCl3")
    assert report["query"]["reaction"]
    assert report["target"]["smiles"] == "Cc1ccccc1"
    roles = [(c["role"], c.get("smiles", c.get("input")))
             for c in report["components"]]
    assert ("reactant", "c1ccccc1") in roles
    assert [r for r, _ in roles].count("reagent") == 1
    assert ("impurity", "bad smiles") in roles
    benzene = next(m for m in report["mixture"]
                   if m["species"] == "c1ccccc1" and m["nucleus"] == "13C")
    assert (128.2, 128.5) in [tuple(x) for x in benzene["overlapping"]]
    solvent = [m for m in report["mixture"] if m["role"] == "solvent"]
    assert len(solvent) == 2


def test_unknown_solvent_is_logged(index):
    report = run(index, "c1ccccc1", solvent="unobtainium-d4")
    assert "solvent_unknown" in codes(report)


def test_renderers(index):
    report = run(index, "c1ccccc1.CCl>>Cc1ccccc1", solvent="CDCl3")
    md = to_markdown(report)
    assert "## Verdict" in md and "## Degradation log" in md
    assert "## Limitations" in md and "## Data sources" in md
    page = to_html(report)
    assert "<table>" in page and "<h2>Verdict</h2>" in page
    assert '"status": "reported"' in to_json(report)


def test_log_deduplicates():
    log = Log()
    log("a", "b", "c")
    log("a", "b", "c")
    assert len(log.entries) == 1
