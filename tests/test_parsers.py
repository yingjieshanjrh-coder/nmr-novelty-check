from nmrcheck.parsers import (filename_to_doi, nmrexp_peaks, norm_nucleus,
                              norm_solvent, read_nmrshiftdb2, read_table)
from nmrcheck.sources import source_path


def test_doi_from_filename():
    assert filename_to_doi("10.1021_acs.orglett.3c00718") == \
        "10.1021/acs.orglett.3c00718"
    assert filename_to_doi("supporting_info") == ""


def test_nucleus_names():
    assert norm_nucleus("13C NMR") == "13C"
    assert norm_nucleus("P NMR") == "31P"
    assert norm_nucleus("{1H}31P NMR") == "31P"
    assert norm_nucleus("2D NMR") is None


def test_solvent_names():
    assert norm_solvent("Chloroform-D1 (CDCl3)") == "CDCl3"
    assert norm_solvent("Acetone-D6 ((CD3)2CO)") == "CD3COCD3"
    assert norm_solvent("Unreported") == ""
    assert norm_solvent("CD2Cl2") == "CD2Cl2"


def test_nmrexp_carbon_peaks():
    text = "[(180.0, 'd', 300.8), (133.7, None, None)]"
    assert nmrexp_peaks("13C", text) == [180.0, 133.7]


def test_nmrexp_proton_peaks():
    text = "[('m', [], '2H', 7.27, 7.24), ('s', [], '9H', 1.33, 1.33)]"
    peaks = nmrexp_peaks("1H", text)
    assert peaks[0] == [7.255, 7.24, 7.27, "m", 2]
    assert peaks[1] == [1.33, 1.33, 1.33, "s", 9]


def test_nmrexp_garbage_is_ignored():
    assert nmrexp_peaks("13C", "not a list") == []
    assert nmrexp_peaks("1H", "[('m', [], None, None, None)]") == []


def test_nmrshiftdb2_reader(data_dir):
    recs = list(read_nmrshiftdb2(source_path(data_dir, "nmrshiftdb2")))
    benzene = [r for r in recs if r["record_id"].startswith("1/")]
    assert {r["nucleus"] for r in benzene} == {"13C", "1H"}
    carbon = next(r for r in benzene if r["nucleus"] == "13C")
    assert carbon["peaks"] == [128.5]
    assert carbon["solvent"] == "CDCl3"
    assert carbon["frequency"] == "400 MHz"
    assert len(carbon["assigned"]) == 6
    proton = next(r for r in benzene if r["nucleus"] == "1H")
    assert proton["peaks"] == [[7.36, 7.36, 7.36, "", 6]]
    acetone = next(r for r in recs if r["record_id"].startswith("7/"))
    assert acetone["origin"] == "calculated"


def test_generic_table(tmp_path):
    path = tmp_path / "extra.csv"
    path.write_text(
        "IUPAC,SMILES,13C NMR chemical shifts,DOI\n"
        'anisole,COc1ccccc1,"159.6, 129.5 (2C), 120.7, 114.0, 55.1",10.1/x\n'
        "empty,CC,,10.1/y\n")
    recs = list(read_table(path))
    assert len(recs) == 1
    assert recs[0]["peaks"] == [159.6, 129.5, 120.7, 114.0, 55.1]
    assert recs[0]["doi"] == "10.1/x"
