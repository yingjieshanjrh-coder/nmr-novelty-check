import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem

from nmrcheck.index import Index, build
from nmrcheck.sources import source_path

# smiles, id, {nucleus: [(atom index, ppm), ...]}
SD_RECORDS = [
    ("c1ccccc1", 1, {"13C": [(i, 128.5) for i in range(6)],
                     "1H": [(i, 7.36) for i in range(6)]}),
    ("Cc1ccccc1", 2, {"13C": [(0, 21.4), (1, 137.8), (2, 129.0), (6, 129.0),
                              (3, 128.2), (5, 128.2), (4, 125.3)],
                      "1H": [(0, 2.34)] * 3 + [(2, 7.17), (6, 7.17), (3, 7.25),
                                               (5, 7.25), (4, 7.17)]}),
    ("CCO", 3, {"13C": [(0, 18.1), (1, 57.8)],
                "1H": [(0, 1.22)] * 3 + [(1, 3.69)] * 2}),
    ("C1CCCCC1", 4, {"13C": [(i, 27.1) for i in range(6)]}),
    ("C[C@H](N)C(=O)O", 5, {"13C": [(0, 16.9), (1, 51.2), (3, 176.5)]}),
    ("CCc1ccccc1", 6, {"13C": [(0, 15.6), (1, 28.9), (2, 144.2), (3, 127.8),
                               (7, 127.8), (4, 128.3), (6, 128.3), (5, 125.6)]}),
]


def sd_block(smiles, mol_id, spectra, calculated=False):
    mol = Chem.MolFromSmiles(smiles)
    lines = [Chem.MolToMolBlock(mol).rstrip("\n")]
    lines += ["> <nmrshiftdb2 ID>", str(mol_id), ""]
    solvents, fields = [], []
    for n, (nucleus, peaks) in enumerate(spectra.items()):
        body = "|".join(f"{ppm};0.0;{atom}" for atom, ppm in peaks) + "|"
        lines += [f"> <Spectrum {nucleus} {n}>", body, ""]
        solvents.append(f"{n}:Chloroform-D1 (CDCl3)")
        fields.append(f"{n}:400")
    lines += ["> <Solvent>", " ".join(solvents) + " ", ""]
    lines += ["> <Field Strength [MHz]>", " ".join(fields) + " ", ""]
    if calculated:
        lines += ["> <NMRProgram>", "0:Gaussian ", ""]
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


NMREXP_ROWS = [
    # a correct pair, a wrongly paired paragraph, an unparsable structure
    ("10.1000_demo.0001", "COc1ccccc1", "13C NMR", "101 MHz", "CDCl3",
     "[(159.6, None, None), (129.5, None, None), (120.7, None, None), "
     "(114.0, None, None), (55.1, None, None)]"),
    ("10.1000_demo.0001", "COc1ccccc1", "1H NMR", "400 MHz", "CDCl3",
     "[('m', [], '2H', 7.32, 7.26), ('m', [], '3H', 6.97, 6.89), "
     "('s', [], '3H', 3.81, 3.81)]"),
    ("10.1000_demo.0002", "c1ccccc1", "13C NMR", "101 MHz", "CDCl3",
     "[(128.4, None, None)]"),
    ("10.1000_demo.0003", "c1ccccc1", "1H NMR", "400 MHz", "CDCl3",
     "[('t', ['7.0Hz'], '3H', 0.82, 0.82), ('d', ['6.5Hz'], '3H', 0.97, 0.97)]"),
    ("10.1000_demo.0004", "not-a-smiles", "13C NMR", "101 MHz", "CDCl3",
     "[(10.0, None, None)]"),
    ("10.1000_demo.0005", "FC(F)(F)c1ccccc1", "19F NMR", "376 MHz", "CDCl3",
     "[(-62.8, 's', None)]"),
]


def write_sources(data_dir, extra_sd=""):
    sd = source_path(data_dir, "nmrshiftdb2")
    sd.parent.mkdir(parents=True, exist_ok=True)
    sd.write_text("".join(sd_block(*r) for r in SD_RECORDS) + extra_sd)

    cols = ["Filename", "SMILES", "NMR_type", "NMR_frequency", "NMR_solvent",
            "NMR_processed"]
    table = pa.table({c: [row[i] for row in NMREXP_ROWS]
                      for i, c in enumerate(cols)})
    pqf = source_path(data_dir, "nmrexp")
    pqf.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, pqf)
    return {"nmrshiftdb2": sd, "nmrexp": pqf}


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    path = tmp_path_factory.mktemp("data")
    calc = sd_block("CC(C)=O", 7, {"13C": [(0, 30.0), (2, 30.0), (1, 205.0)]},
                    calculated=True)
    files = write_sources(path, extra_sd=calc)
    build(path, files, workers=1, quiet=True)
    return path


@pytest.fixture(scope="session")
def index(data_dir):
    idx = Index(data_dir)
    yield idx
    idx.close()
