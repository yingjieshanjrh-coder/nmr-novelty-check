"""Readers that turn each raw source into a common record stream.

Every reader yields dicts with the keys
    smiles / mol     structure (one of the two)
    record_id        identifier inside the source
    nucleus          '13C', '1H', ...
    solvent, frequency, doi
    origin           'experimental' or 'calculated'
    peaks            13C-like: [ppm, ...]
                     1H: [[centre, low, high, multiplicity, n_h], ...]
    assigned         optional [(atom_index, ppm), ...] on the heavy-atom graph
"""
import ast
import csv
import re

from rdkit import Chem

NUCLEUS_ALIASES = {
    "13C": "13C", "C": "13C", "1H": "1H", "H": "1H", "19F": "19F", "F": "19F",
    "31P": "31P", "P": "31P", "{1H}31P": "31P", "11B": "11B", "B": "11B",
    "10B": "10B", "29SI": "29Si", "SI": "29Si", "15N": "15N", "17O": "17O",
}

SOLVENT_ALIASES = {
    "(cd3)2so": "DMSO-d6", "dmso-d6, c2d6so": "DMSO-d6", "dmso": "DMSO-d6",
    "(cd3)2co": "CD3COCD3", "acetone-d6": "CD3COCD3",
    "meod": "CD3OD", "methanol-d4": "CD3OD",
    "unreported": "", "unknown": "", "not_known": "",
}
_TRAILING_FORMULA = re.compile(r"\(((?:[^()]|\([^()]*\))*)\)\s*$")


def norm_nucleus(text):
    key = text.replace("NMR", "").strip().upper()
    return NUCLEUS_ALIASES.get(key)


def norm_solvent(text):
    """Reduce 'Chloroform-D1 (CDCl3)' style names to the short formula."""
    text = (text or "").strip()
    if text.lower() in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[text.lower()]
    m = _TRAILING_FORMULA.search(text)
    if m and m.start() > 0:
        text = m.group(1).strip()
    return SOLVENT_ALIASES.get(text.lower(), text)


def filename_to_doi(name):
    """NMRexp stores DOIs as file stems with the slash replaced by '_'."""
    if not name or not name.startswith("10."):
        return ""
    prefix, sep, rest = name.partition("_")
    return f"{prefix}/{rest}" if sep else name


def _num(value):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if f == f else None


def nmrexp_peaks(nucleus, processed):
    """Decode the NMR_processed column of NMRexp."""
    try:
        items = ast.literal_eval(processed)
    except (ValueError, SyntaxError, MemoryError, RecursionError):
        return []
    peaks = []
    if nucleus == "1H":
        for item in items:
            if not isinstance(item, tuple) or len(item) < 5:
                continue
            mult, _, count, high, low = item[:5]
            high, low = _num(high), _num(low)
            if high is None and low is None:
                continue
            high = low if high is None else high
            low = high if low is None else low
            low, high = min(low, high), max(low, high)
            m = re.match(r"\s*(\d+)", str(count or ""))
            n_h = int(m.group(1)) if m else 0
            peaks.append([round((low + high) / 2, 3), low, high, mult or "", n_h])
    else:
        for item in items:
            shift = _num(item[0]) if isinstance(item, tuple) and item else _num(item)
            if shift is not None:
                peaks.append(shift)
    return peaks


def read_nmrexp(path, batch_size=50000):
    import pyarrow.parquet as pq

    columns = ["Filename", "SMILES", "NMR_type", "NMR_frequency", "NMR_solvent",
               "NMR_processed"]
    pf = pq.ParquetFile(path)
    row = 0
    for batch in pf.iter_batches(batch_size=batch_size, columns=columns):
        cols = [batch.column(c).to_pylist() for c in columns]
        for fname, smiles, ntype, freq, solvent, processed in zip(*cols):
            row += 1
            nucleus = norm_nucleus(ntype or "")
            if not nucleus or not smiles or not processed:
                continue
            peaks = nmrexp_peaks(nucleus, processed)
            if not peaks:
                continue
            yield {
                "smiles": smiles,
                "record_id": str(row - 1),
                "nucleus": nucleus,
                "solvent": norm_solvent(solvent),
                "frequency": (freq or "").strip(),
                "doi": filename_to_doi(fname),
                "origin": "experimental",
                "peaks": peaks,
            }


def _indexed_field(value):
    """nmrshiftdb2 metadata looks like '0:298 1:300 '; map index -> text."""
    out = {}
    for m in re.finditer(r"(\d+):(.*?)(?=\s\d+:|$)", value or "", flags=re.S):
        out[m.group(1)] = m.group(2).strip()
    return out


def _heavy_index_map(mol):
    """Atom index in the file -> index after hydrogens are removed."""
    mapping, n = {}, 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            mapping[atom.GetIdx()] = n
            n += 1
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            nbrs = atom.GetNeighbors()
            if nbrs and nbrs[0].GetIdx() in mapping:
                mapping[atom.GetIdx()] = mapping[nbrs[0].GetIdx()]
    return mapping


def read_nmrshiftdb2(path):
    supplier = Chem.ForwardSDMolSupplier(str(path), removeHs=False, sanitize=True)
    for raw in supplier:
        if raw is None:
            continue
        props = raw.GetPropsAsDict(includePrivate=False, includeComputed=False)
        mol_id = str(props.get("nmrshiftdb2 ID", "")).strip()
        try:
            index_map = _heavy_index_map(raw)
            mol = Chem.RemoveHs(raw)
        except Exception:
            continue
        solvents = _indexed_field(str(props.get("Solvent", "")))
        fields = _indexed_field(str(props.get("Field Strength [MHz]", "")))
        computed = {k for k, v in
                    _indexed_field(str(props.get("NMRProgram", ""))).items() if v}
        for name, value in props.items():
            m = re.match(r"Spectrum (\S+) (\d+)$", name)
            if not m:
                continue
            nucleus = norm_nucleus(m.group(1))
            if nucleus is None:
                continue
            spec_no = m.group(2)
            assigned = []
            for chunk in str(value).split("|"):
                bits = chunk.split(";")
                if len(bits) < 3:
                    continue
                shift = _num(bits[0])
                if shift is None:
                    continue
                try:
                    atom = index_map.get(int(bits[2]))
                except ValueError:
                    atom = None
                assigned.append((atom, shift))
            if not assigned:
                continue
            if nucleus == "1H":
                groups = {}
                for atom, shift in assigned:
                    groups.setdefault((atom, shift), 0)
                    groups[(atom, shift)] += 1
                merged = {}
                for (atom, shift), n in groups.items():
                    merged[shift] = merged.get(shift, 0) + n
                peaks = [[s, s, s, "", n] for s, n in sorted(merged.items(),
                                                              reverse=True)]
            else:
                peaks = sorted({s for _, s in assigned}, reverse=True)
            freq = fields.get(spec_no, "")
            yield {
                "mol": mol,
                "record_id": f"{mol_id}/{spec_no}",
                "nucleus": nucleus,
                "solvent": norm_solvent(solvents.get(spec_no, "")),
                "frequency": f"{freq} MHz" if freq and freq != "Unreported" else "",
                "doi": "",
                "origin": "calculated" if spec_no in computed else "experimental",
                "peaks": peaks,
                "assigned": [(a, s) for a, s in assigned if a is not None],
            }


_SHIFT_NUMBER = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?![\w.])")


def _loose_shift_list(text):
    """Pull ppm values out of free text, ignoring anything in parentheses."""
    bare = re.sub(r"\([^)]*\)", " ", text or "")
    return [float(x) for x in _SHIFT_NUMBER.findall(bare)]


def read_table(path):
    """Generic CSV reader for extra sources such as NMRBank.

    Columns are found by name: one containing 'smiles', one containing both
    '13c' and 'shift', optionally 'doi' and a 13C solvent/condition column.
    Only 13C lists are taken because 1H text mixes shifts with couplings.
    """
    csv.field_size_limit(1 << 30)
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        names = {n.lower(): n for n in reader.fieldnames or []}

        def find(*words):
            for low, orig in names.items():
                if all(w in low for w in words):
                    return orig
            return None

        col_smiles = find("smiles")
        col_c = find("13c", "shift")
        col_doi = find("doi")
        col_cond = find("13c", "condition") or find("solvent")
        if not col_smiles or not col_c:
            raise ValueError(f"{path}: need a SMILES column and a 13C shift column")
        for n, row in enumerate(reader):
            peaks = _loose_shift_list(row.get(col_c))
            if not row.get(col_smiles) or not peaks:
                continue
            yield {
                "smiles": row[col_smiles],
                "record_id": str(n),
                "nucleus": "13C",
                "solvent": norm_solvent(row.get(col_cond, "")) if col_cond else "",
                "frequency": "",
                "doi": (row.get(col_doi) or "").strip() if col_doi else "",
                "origin": "experimental",
                "peaks": peaks,
            }


READERS = {
    "nmrshiftdb2": read_nmrshiftdb2,
    "nmrexp": read_nmrexp,
    "nmrbank": read_table,
}
