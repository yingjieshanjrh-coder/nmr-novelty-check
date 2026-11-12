"""Readers for targets.yaml and measured_peaks.csv."""
import csv
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import impurities


class InputError(ValueError):
    pass


@dataclass
class Compound:
    smiles: str
    name: str = ""


@dataclass
class Targets:
    product: Compound
    reactants: list
    solvent: str
    extra_solvents: list = field(default_factory=list)
    extra_peaks: list = field(default_factory=list)
    reference: dict = None


@dataclass
class Peak:
    ppm: float
    integral: float
    multiplicity: str = ""


def _compound(entry, what):
    if isinstance(entry, str) and entry.strip():
        return Compound(entry.strip())
    if isinstance(entry, dict) and str(entry.get("smiles") or "").strip():
        return Compound(str(entry["smiles"]).strip(),
                        str(entry.get("name") or "").strip())
    raise InputError(f"{what}: expected a SMILES string or "
                     "{smiles: ..., name: ...}")


def load_targets(path):
    try:
        doc = yaml.safe_load(Path(path).read_text(encoding="utf-8-sig"))
    except yaml.YAMLError as exc:
        raise InputError(f"{path}: not valid YAML ({exc})") from None
    if not isinstance(doc, dict):
        raise InputError(f"{path}: expected keys product, reactants, solvent")
    if "product" not in doc:
        raise InputError(f"{path}: 'product' is required")
    if not isinstance(doc.get("reactants"), list):
        raise InputError(f"{path}: 'reactants' must be a list "
                         "(use [] if there are none)")
    solvent = impurities.resolve_solvent(str(doc.get("solvent") or ""))
    if solvent is None:
        raise InputError(
            f"{path}: unknown solvent {doc.get('solvent')!r}; known: "
            + ", ".join(impurities.known_solvents()))

    targets = Targets(
        product=_compound(doc["product"], "product"),
        reactants=[_compound(r, f"reactants[{i}]")
                   for i, r in enumerate(doc["reactants"])],
        solvent=solvent)

    for name in doc.get("extra_solvents") or []:
        found = impurities.resolve_trace(str(name), solvent)
        if found is None:
            raise InputError(
                f"{path}: no tabulated shifts for {name!r} in {solvent}; "
                "list its signals under extra_peaks instead")
        if found not in targets.extra_solvents:
            targets.extra_solvents.append(found)

    for i, item in enumerate(doc.get("extra_peaks") or []):
        try:
            entry = {"name": str(item["name"]), "ppm": float(item["ppm"])}
        except (TypeError, KeyError, ValueError):
            raise InputError(f"{path}: extra_peaks[{i}] needs name and ppm") \
                from None
        targets.extra_peaks.append(entry)

    ref = doc.get("reference")
    if ref is not None:
        try:
            targets.reference = {"ppm": float(ref["ppm"]),
                                 "protons": float(ref["protons"])}
        except (TypeError, KeyError, ValueError):
            raise InputError(f"{path}: reference needs ppm and protons") \
                from None
        if targets.reference["protons"] <= 0:
            raise InputError(f"{path}: reference protons must be positive")
    return targets


PPM_COLUMNS = ("ppm", "shift", "delta")
INTEGRAL_COLUMNS = ("integral", "area", "integration")
MULT_COLUMNS = ("multiplicity", "mult")


def _column(header, names):
    for name in names:
        if name in header:
            return header.index(name)
    return None


def load_peaks(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = [r for r in csv.reader(handle) if any(c.strip() for c in r)]
    if not rows:
        raise InputError(f"{path}: file is empty")
    header = [c.strip().lower() for c in rows[0]]
    ppm_col = _column(header, PPM_COLUMNS)
    int_col = _column(header, INTEGRAL_COLUMNS)
    mult_col = _column(header, MULT_COLUMNS)
    if ppm_col is None or int_col is None:
        raise InputError(f"{path}: header must contain 'ppm' and 'integral'")

    peaks = []
    for n, row in enumerate(rows[1:], start=2):
        try:
            ppm = float(row[ppm_col])
            integral = float(row[int_col])
        except (ValueError, IndexError):
            raise InputError(f"{path}, row {n}: ppm and integral must be "
                             "numbers") from None
        if integral <= 0:
            raise InputError(f"{path}, row {n}: integral must be positive")
        mult = ""
        if mult_col is not None and mult_col < len(row):
            mult = row[mult_col].strip()
        peaks.append(Peak(ppm, integral, mult))
    if not peaks:
        raise InputError(f"{path}: no peaks found")
    peaks.sort(key=lambda p: -p.ppm)
    return peaks
