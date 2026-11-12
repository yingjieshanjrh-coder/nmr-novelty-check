"""Chemical shift prediction from atom environments.

The predictor is a lookup table trained on the assigned spectra in
nmrshiftdb2. Each atom is described by circular environments of radius
0..MAX_RADIUS bonds (Morgan identifiers); the widest environment seen in
training supplies the predicted shift. This is the same idea as HOSE-code
prediction, with RDKit environments instead of HOSE strings.
"""
from statistics import median

from rdkit.Chem import rdFingerprintGenerator

from .chem import symmetry_classes

MAX_RADIUS = 4
TARGET_ELEMENT = {"13C": 6, "1H": None}

_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MAX_RADIUS, includeRedundantEnvironments=True)


def atom_environments(mol):
    """Return {atom_index: {radius: identifier}} for a heavy-atom molecule."""
    extra = rdFingerprintGenerator.AdditionalOutput()
    extra.AllocateBitInfoMap()
    _generator.GetSparseCountFingerprint(mol, additionalOutput=extra)
    envs = {}
    for ident, hits in extra.GetBitInfoMap().items():
        for atom, radius in hits:
            # identifiers are unsigned 32 bit; keep them inside sqlite's range
            envs.setdefault(atom, {})[radius] = int(ident)
    return envs


def nucleus_atoms(mol, nucleus):
    """Heavy atoms that carry the observed nucleus, with their multiplicity."""
    out = {}
    for atom in mol.GetAtoms():
        if nucleus == "13C" and atom.GetAtomicNum() == 6:
            out[atom.GetIdx()] = 1
        elif nucleus == "1H" and atom.GetTotalNumHs() > 0:
            out[atom.GetIdx()] = atom.GetTotalNumHs()
    return out


class ShiftTable:
    """Accumulates shifts per environment while the index is being built."""

    def __init__(self):
        self.values = {}

    def add(self, mol, nucleus, assigned):
        wanted = nucleus_atoms(mol, nucleus)
        per_atom = {}
        for atom, shift in assigned:
            if atom in wanted:
                per_atom.setdefault(atom, []).append(shift)
        if not per_atom:
            return
        envs = atom_environments(mol)
        for atom, shifts in per_atom.items():
            shift = sum(shifts) / len(shifts)
            for radius, ident in envs.get(atom, {}).items():
                self.values.setdefault((nucleus, radius, ident), []).append(shift)

    def rows(self):
        for (nucleus, radius, ident), shifts in self.values.items():
            n = len(shifts)
            mid = median(shifts)
            spread = median(abs(s - mid) for s in shifts) if n > 1 else 0.0
            yield nucleus, radius, ident, n, mid, spread

    def lookup(self, nucleus, radius, ident):
        shifts = self.values.get((nucleus, radius, ident))
        if not shifts:
            return None
        return len(shifts), median(shifts)


FALLBACK = {
    # used only when not even the bare atom type is in the table
    "13C": {"aromatic": 128.0, "sp2": 130.0, "sp": 80.0, "sp3": 30.0},
    "1H": {"aromatic": 7.3, "sp2": 5.5, "sp": 2.5, "sp3": 1.5},
}


def _fallback(atom, nucleus):
    from rdkit.Chem import HybridizationType as H
    if atom.GetIsAromatic():
        kind = "aromatic"
    else:
        kind = {H.SP2: "sp2", H.SP: "sp"}.get(atom.GetHybridization(), "sp3")
    return FALLBACK[nucleus][kind]


def predict(mol, nucleus, lookup, min_radius=0):
    """Predict shifts for every atom carrying `nucleus`.

    `lookup(nucleus, radius, ident)` returns (n, shift) or None. The result is
    one entry per symmetry class:
        {"atoms": [...], "shift": ppm, "radius": r, "support": n, "count": k}
    radius -1 marks a table miss that fell back to a fixed atom-type value.
    """
    if nucleus not in TARGET_ELEMENT:
        return []
    wanted = nucleus_atoms(mol, nucleus)
    if not wanted:
        return []
    envs = atom_environments(mol)
    per_atom = {}
    for atom in wanted:
        hit = None
        for radius in range(MAX_RADIUS, min_radius - 1, -1):
            ident = envs.get(atom, {}).get(radius)
            if ident is None:
                continue
            found = lookup(nucleus, radius, ident)
            if found:
                hit = (found[1], radius, found[0])
                break
        if hit is None:
            hit = (_fallback(mol.GetAtomWithIdx(atom), nucleus), -1, 0)
        per_atom[atom] = hit

    out = []
    for group in symmetry_classes(mol):
        atoms = [a for a in group if a in per_atom]
        if not atoms:
            continue
        shifts = [per_atom[a][0] for a in atoms]
        out.append({
            "atoms": atoms,
            "shift": round(sum(shifts) / len(shifts), 2),
            "radius": min(per_atom[a][1] for a in atoms),
            "support": min(per_atom[a][2] for a in atoms),
            "count": sum(wanted[a] for a in atoms),
        })
    out.sort(key=lambda p: -p["shift"])
    return out
