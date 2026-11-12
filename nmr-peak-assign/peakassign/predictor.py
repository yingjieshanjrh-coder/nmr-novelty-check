"""Predicted 1H signals for a SMILES, using the parent project's predictor."""
from .assign import Group, Species
from .parent import ensure_importable

# Extra search width on top of the tolerance, by the radius of the atom
# environment the prediction came from. The parent index reports hold-out
# errors of 0.14 ppm (radius 4) to 0.63 ppm (radius 1) for 1H; radius -1
# means the table had nothing and a fixed atom-type value was used.
SLACK = {4: 0.12, 3: 0.20, 2: 0.30, 1: 0.45, 0: 0.60, -1: 0.90}
LABILE = {"O", "N", "S"}


def _kind(atom):
    hydrogens = atom.GetTotalNumHs()
    if atom.GetSymbol() != "C":
        return atom.GetSymbol() + "H" + (str(hydrogens) if hydrogens > 1 else "")
    if atom.GetIsAromatic():
        return "ArH"
    return {3: "CH3", 2: "CH2"}.get(hydrogens, "CH")


def _atom_map(mol):
    """SMILES with 1-based atom numbers, so labels like C8 can be looked up."""
    from rdkit import Chem

    tagged = Chem.Mol(mol)
    for atom in tagged.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return Chem.MolToSmiles(tagged)


def build_species(smiles, key, role, name="", index=None, tol=0.08):
    """Predict the 1H signals of one compound.

    `index` is the parent's Index (or None, in which case every shift falls
    back to a fixed atom-type value and the windows are very wide).
    """
    ensure_importable()
    from nmrcheck.chem import parse_smiles
    from nmrcheck.predict import predict

    mol = parse_smiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse SMILES: {smiles}")
    lookup = index.lookup_env if index is not None else (lambda *args: None)

    groups = []
    for signal in predict(mol, "1H", lookup):
        atom = mol.GetAtomWithIdx(signal["atoms"][0])
        numbers = ",".join(str(a + 1) for a in signal["atoms"])
        window = tol + SLACK.get(signal["radius"], SLACK[-1])
        groups.append(Group(
            label=f"{_kind(atom)}（{atom.GetSymbol()}{numbers}）",
            shift=signal["shift"], protons=signal["count"],
            window=round(window, 3), tight=round(max(tol, window / 2), 3),
            radius=signal["radius"], labile=atom.GetSymbol() in LABILE,
            atoms=tuple(signal["atoms"])))
    return Species(key=key, role=role, name=name or smiles, groups=groups,
                   smiles=smiles, atom_map=_atom_map(mol))
