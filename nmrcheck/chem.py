"""Structure normalisation helpers built on RDKit."""
from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class Identity:
    inchikey: str
    smiles: str
    formula: str

    @property
    def skeleton(self):
        """First InChIKey block: connectivity only, no stereo or isotopes."""
        return self.inchikey[:14]


def parse_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return mol


def identify(mol):
    key = Chem.MolToInchiKey(mol)
    if not key:
        return None
    return Identity(
        inchikey=key,
        smiles=Chem.MolToSmiles(mol),
        formula=rdMolDescriptors.CalcMolFormula(mol),
    )


def identify_smiles(smiles):
    """Return (inchikey, canonical smiles, formula) or None. Picklable for pools."""
    try:
        mol = parse_smiles(smiles)
        if mol is None:
            return None
        ident = identify(mol)
    except Exception:
        return None
    if ident is None:
        return None
    return ident.inchikey, ident.smiles, ident.formula


def symmetry_classes(mol):
    """Group atom indices that are topologically equivalent."""
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    groups = {}
    for idx, rank in enumerate(ranks):
        groups.setdefault(rank, []).append(idx)
    return list(groups.values())


def parse_reaction(text):
    """Split 'reactants>agents>products' into three lists of SMILES."""
    parts = text.split(">")
    if len(parts) != 3:
        raise ValueError("reaction SMILES must look like reactants>agents>products")
    return [[s for s in p.split(".") if s] for p in parts]
