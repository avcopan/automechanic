""" rdkit interface
"""
from rdkit import RDLogger
import rdkit.Chem as _rd_chem
import rdkit.Chem.AllChem as _rd_all_chem

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


def inchi_to_inchi_key(ich):
    """ InChI-Key from an InChI string
    """
    ick = _rd_chem.inchi.InchiToInchiKey(ich)
    return ick


def from_smiles(smi):
    """ rdkit molecule object from a SMILES string
    """
    rdm = _rd_chem.MolFromSmiles(smi)
    assert rdm is not None
    return rdm


def from_inchi(ich):
    """ rdkit molecule object from an InChI string
    """
    rdm = _rd_chem.inchi.MolFromInchi(ich, treatWarningAsError=False)
    assert rdm is not None
    return rdm


def to_smiles(rdm):
    """ SMILES string from an rdkit molecule object
    """
    smi = _rd_chem.MolToSmiles(rdm)
    return smi


def to_inchi(rdm, options='', with_aux_info=False):
    """ InChI string from an rdkit molecule object
    """
    if with_aux_info:
        ret = _rd_chem.inchi.MolToInchiAndAuxInfo(rdm, options=options)
    else:
        ret = _rd_chem.inchi.MolToInchi(rdm, options=options)
    return ret


def geometry(rdm):
    """ cartesian geometry from an rdkit molecule object
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    natms = len(rdm.GetAtoms())
    if natms == 1:
        asb = atms[0].GetSymbol()
        xyz = (0., 0., 0.)
        geo = ((asb, xyz),)
    else:
        _rd_all_chem.EmbedMolecule(rdm)
        _rd_all_chem.MMFFOptimizeMolecule(rdm)
        asbs = tuple(rda.GetSymbol() for rda in atms)
        xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
        geo = tuple(zip(asbs, xyzs))
    return geo


def connectivity_graph(rdm):
    """ connection graph from an rdkit molecule object
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    bnds = rdm.GetBonds()
    asbs = dict(enumerate((rda.GetSymbol(), 0, None) for rda in atms))
    cnns = {frozenset([rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx()]): (1, None)
            for rdb in bnds}
    return (asbs, cnns)
