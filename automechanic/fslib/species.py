""" functions for generating the species filesystem
"""
from .. import mol
from .. import params as par


def assert_complete_species_id(ich, mult):
    """ is this species addressable, given the available information?
    """
    ick = mol.inchi.inchi_key(ich)
    cgr = mol.inchi.connectivity_graph(ich)
    assert mol.inchi.key.is_standard_neutral(ick)
    assert mol.inchi.has_unknown_stereo_elements(ich) is False
    assert mol.inchi.is_closed(ich)
    assert mult in mol.graph.possible_spin_multiplicities(cgr)


def branch_segments(ich, mult):
    """ get the species address from its InChI string and multiplicity
    """
    assert_complete_species_id(ich, mult)
    sgms = (_base_segment(),
            _formula_segment(ich),
            _connectivity_segment(ich),
            _multiplicity_segment(mult),
            _stereochemistry_segment(ich))
    return sgms


def _base_segment():
    dir_name = par.SPC.FILESYSTEM_DIR_NAME
    info = None
    return (dir_name, info)


def _formula_segment(ich):
    dir_name = mol.inchi.formula_layer(ich)
    info = None
    return (dir_name, info)


def _connectivity_segment(ich):
    ick = mol.inchi.inchi_key(ich)
    dir_name = mol.inchi.key.first_hash(ick)
    inf_dct = _inchi_hash_information(mol.inchi.core_parent(ich))
    return (dir_name, inf_dct)


def _multiplicity_segment(mult):
    dir_name = '{:d}'.format(mult)
    inf_dct = None
    return (dir_name, inf_dct)


def _stereochemistry_segment(ich):
    ick = mol.inchi.inchi_key(ich)
    dir_name = mol.inchi.key.second_hash(ick)
    inf_dct = _inchi_hash_information(ich)
    return (dir_name, inf_dct)


def _inchi_hash_information(ich):
    smi = mol.inchi.smiles(ich)
    inf_dct = {par.SPC.ID_ICH_KEY: ich,
               par.SPC.ID_SMI_KEY: smi}
    return inf_dct
