""" specific connnectivity graph functions
"""
from .to_inchi import (with_atom_inchi_numbers as
                       _inchi_with_atom_inchi_numbers)


def inchi(xgr):
    """ InChI string of this connectivity graph
    """
    ich, _ = _inchi_with_atom_inchi_numbers(xgr)
    return ich


def atom_inchi_numbers(xgr):
    """ InChI numbers, by atom
    """
    _, atm_ich_num_dct = _inchi_with_atom_inchi_numbers(xgr)
    return atm_ich_num_dct
