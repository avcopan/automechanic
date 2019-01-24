""" stereo graph -> InChI string conversion
"""
from ._intco import atom_stereo_coordinates
from .._base import explicit_stereo_sites
from ..to_inchi import (with_atom_inchi_numbers as
                        _inchi_with_atom_inchi_numbers)


def inchi(sgr):
    """ InChI string of this stereo graph
    """
    sgr = explicit_stereo_sites(sgr)
    atm_xyz_dct = atom_stereo_coordinates(sgr)
    ich, _ = _inchi_with_atom_inchi_numbers(sgr, atm_xyz_dct=atm_xyz_dct)
    return ich
