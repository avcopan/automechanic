""" functions operating on cartesian geometries (angstroms)
"""
import itertools
import numpy
from .graph import (stereo_inchi_from_coordinates as
                    _stereo_inchi_from_coordinates)


def stereo_inchi(geo):
    """ InChI string of a cartesian geometry
    """
    cgr, atm_xyz_dct = _connectivity_graph_and_atom_coordinates(geo)
    ich = _stereo_inchi_from_coordinates(cgr, atm_xyz_dct)
    return ich


def connectivity_graph(geo):
    """ connectivity graph from a cartesian geometry
    """
    cgr, _ = _connectivity_graph_and_atom_coordinates(geo)
    return cgr


def _connectivity_graph_and_atom_coordinates(geo):
    # using the same cut-offs as x2z:
    xy_bond_max = 3.5 / 1.8897259886
    xh_bond_max = 2.5 / 1.8897259886

    syms, xyzs = zip(*geo)

    atm_sym_dct = dict(enumerate(syms))
    atm_xyz_dct = dict(enumerate(xyzs))

    def _are_bonded(atm_key1, atm_key2):
        atm_sym1 = atm_sym_dct[atm_key1]
        atm_sym2 = atm_sym_dct[atm_key2]
        atm_xyz1 = atm_xyz_dct[atm_key1]
        atm_xyz2 = atm_xyz_dct[atm_key2]
        dist = numpy.linalg.norm(numpy.subtract(atm_xyz1, atm_xyz2))

        ret = False
        if 'H' in (atm_sym1, atm_sym2):
            ret = (dist < xh_bond_max)
        else:
            ret = (dist < xy_bond_max)

        return ret

    atm_dct = dict(enumerate((sym, 0, None) for sym in syms))
    bnd_dct = {frozenset({atm_key1, atm_key2}): (1, None)
               for (atm_key1, atm_key2)
               in itertools.combinations(atm_dct.keys(), r=2)
               if _are_bonded(atm_key1, atm_key2)}
    cgr = (atm_dct, bnd_dct)
    return cgr, atm_xyz_dct
