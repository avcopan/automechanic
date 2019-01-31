""" core functions defining the graph structure

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}
bnd_dct: {bnd_key: (bnd_ord, bnd_ste_par), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from numbers import Integral as _Integer
from itertools import starmap as _starmap
import numpy
from ._dict import right_update as _right_update
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict.multi import by_key_by_position as _by_key_by_position
from ._dict.multi import set_by_key_by_position as _set_by_key_by_position
from ..atom import SYMBOLS as _ATM_SYMS

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1


# constructors
def from_atoms_and_bonds(atm_dct, bnd_dct):
    """ molecular graph from atoms and bonds
    """
    return (atm_dct, bnd_dct)


def from_data(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ molecular graph (any type) from data
    """
    xgr = empty_graph()
    xgr = add_atoms(xgr, atm_sym_dct, atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
                    atm_ste_par_dct=atm_ste_par_dct)
    xgr = add_bonds(xgr, bnd_keys, bnd_ord_dct, bnd_ste_par_dct)
    return xgr


def add_atoms(xgr, atm_sym_dct, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None):
    """ add atoms to this molecular graph
    """
    atm_imp_hyd_vlc_dct = ({} if atm_imp_hyd_vlc_dct is None else
                           atm_imp_hyd_vlc_dct)
    atm_ste_par_dct = {} if atm_ste_par_dct is None else atm_ste_par_dct

    atm_keys = atm_sym_dct.keys()
    atm_keys = [_atom_key(xgr, atm_key) for atm_key in atm_keys]
    assert set(atm_imp_hyd_vlc_dct.keys()) <= set(atm_keys)
    assert set(atm_ste_par_dct.keys()) <= set(atm_keys)

    atm_vals_lst = tuple(_starmap(
        _atom_values,
        zip(*(_values_by_key(atm_sym_dct, atm_keys),
              _values_by_key(atm_imp_hyd_vlc_dct, atm_keys, fill_val=0),
              _values_by_key(atm_ste_par_dct, atm_keys, fill_val=None)))))

    atm_dct = _right_update(atoms(xgr), dict(zip(atm_keys, atm_vals_lst)))
    return from_atoms_and_bonds(atm_dct, bonds(xgr))


def add_bonds(xgr, bnd_keys, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ add bonds to this molecular graph
    """
    bnd_ord_dct = {} if bnd_ord_dct is None else bnd_ord_dct
    bnd_ste_par_dct = {} if bnd_ste_par_dct is None else bnd_ste_par_dct

    bnd_keys = [_bond_key(xgr, *bnd_key) for bnd_key in bnd_keys]
    assert set(bnd_ord_dct.keys()) <= set(bnd_keys)
    assert set(bnd_ste_par_dct.keys()) <= set(bnd_keys)

    bnd_vals_lst = tuple(_starmap(
        _bond_values,
        zip(*(_values_by_key(bnd_ord_dct, bnd_keys, fill_val=1),
              _values_by_key(bnd_ste_par_dct, bnd_keys, fill_val=None)))))

    bnd_dct = _right_update(bonds(xgr), dict(zip(bnd_keys, bnd_vals_lst)))
    return from_atoms_and_bonds(atoms(xgr), bnd_dct)


def empty_graph():
    """ a molecular graph with no atoms
    """
    return from_atoms_and_bonds(dict(), dict())


def frozen(xgr):
    """ hashable, sortable, immutable container of graph data
    """
    atm_keys = sorted(atom_keys(xgr))
    bnd_keys = sorted(bond_keys(xgr), key=sorted)

    # make it sortable by replacing Nones with -infinity
    atm_vals = numpy.array(_values_by_key(atoms(xgr), atm_keys))
    bnd_vals = numpy.array(_values_by_key(bonds(xgr), bnd_keys))
    atm_vals[numpy.equal(atm_vals, None)] = -numpy.inf
    bnd_vals[numpy.equal(bnd_vals, None)] = -numpy.inf

    frz_atms = tuple(zip(atm_keys, map(tuple, atm_vals)))
    frz_bnds = tuple(zip(bnd_keys, map(tuple, bnd_vals)))
    return (frz_atms, frz_bnds)


def _atom_key(xgr, atm_key):
    assert isinstance(atm_key, _Integer)
    assert atm_key not in atom_keys(xgr)
    return int(atm_key)


def _atom_values(atm_sym, atm_imp_hyd_vlc=0, atm_ste_par=None):
    assert atm_sym in _ATM_SYMS
    assert isinstance(atm_imp_hyd_vlc, _Integer)
    assert atm_ste_par in (None, False, True)
    return (atm_sym, atm_imp_hyd_vlc, atm_ste_par)


def _bond_key(xgr, atm1_key, atm2_key):
    assert atm1_key in atom_keys(xgr)
    assert atm2_key in atom_keys(xgr)
    return frozenset({int(atm1_key), int(atm2_key)})


def _bond_values(bnd_ord=1, bnd_ste_par=None):
    assert isinstance(bnd_ord, _Integer)
    assert bnd_ste_par in (None, False, True)
    return (bnd_ord, bnd_ste_par)


# value getters
def atoms(xgr):
    """ atoms, as a dictionary
    """
    atm_dct, _ = xgr
    return atm_dct


def bonds(xgr):
    """ bonds, as a dictionary
    """
    _, bnd_dct = xgr
    return bnd_dct


def atom_keys(xgr):
    """ sorted atom keys
    """
    return frozenset(atoms(xgr).keys())


def bond_keys(xgr):
    """ sorted bond keys
    """
    return frozenset(bonds(xgr).keys())


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_SYM_POS)


def atom_implicit_hydrogen_valences(xgr):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr),
                               ATM_IMP_HYD_VLC_POS)


def atom_stereo_parities(sgr):
    """ atom parities, as a dictionary
    """
    return _by_key_by_position(atoms(sgr), atom_keys(sgr), ATM_STE_PAR_POS)


def bond_orders(rgr):
    """ bond orders, as a dictionary
    """
    return _by_key_by_position(bonds(rgr), bond_keys(rgr), BND_ORD_POS)


def bond_stereo_parities(sgr):
    """ bond parities, as a dictionary
    """
    return _by_key_by_position(bonds(sgr), bond_keys(sgr), BND_STE_PAR_POS)


# value setters
def set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = _set_by_key_by_position(atoms(xgr), atm_imp_hyd_vlc_dct,
                                      ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(xgr)
    return from_atoms_and_bonds(atm_dct, bnd_dct)


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = _set_by_key_by_position(atoms(sgr), atm_par_dct, ATM_STE_PAR_POS)
    return from_atoms_and_bonds(atm_dct, bonds(sgr))


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = _set_by_key_by_position(bonds(rgr), bnd_ord_dct, BND_ORD_POS)
    return from_atoms_and_bonds(atoms(rgr), bnd_dct)


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = _set_by_key_by_position(bonds(sgr), bnd_par_dct, BND_STE_PAR_POS)
    return from_atoms_and_bonds(atoms(sgr), bnd_dct)


# transformations
def without_bond_orders(xgr):
    """ resonance graph with maximum spin (i.e. no pi bonds)
    """
    bnd_ord_dct = _by_key({}, bond_keys(xgr), fill_val=1)
    return set_bond_orders(xgr, bnd_ord_dct)


def without_stereo_parities(xgr):
    """ graph with stereo assignments wiped out
    """
    atm_ste_par_dct = _by_key({}, atom_keys(xgr), fill_val=None)
    bnd_ste_par_dct = _by_key({}, bond_keys(xgr), fill_val=None)
    xgr = set_atom_stereo_parities(xgr, atm_ste_par_dct)
    xgr = set_bond_stereo_parities(xgr, bnd_ste_par_dct)
    return xgr
