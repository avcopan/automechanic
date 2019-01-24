""" base graph library; depending only on connectivity

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, ...), ...}
bnd_dct: {bnd_key: (1 or bnd_ord, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from numbers import Integral as _Integer
from itertools import chain as _chain
from itertools import starmap as _starmap
import numpy
from ._seq import remove as _remove
from ._seq import filter_ as _filter
from ._seq import filterfalse as _filterfalse
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict import keys_by_value as _keys_by_value
from ._dict import transform_keys as _transform_keys
from ._dict import transform_values as _transform_values
from ._dict import filter_by_value as _filter_by_value
from ._tdict import by_key_by_position as _by_key_by_position
from ._tdict import set_by_key_by_position as _set_by_key_by_position
from ._networkx import from_graph as _nxg_from_graph
from ._networkx import ring_keys_list as _nxg_ring_keys_list
from ._networkx import isomorphism as _nxg_isomorphism
from ..atom import nuclear_charge as _atom_nuclear_charge
from ..atom import valence as _atom_valence
from ..atom import SYMBOLS as _ATM_SYMS

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1
ATM_STE_PAR_POS = 2

BND_ORD_POS = 0
BND_STE_PAR_POS = 1


# constructors
def from_data(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ molecular graph (any type) from data
    """
    xgr = empty_graph()
    xgr = add_atoms(xgr, atm_sym_dct, atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
                    atm_ste_par_dct=atm_ste_par_dct)
    xgr = add_bonds(xgr, bnd_keys, bnd_ord_dct, bnd_ste_par_dct)
    return xgr


def empty_graph():
    """ a molecular graph with no atoms
    """
    return (dict(), dict())


def add_atoms(xgr, atm_sym_dct, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None):
    """ add atoms to this molecular graph
    """
    atm_imp_hyd_vlc_dct = ({} if atm_imp_hyd_vlc_dct is None else
                           atm_imp_hyd_vlc_dct)
    atm_ste_par_dct = {} if atm_ste_par_dct is None else atm_ste_par_dct

    assert hasattr(atm_sym_dct, 'keys')
    assert hasattr(atm_imp_hyd_vlc_dct, 'keys')
    assert hasattr(atm_ste_par_dct, 'keys')

    atm_keys = atm_sym_dct.keys()
    atm_keys = [_atom_key(xgr, atm_key) for atm_key in atm_keys]
    assert set(atm_imp_hyd_vlc_dct.keys()) <= set(atm_keys)
    assert set(atm_ste_par_dct.keys()) <= set(atm_keys)

    atm_vals_lst = tuple(_starmap(
        _atom_values,
        zip(*(
            _values_by_key(atm_sym_dct, atm_keys),
            _values_by_key(atm_imp_hyd_vlc_dct, atm_keys, fill_val=0),
            _values_by_key(atm_ste_par_dct, atm_keys, fill_val=None)))))

    atm_dct = dict.copy(atoms(xgr))
    atm_dct.update(dict(zip(atm_keys, atm_vals_lst)))
    xgr = (atm_dct, bonds(xgr))
    return xgr


def add_bonds(xgr, bnd_keys, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ add bonds to this molecular graph
    """
    bnd_ord_dct = {} if bnd_ord_dct is None else bnd_ord_dct
    bnd_ste_par_dct = {} if bnd_ste_par_dct is None else bnd_ste_par_dct

    assert hasattr(bnd_ord_dct, 'keys')
    assert hasattr(bnd_ste_par_dct, 'keys')

    bnd_keys = [_bond_key(xgr, *bnd_key) for bnd_key in bnd_keys]
    assert set(bnd_ord_dct.keys()) <= set(bnd_keys)
    assert set(bnd_ste_par_dct.keys()) <= set(bnd_keys)

    bnd_vals_lst = tuple(_starmap(
        _bond_values,
        zip(*(_values_by_key(bnd_ord_dct, bnd_keys, fill_val=1),
              _values_by_key(bnd_ste_par_dct, bnd_keys, fill_val=None)))))
    bnd_dct = dict.copy(bonds(xgr))
    bnd_dct.update(dict(zip(bnd_keys, bnd_vals_lst)))
    xgr = (atoms(xgr), bnd_dct)
    return xgr


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
    return tuple(sorted(atoms(xgr).keys()))


def bond_keys(xgr):
    """ sorted bond keys
    """
    return tuple(sorted(bonds(xgr).keys(), key=sorted))


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_SYM_POS)


def atom_implicit_hydrogen_valences(xgr):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr),
                               ATM_IMP_HYD_VLC_POS)


def atom_stereo_keys(sgr):
    """ keys to atom stereo-centers
    """
    atm_ste_keys = _keys_by_value(atom_stereo_parities(sgr), [True, False])
    return atm_ste_keys


def atom_stereo_parities(sgr):
    """ atom parities, as a dictionary
    """
    return _by_key_by_position(atoms(sgr), atom_keys(sgr), ATM_STE_PAR_POS)


def bond_orders(rgr):
    """ bond orders, as a dictionary
    """
    return _by_key_by_position(bonds(rgr), bond_keys(rgr), BND_ORD_POS)


def bond_stereo_keys(sgr):
    """ keys to bond stereo-centers
    """
    bnd_ste_keys = _keys_by_value(bond_stereo_parities(sgr), [True, False])
    return bnd_ste_keys


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
    xgr = (atm_dct, bnd_dct)
    return xgr


def set_atom_stereo_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = _set_by_key_by_position(atoms(sgr), atm_par_dct, ATM_STE_PAR_POS)
    sgr = (atm_dct, bonds(sgr))
    return sgr


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    bnd_dct = _set_by_key_by_position(bonds(rgr), bnd_ord_dct, BND_ORD_POS)
    rgr = (atoms(rgr), bnd_dct)
    return rgr


def set_bond_stereo_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    bnd_dct = _set_by_key_by_position(bonds(sgr), bnd_par_dct, BND_STE_PAR_POS)
    sgr = (atoms(sgr), bnd_dct)
    return sgr


# derived values
def is_chiral(sgr):
    """ is this stereo graph chiral?
    """
    return not backbone_isomorphic(sgr, reflection(sgr))


def ring_keys_list(xgr):
    """ a series of key-sets for each ring in the graph
    """
    nxg = _nxg_from_graph(xgr)
    rng_keys_lst = _nxg_ring_keys_list(nxg)
    return rng_keys_lst


def backbone_keys(xgr):
    """ backbone atom keys
    """
    bbn_keys = _remove(atom_keys(xgr), explicit_hydrogen_keys(xgr))
    return bbn_keys


def explicit_hydrogen_keys(xgr):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = _keys_by_value(atom_symbols(xgr), ('H',))
    atm_ngb_keys_dct = atom_neighbor_keys(xgr)

    def _is_backbone(hyd_key):
        return all(ngb_key in hyd_keys and hyd_key < ngb_key
                   for ngb_key in atm_ngb_keys_dct[hyd_key])

    exp_hyd_keys = _filterfalse(_is_backbone, hyd_keys)
    return exp_hyd_keys


def atom_nuclear_charges(xgr):
    """ nuclear charges, by atom (connectivity-independent)
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_nuc_chg_dct = _transform_values(atm_sym_dct, func=_atom_nuclear_charge)
    return atm_nuc_chg_dct


def atom_total_valences(xgr):
    """ total valence electron counts, by atom (connectivity-independent)
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_tot_vlc_dct = _transform_values(atm_sym_dct, func=_atom_valence)
    return atm_tot_vlc_dct


def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    atm_ngb_keys_dct = {
        atm_key: _remove(atom_keys(atm_nbh), [atm_key])
        for atm_key, atm_nbh in atom_neighborhoods(xgr).items()}
    return atm_ngb_keys_dct


def atom_explicit_hydrogen_keys(xgr):
    """ explicit hydrogen valences, by atom
    """
    atm_exp_hyd_keys_dct = {
        atm_key: _remove(explicit_hydrogen_keys(atm_nbh), [atm_key])
        for atm_key, atm_nbh in atom_neighborhoods(xgr).items()}
    return atm_exp_hyd_keys_dct


def atom_bond_keys(xgr):
    """ bond keys, by atom
    """
    return _transform_values(atom_neighborhoods(xgr), bond_keys)


def atom_neighborhoods(xgr):
    """ bonded neighbor subgraphs, by atom
    """
    bnd_keys = bond_keys(xgr)

    def _neighborhood(atm_key):
        atm_bnd_keys = _filter(lambda x: atm_key in x, bnd_keys)
        return subgraph_by_bonds(xgr, atm_bnd_keys)

    atm_keys = atom_keys(xgr)
    atm_nbhs = tuple(map(_neighborhood, atm_keys))
    atm_nbh_dct = dict(zip(atm_keys, atm_nbhs))
    return atm_nbh_dct


# transformations
def implicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms implicit
    """
    bbn_keys = backbone_keys(xgr)
    atm_keys = bbn_keys if atm_keys is None else atm_keys
    assert set(atm_keys) <= set(bbn_keys)

    atm_imp_hyd_vlcs = _values_by_key(
        atom_implicit_hydrogen_valences(xgr), atm_keys)

    atm_exp_hyd_keys = _values_by_key(
        atom_explicit_hydrogen_keys(xgr), atm_keys)
    atm_exp_hyd_vlcs = tuple(map(len, atm_exp_hyd_keys))
    atm_tot_hyd_vlcs = numpy.add(atm_imp_hyd_vlcs, atm_exp_hyd_vlcs)

    exp_hyd_keys = tuple(_chain(*atm_exp_hyd_keys))

    xgr = set_atom_implicit_hydrogen_valences(
        xgr, dict(zip(atm_keys, atm_tot_hyd_vlcs)))
    xgr = delete_atoms(xgr, exp_hyd_keys)
    return xgr


def explicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms explicit
    """
    bbn_keys = backbone_keys(xgr)
    atm_keys = bbn_keys if atm_keys is None else atm_keys
    assert set(atm_keys) <= set(bbn_keys)

    atm_imp_hyd_vlcs = _values_by_key(
        atom_implicit_hydrogen_valences(xgr), atm_keys)

    xgr = set_atom_implicit_hydrogen_valences(
        xgr, _by_key({}, atm_keys, fill_val=0))
    xgr = add_explicit_hydrogens(
        xgr, dict(zip(atm_keys, atm_imp_hyd_vlcs)))
    return xgr


def explicit_stereo_sites(sgr):
    """ make the hydrogens at atom and bond stereo sites explicit
    """
    atm_ste_keys = atom_stereo_keys(sgr)
    bnd_ste_keys = bond_stereo_keys(sgr)
    bnd_ste_atm_keys = tuple(set(_chain(*bnd_ste_keys)))
    ste_atm_keys = atm_ste_keys + bnd_ste_atm_keys
    return explicit(sgr, atm_keys=ste_atm_keys)


def delete_atoms(xgr, atm_keys):
    """ delete atoms from the molecular graph
    """
    all_atm_keys = set(atom_keys(xgr))
    atm_keys = set(atm_keys)
    assert atm_keys <= all_atm_keys
    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(xgr, atm_keys_left)


def add_explicit_hydrogens(xgr, atm_exp_hyd_vlc_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_vlc_dct.keys()) <= set(atom_keys(xgr))
    for atm_key, atm_exp_hyd_vlc in atm_exp_hyd_vlc_dct.items():
        next_atm_key = max(atom_keys(xgr)) + 1
        atm_exp_hyd_keys = tuple(range(
            next_atm_key, next_atm_key + atm_exp_hyd_vlc))
        xgr = add_atoms(
            xgr, _by_key({}, atm_exp_hyd_keys, fill_val='H'))
        xgr = add_bonds(
            xgr, [frozenset({atm_key, atm_exp_hyd_key})
                  for atm_exp_hyd_key in atm_exp_hyd_keys])
    return xgr


def subgraph(xgr, atm_keys):
    """ the subgraph induced by a subset of the atoms
    """
    assert set(atm_keys) <= set(atom_keys(xgr))
    bnd_keys = [bnd_key for bnd_key in bond_keys(xgr)
                if set(bnd_key) <= set(atm_keys)]
    atm_dct = _by_key(atoms(xgr), atm_keys)
    bnd_dct = _by_key(bonds(xgr), bnd_keys)
    sub_xgr = (atm_dct, bnd_dct)
    return sub_xgr


def subgraph_by_bonds(xgr, bnd_keys):
    """ the subgraph induced by a subset of the bonds
    """
    atm_keys = set(_chain(*bnd_keys))
    atm_dct = _by_key(atoms(xgr), atm_keys)
    bnd_dct = _by_key(bonds(xgr), bnd_keys)
    sub_xgr = (atm_dct, bnd_dct)
    return sub_xgr


def relabel(xgr, atm_key_dct):
    """ relabel the graph with new atom keys
    """
    assert set(atm_key_dct.keys()) == set(atom_keys(xgr))

    _relabel_atom_key = atm_key_dct.__getitem__

    def _relabel_bond_key(bnd_key):
        return frozenset(map(_relabel_atom_key, bnd_key))

    atm_dct = _transform_keys(atoms(xgr), _relabel_atom_key)
    bnd_dct = _transform_keys(bonds(xgr), _relabel_bond_key)
    xgr = (atm_dct, bnd_dct)
    return xgr


def reflection(sgr):
    """ stereo graph reflection (inverts atom parities)
    """
    atm_pars = _filter_by_value(atom_stereo_parities(sgr),
                                lambda val: val is not None)
    refl_atm_pars = _transform_values(atm_pars, lambda val: not val)
    return set_atom_stereo_parities(sgr, refl_atm_pars)


# comparisons
def backbone_isomorphic(xgr1, xgr2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(xgr1, xgr2) is not None


def backbone_isomorphism(xgr1, xgr2):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `xgr1` to produce `xgr2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    xgr1 = implicit(xgr1)
    xgr2 = implicit(xgr2)
    nxg1 = _nxg_from_graph(xgr1)
    nxg2 = _nxg_from_graph(xgr2)
    iso_dct = _nxg_isomorphism(nxg1, nxg2)
    return iso_dct
