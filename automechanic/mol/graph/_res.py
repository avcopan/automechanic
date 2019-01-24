""" specific resonance graph functions
"""
from functools import partial as _partial
from itertools import product as _product
import numpy
from ._dict import values_by_key as _values_by_key
from ._base import atom_keys as _atom_keys
from ._base import bond_keys as _bond_keys
from ._base import bond_orders as _bond_orders
from ._base import set_bond_orders as _set_bond_orders
from ._base import atom_bond_keys as _atom_bond_keys
from ._base import atom_neighborhoods as _atom_neighborhoods
from ._base import atom_total_valences as _atom_total_valences
from ._base import (atom_implicit_hydrogen_valences as
                    _atom_implicit_hydrogen_valences)


def increment_bond_orders(rgr, bnd_ord_inc_dct):
    """ increment bond ordersby a given amount
    """
    bnd_keys = _bond_keys(rgr)
    bnd_ords = _values_by_key(_bond_orders(rgr), bnd_keys)
    bnd_ord_incs = _values_by_key(bnd_ord_inc_dct, bnd_keys, fill_val=0)
    new_bnd_ords = numpy.add(bnd_ords, bnd_ord_incs)

    bnd_ord_dct = dict(zip(bnd_keys, new_bnd_ords))
    rgr = _set_bond_orders(rgr, bnd_ord_dct)
    return rgr


def atom_bond_valences(rgr):
    """ bond valences, by atom
    """
    atm_keys = _atom_keys(rgr)
    atm_nbhs = _values_by_key(_atom_neighborhoods(rgr), atm_keys)
    atm_exp_bnd_vlcs = [sum(_bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_imp_hyd_vlcs = _values_by_key(
        _atom_implicit_hydrogen_valences(rgr), atm_keys)
    atm_bnd_vlcs = numpy.add(atm_exp_bnd_vlcs, atm_imp_hyd_vlcs)
    atm_bnd_vlc_dct = dict(zip(atm_keys, atm_bnd_vlcs))
    return atm_bnd_vlc_dct


def atom_radical_valences(rgr):
    """ radical valences, by atom
    """
    atm_keys = _atom_keys(rgr)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_tot_vlcs = _values_by_key(_atom_total_valences(rgr), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    return dict(zip(atm_keys, atm_rad_vlcs))


def maximum_spin_multiplicity(rgr):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), _atom_keys(rgr))
    return sum(atm_rad_vlcs) + 1


def possible_spin_multiplicities(rgr):
    """ possible spin multiplicities for this molecular graph
    """
    mult_max = maximum_spin_multiplicity(rgr)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


def subresonances(rgr):
    """ this graph and its lower-spin resonances
    """
    def _inc_range(max_bnd_ord_inc):
        return tuple(range(0, max_bnd_ord_inc+1))

    _incremented_graph = _partial(increment_bond_orders, rgr)

    atm_keys = _atom_keys(rgr)
    bnd_keys = _bond_keys(rgr)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), atm_keys)
    atm_bnd_keys_lst = _values_by_key(_atom_bond_keys(rgr), atm_keys)
    max_bnd_ord_incs = _values_by_key(_maximum_bond_increments(rgr), bnd_keys)

    def _is_valid(bnd_ord_inc_dct):
        # check if radical decrements are less than radical valences
        def __tally(atm_bnd_keys):
            return sum(_values_by_key(bnd_ord_inc_dct, atm_bnd_keys))
        atm_rad_vlc_decs = tuple(map(__tally, atm_bnd_keys_lst))
        ret = numpy.all(numpy.less_equal(atm_rad_vlc_decs, atm_rad_vlcs))
        return ret

    def _bond_value_dictionary(bnd_vals):
        return dict(zip(bnd_keys, bnd_vals))

    bnd_ord_incs_itr = _product(*map(_inc_range, max_bnd_ord_incs))
    bnd_ord_inc_dct_itr = map(_bond_value_dictionary, bnd_ord_incs_itr)
    bnd_ord_inc_dct_itr = filter(_is_valid, bnd_ord_inc_dct_itr)
    rgrs = tuple(map(_incremented_graph, bnd_ord_inc_dct_itr))
    return rgrs


def lowspin_resonance(rgr):
    """ get the resonance graph with the lowest maximum spin
    """
    return min(subresonances(rgr), key=maximum_spin_multiplicity)


def _maximum_bond_increments(rgr):
    """ maximum possible order increments for each bond in the graph
    """
    atm_rad_vlc_dct = atom_radical_valences(rgr)

    def _max_increment(bnd_key):
        return min(map(atm_rad_vlc_dct.__getitem__, bnd_key))

    bnd_keys = _bond_keys(rgr)
    max_bnd_ord_incs = tuple(map(_max_increment, bnd_keys))
    max_bnd_ord_inc_dct = dict(zip(bnd_keys, max_bnd_ord_incs))
    return max_bnd_ord_inc_dct


# def _filter_repeats(rgrs):
#     return tuple(map(_unfreeze, sorted(set(map(_freeze, rgrs)))))
#
#
# def _freeze(rgr):
#     atm_keys = _atom_keys(rgr)
#     bnd_keys = _bond_keys(rgr)
#     atm_vals = _values_by_key(_atoms(rgr), atm_keys)
#     bnd_vals = _values_by_key(_bonds(rgr), bnd_keys)
#     atm_itms = tuple(zip(atm_keys, atm_vals))
#     bnd_itms = tuple(zip(bnd_keys, bnd_vals))
#     frz_rgr = (atm_itms, bnd_itms)
#     return frz_rgr
#
#
# def _unfreeze(frz_rgr):
#     atm_itms, bnd_itms = frz_rgr
#     atm_dct = dict(atm_itms)
#     bnd_dct = dict(bnd_itms)
#     rgr = (atm_dct, bnd_dct)
#     return rgr
