""" tuple dictionary helpers

dictionary values must all be iterables of the same length
"""
import numpy
from ._dict import transform_values as _transform_values


def position_count(dct):
    """ the number of positions in this tuple dictionary
    """
    pos_cnt = numpy.shape(list(dct.values()))[1] if dct else None
    return pos_cnt


def by_key_by_position(dct, keys, pos):
    """ position values, as a dictionary with these keys
    """
    pos_dct = {}
    if keys:
        assert set(keys) <= set(dct.keys())
        assert pos <= position_count(dct)
        vals = tuple(map(dct.__getitem__, keys))
        pos_vals = list(zip(*vals))[pos]
        pos_dct = dict(zip(keys, pos_vals))
    return pos_dct


def set_by_key_by_position(dct, pos_dct, pos):
    """ set values by position and key
    """
    if pos_dct:
        assert set(pos_dct.keys()) <= set(dct.keys())
        assert pos <= position_count(dct)
        dct = _transform_values(dct, func=list)
        for key, pos_val in pos_dct.items():
            dct[key][pos] = pos_val
        dct = _transform_values(dct, func=tuple)
    return dct
