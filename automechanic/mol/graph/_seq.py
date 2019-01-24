""" sequence helpers
"""
from itertools import filterfalse as _filterfalse


def filter_(func, seq):
    """ imitates built-in `filter`, returning tuples
    """
    return tuple(filter(func, seq))


def filterfalse(func, seq):
    """ imitates `itertools.filterfalse`, returning tuples
    """
    return tuple(_filterfalse(func, seq))


def remove(seq1, seq2):
    """ remove the elements in `seq2` from `seq1`
    """
    return tuple(elem for elem in seq1 if elem not in seq2)
