""" interface to pandas for tabular data
"""
from collections import OrderedDict as _OrderedDict
from collections import Sequence as _Sequence
from itertools import starmap as _starmap
from itertools import chain as _chain
import numpy
from pandas import Series as _Series
from pandas import merge as _merge
from pandas import read_csv as _read_csv
from pandas import DataFrame as _DataFrame
from pandas import Int64Index as _Int64Index
from pandas import RangeIndex as _RangeIndex
from autoparse.find import first_capture as _first_capture
from autoparse.pattern import escape as _escape
from autoparse.pattern import capturing as _capturing
from autoparse.pattern import STRING_START as _STRING_START
from autoparse.pattern import STRING_END as _STRING_END
from autoparse.pattern import UNSIGNED_INTEGER as _UNSIGNED_INTEGER

IDX_KEY = 'i_'
IDX_TYP = numpy.dtype('int64')

IDX_SAVE_KEY_FORMAT = 'i{:d}_'
NAN = numpy.nan


# table creation
def from_records(vals, keys, typs=None, idxs=None):
    """ create a table from a series of records
    """
    if any(map(_is_fancy, keys)):
        ret = _from_fancy_records(vals, keys, typs, idxs)
    else:
        ret = _from_records(vals, keys, typs, idxs)
    return ret


def from_starmap(tbl, func, arg_keys, keys, typs=None):
    """ starmap columns and create a table from the resulting values
    """
    idxs_itr, fancy_vals_itr = zip(*enum_(tbl, arg_keys))
    idxs = list(idxs_itr)
    vals = list(_starmap(func, fancy_vals_itr))
    return from_records(vals, keys, typs, idxs)


def enforce_schema(tbl, keys, typs=None):
    """ make sure the table adheres to a specific schema
    """
    if keys is not None:
        assert has_keys(tbl, keys=keys)
    if typs is not None:
        assert keys is not None
        tbl = set_typs(tbl, keys=keys, typs=typs)
    return tbl


# table I/O
def read_csv(file_pth):
    """ read table from a CSV file
    """
    tbl = _read_csv(file_pth)
    if IDX_KEY in keys_(tbl):
        tbl.set_index(IDX_KEY, inplace=True)
    else:
        tbl.index.rename(IDX_KEY, inplace=True)
    return tbl


def write_csv(file_pth, tbl, float_format=None):
    """ write table to a CSV file
    """
    tbl.to_csv(file_pth, float_format=float_format)


# table properties
def idxs_(tbl):
    """ table indices
    """
    return tuple(tbl.index)


def keys_(tbl):
    """ table keys
    """
    return tuple(tbl.columns)


def vals_(tbl):
    """ table values
    """
    return numpy.array(tbl.values)


def typs_(tbl):
    """ table column types
    """
    return tuple(tbl.dtypes)


# data types
def dt_(typ):
    """ return data type conversion
    """
    np_typ = numpy.dtype(typ)
    if numpy.issubdtype(np_typ, numpy.bool_):
        typ = numpy.dtype('bool')
    elif numpy.issubdtype(np_typ, numpy.integer):
        typ = numpy.dtype('int64')
    elif numpy.issubdtype(np_typ, numpy.number):
        typ = numpy.dtype('float64')
    else:
        typ = numpy.dtype('O')
    return typ


def set_typs(tbl, keys, typs):
    """ set one or more column types
    """
    assert has_keys(tbl, keys)
    ukeys = _unfancy_keys(keys)
    utyps = _unfancy_group_types(typs, keys)
    tbl = tbl.astype(dict(zip(ukeys, utyps)))
    return tbl


# iterators
def iter_(tbl, keys):
    """ yield records from one or more columns (fancy)
    """
    itrs = map(_iter_vals, _fancy_select(tbl, keys))
    for fancy_vals in zip(*itrs):
        yield fancy_vals


def enum_(tbl, keys):
    """ yield records from one or more columns (fancy), with index
    """
    itrs = map(_iter_vals, _fancy_select(tbl, keys))
    for idx, fancy_vals in zip(_iter_idxs(tbl), zip(*itrs)):
        yield idx, fancy_vals


# misc
def equal(tbl1, tbl2):
    """ returns two if tables are the same, including datatypes
    """
    return (numpy.array_equal(keys_(tbl1), keys_(tbl2)) and
            numpy.array_equal(typs_(tbl1), typs_(tbl2)) and
            numpy.array_equal(idxs_(tbl1), idxs_(tbl2)) and
            numpy.array_equal(vals_(tbl1), vals_(tbl2)))


def has_keys(tbl, keys):
    """ does the table have these keys?
    """
    ukeys = _unfancy_keys(keys)
    return all(ukey in keys_(tbl) for ukey in ukeys)


def update(tbl1, tbl2):
    """ write the values in table 2 over those in table 1

    indices of table 2 must be a subset of those in table 1, but it may
    contain new columns to be added to table 1
    """
    assert set(idxs_(tbl1)) >= set(idxs_(tbl2))
    tbl = tbl1.copy(deep=True)
    for key in tbl2:
        tbl.loc[list(idxs_(tbl2)), key] = tbl2[key]
    return tbl


def left_join(tbl1, tbl2, key=None):
    """ SQL LEFT OUTER join of two tables

    Keeps keys from the left table and preserves their order
    """
    return _join(tbl1, tbl2, key=key, right=False)


def right_join(tbl1, tbl2, key=None):
    """ SQL RIGHT OUTER join of two tables

    Keeps keys from the right table and preserves their order
    """
    return _join(tbl1, tbl2, key=key, right=True)


def _join(tbl1, tbl2, key=None, right=False):
    key = IDX_KEY if key is None else key
    how = 'right' if right else 'left'
    tbl = _merge(tbl1, tbl2, how=how, left_on=key, right_on=key)
    tbl.index.rename(IDX_KEY, inplace=True)
    return tbl


def save_index(tbl):
    """ save the index as a regular column
    """
    tbl = _DataFrame.copy(tbl)
    idx_save_key = next_index_save_key(tbl)
    tbl.insert(loc=0, column=idx_save_key, value=idxs_(tbl))
    return tbl


def next_index_save_key(tbl):
    """ get the value of the next saved index key
    """
    sids = _index_save_ids(tbl)
    sid = max(sids, default=-1) + 1
    idx_save_key = IDX_SAVE_KEY_FORMAT.format(sid)
    return idx_save_key


def _index_save_ids(tbl):
    """ get the values of previously saved index keys
    """
    _pattern = (_STRING_START + 'i' + _capturing(_UNSIGNED_INTEGER) +
                _escape('_') + _STRING_END)
    keys = keys_(tbl)
    caps = tuple(_first_capture(_pattern, key) for key in keys)
    sids = tuple(int(cap) for cap in caps if cap is not None)
    return sids


def group_dictionary(tbl, key1, key2):
    """ get a dictionary from some column(s) onto others
    """
    grps = tbl.groupby(_fancy_selector(key1))
    dct = {grp_key: tuple(iter_(grp_tbl, key2)) for grp_key, grp_tbl in grps}
    return dct


# helpers
def _iter_vals(tbl):
    """ iterate over values in a table or series
    """
    if hasattr(tbl, 'itertuples'):
        for vals in tbl.itertuples(index=False, name=None):
            yield vals
    elif hasattr(tbl, 'iteritems'):
        for _, val in tbl.iteritems():
            yield val


def _iter_idxs(tbl):
    """ iterate over indices in a table or series
    """
    if hasattr(tbl, 'itertuples'):
        for vals in tbl.itertuples(index=True, name=None):
            idx = vals[0]
            yield idx
    elif hasattr(tbl, 'iteritems'):
        for idx, _ in tbl.iteritems():
            yield idx


def _from_records(vals, keys, typs=None, idxs=None):
    """ create a table from a series of records
    """
    if numpy.size(vals):
        assert numpy.ndim(vals) in (1, 2)
        vals = (numpy.array(vals) if numpy.ndim(vals) == 2 else
                numpy.reshape(vals, (-1, 1)))
    else:
        vals = numpy.reshape([], (0, len(keys)))
    nrows, ncols = numpy.shape(vals)

    idxs = (_RangeIndex(stop=nrows, name=IDX_KEY) if idxs is None
            else _Int64Index(idxs, name=IDX_KEY))
    typs = (object,) * ncols if typs is None else typs
    typs = tuple(map(dt_, typs))
    assert len(keys) == len(typs) == ncols and len(idxs) == nrows
    assert IDX_KEY not in keys

    cols = numpy.transpose(vals)
    data = _OrderedDict((key, _Series(data=col, dtype=typ, index=idxs))
                        for col, key, typ in zip(cols, keys, typs))

    return _DataFrame(data=data, index=idxs)


# fancy helpers
def _fancy_select(tbl, fancy_keys):
    """ select individual columns and groups of columns
    """
    return tuple(tbl[_fancy_selector(fancy_key)] for fancy_key in fancy_keys)


def _fancy_selector(fancy_key):
    """ if this is a sequency of keys, return a list
    """
    return list(fancy_key) if _is_fancy(fancy_key) else fancy_key


def _is_fancy(key):
    assert isinstance(key, _Sequence)
    return (isinstance(key, _Sequence)
            and not isinstance(key, (str, bytes, bytearray)))


def _from_fancy_records(fancy_vals, fancy_keys, group_typs=None, idxs=None):
    """ flatten fancily-shaped values and keys to create a table from records
    """
    ngroups = len(fancy_keys)

    fancy_cols = tuple(zip(*fancy_vals))
    group_typs = (object,) * ngroups if group_typs is None else group_typs
    keys = _unfancy_keys(fancy_keys)
    typs = _unfancy_group_types(group_typs, fancy_keys)
    cols = _unfancy_columns(fancy_cols, fancy_keys) if fancy_cols else ()
    vals = numpy.transpose(cols)
    return _from_records(vals, keys, typs, idxs)


def _unfancy_keys(fancy_keys):
    keys = tuple(_chain(*(fkey if _is_fancy(fkey) else [fkey]
                          for fkey in fancy_keys)))
    return keys


def _unfancy_group_types(group_typs, fancy_keys):
    assert len(group_typs) == len(fancy_keys)
    typs = tuple(_chain(*([typ] * len(fkey) if _is_fancy(fkey) else [typ]
                          for typ, fkey in zip(group_typs, fancy_keys))))
    return typs


def _unfancy_columns(fancy_cols, fancy_keys):
    assert len(fancy_cols) == len(fancy_keys)
    cols = tuple(_chain(*(zip(*fcol) if _is_fancy(fkey) else [fcol]
                          for fcol, fkey in zip(fancy_cols, fancy_keys))))
    return cols
