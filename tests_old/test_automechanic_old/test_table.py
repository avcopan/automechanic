""" test the automechanic.table module
"""
import pytest
import numpy
from pandas.testing import assert_frame_equal
from automechanic_old import table2 as table


def test__column_types():
    """ test table.column_types()
    """
    keys = ('a', 'b', 'c')
    typs = (table.INT_DTYPE, table.FLOAT_DTYPE, table.STR_DTYPE)
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    tbl = table.from_data(rows, keys, typs)
    assert table.column_types(tbl) == typs


def test__column_keys():
    """ test table.column_keys()
    """
    keys = ('a', 'b', 'c')
    typs = (table.INT_DTYPE, table.FLOAT_DTYPE, table.STR_DTYPE)
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    tbl = table.from_data(rows, keys, typs)
    assert table.column_keys(tbl) == keys


def test__row_indices():
    """ test table.column_keys()
    """
    keys = ('a', 'b', 'c')
    typs = (table.INT_DTYPE, table.FLOAT_DTYPE, table.STR_DTYPE)
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    tbl1 = table.from_data(rows, keys, typs)
    assert table.row_indices(tbl1) == (0, 1, 2)

    row_idxs = (5, 8, 17)
    tbl2 = table.from_data(rows, keys, typs, row_idxs)
    assert table.row_indices(tbl2) == row_idxs


def test__sql_where_in():
    """ test table.sql_where_in
    """
    keys = ('a', 'b', 'c')
    typs = (table.INT_DTYPE, table.FLOAT_DTYPE, table.STR_DTYPE)
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    row_idxs = (5, 6, 17)
    tbl = table.from_data(rows, keys, typs, row_idxs)
    print(tbl)
    print(table.sql_where_in(tbl, 'a', (1, 3)))


if __name__ == '__main__':
    test__column_types()
    test__column_keys()
    test__row_indices()
    test__sql_where_in()
