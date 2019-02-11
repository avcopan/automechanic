""" test the automechanic.dxyz module
"""
import pytest
from automechanic_old import dotxyz


def test__number_of_atoms():
    """ test dotxyz.number_of_atoms()
    """
    dxyz = ('4\n'
            '\n'
            'C          1.19654        0.06238        0.04613\n'
            'C          0.71372        0.14642       -1.19256\n'
            'C          2.10909        0.13478       -1.02095\n'
            'H          1.86300        0.00588        0.87893\n'
            'H          0.64110        0.21855       -2.25578\n')
    assert dotxyz.number_of_atoms(dxyz) == 4


def test__geometry():
    """ test dotxyz.geometry()
    """
    dxyz = ('5\n'
            '\n'
            'C          1.19654        0.06238        0.04613\n'
            'C          0.71372        0.14642       -1.19256\n'
            'C          2.10909        0.13478       -1.02095\n'
            'H          1.86300        0.00588        0.87893\n'
            'H          0.64110        0.21855       -2.25578\n')

    assert (dotxyz.geometry(dxyz)
            == (('C', (1.19654, 0.06238, 0.04613)),
                ('C', (0.71372, 0.14642, -1.19256)),
                ('C', (2.10909, 0.13478, -1.02095)),
                ('H', (1.863, 0.00588, 0.87893)),
                ('H', (0.6411, 0.21855, -2.25578))))

    bad_dxyz = ('4\n'
                '\n'
                'C          1.19654        0.06238        0.04613\n'
                'C          0.71372        0.14642       -1.19256\n'
                'C          2.10909        0.13478       -1.02095\n'
                'H          1.86300        0.00588        0.87893\n'
                'H          0.64110        0.21855       -2.25578\n')

    with pytest.raises(ValueError):
        dotxyz.geometry(bad_dxyz)


if __name__ == '__main__':
    test__geometry()
