""" .xyz-based functions
"""
import warnings
import re
from .parse import WHITESPACES
from .parse import INTEGER
from .parse import STRING_START
from .parse import LINE_START
from .parse import LINE_END
from .parse import LETTER
from .parse import FLOAT
from .parse import maybe
from .parse import named_capture


def number_of_atoms(dxyz):
    """ number of atoms from a .xyz string
    """
    natms_pattern = maybe(WHITESPACES).join(
        [STRING_START, named_capture(INTEGER, 'natms'), LINE_END])
    match = re.search(natms_pattern, dxyz, re.MULTILINE)
    assert match
    gdct = match.groupdict()
    natms = int(gdct['natms'])
    return natms


def geometry(dxyz):
    """ geometry from a .xyz string
    """
    natms = number_of_atoms(dxyz)
    atomic_symbol = LETTER + maybe(LETTER)
    atom_pattern = WHITESPACES.join(
        [named_capture(atomic_symbol, 'asymb'), named_capture(FLOAT, 'x'),
         named_capture(FLOAT, 'y'), named_capture(FLOAT, 'z')])
    line_pattern = LINE_START + atom_pattern + LINE_END

    mgeo = []
    for match in re.finditer(line_pattern, dxyz, re.MULTILINE):
        gdct = match.groupdict()
        asymb = gdct['asymb']
        xyz = tuple(map(float, [gdct['x'], gdct['y'], gdct['z']]))
        mgeo.append((asymb, xyz))

    if len(mgeo) != natms:
        warnings.warn("\nThis .xyz file is inconsistent: {:s}".format(dxyz))

    return tuple(mgeo)


if __name__ == '__main__':
    DXYZ = """4

C          1.19654        0.06238        0.04613
C          0.71372        0.14642       -1.19256
C          2.10909        0.13478       -1.02095
H          1.86300        0.00588        0.87893
H          0.64110        0.21855       -2.25578
    """
    print geometry(DXYZ)
    print repr(DXYZ)