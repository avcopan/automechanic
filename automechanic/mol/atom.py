""" functions operating on atomic symbols
"""

SYMBOLS = ['H', 'HE',
           'C',
           'N',
           'O', 'S',
           'F', 'CL',
           'NE', 'AR']


def valence(sym):
    """ bonding valence
    """
    return {'H': 1, 'HE': 0,
            'C': 4,
            'N': 3,
            'O': 2, 'S': 2,
            'F': 1, 'CL': 1,
            'NE': 0, 'AR': 0}[sym.upper()]


def nuclear_charge(sym):
    """ nuclear charge
    """
    return {'H': 1, 'HE': 2,
            'C': 6,
            'N': 7,
            'O': 8, 'S': 16,
            'F': 9, 'CL': 17,
            'NE': 10, 'AR': 18}[sym.upper()]


def lone_pair_count(sym):
    """ lone pair count
    """
    return {'H': 0, 'HE': 1,
            'C': 0,
            'N': 1,
            'O': 2, 'S': 2,
            'F': 3, 'CL': 3,
            'NE': 4, 'AR': 4}[sym.upper()]
