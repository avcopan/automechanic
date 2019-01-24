""" dictionary helpers
"""


def by_key(dct, keys, fill_val=None):
    """ dictionary on a set of keys, filling missing entries
    """
    return dict(zip(keys, values_by_key(dct, keys, fill_val=fill_val)))


def values_by_key(dct, keys, fill_val=None):
    """ return dictionary values for specific keys, filling missing entries
    """
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def keys_by_value(dct, vals):
    """ return dictionary keys for specific values
    """
    return tuple(key for key, val in dct.items() if val in vals)


def transform_keys(dct, func):
    """ apply a function to each key
    """
    return dict(zip(map(func, dct.keys()), dct.values()))


def transform_values(dct, func):
    """ apply a function to each value
    """
    return dict(zip(dct.keys(), map(func, dct.values())))


def keys_sorted_by_value(dct):
    """ dictionary keys sorted by their associated values
    """
    return tuple(key for key, _ in sorted(dct.items(), key=lambda x: x[1]))


# def filter_by_key(dct, func):
#     """ filter dictionary entries by their keys
#     """
#     return {key: val for key, val in dct.items() if func(key)}


def filter_by_value(dct, func=lambda val: val):
    """ filter dictionary entries by their values
    """
    return {key: val for key, val in dct.items() if func(val)}
