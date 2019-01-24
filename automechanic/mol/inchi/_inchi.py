""" functions operating on InChI strings
"""
import itertools
from string import ascii_lowercase as _ascii_lowercase
from ._rdkit import from_inchi as _rdm_from_inchi
from ._rdkit import to_inchi as _rdm_to_inchi
from ._rdkit import to_smiles as _rdm_to_smiles
from ._rdkit import inchi_to_inchi_key as _inchi_to_inchi_key
from ._rdkit import geometry as _rdm_to_geometry
from ._rdkit import connectivity_graph as _rdm_to_connectivity_graph
from ._pybel import from_inchi as _pbm_from_inchi
from ._pybel import geometry as _pbm_to_geometry
from ..geom import inchi as _inchi_from_geometry
from ..graph import inchi as _inchi_from_graph
from ..graph import stereo_inchi as _inchi_from_stereo_graph
from ...rere.pattern import escape as _escape
from ...rere.pattern import named_capturing as _named_capturing
from ...rere.pattern import one_or_more as _one_or_more
from ...rere.pattern import one_of_these as _one_of_these
from ...rere.pattern import not_followed_by as _not_followed_by
from ...rere.pattern_lib import LOWERCASE_LETTER as _LOWERCASE_LETTER
from ...rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ...rere.pattern_lib import NONWHITESPACE as _NONWHITESPACE
from ...rere.pattern_lib import STRING_START as _STRING_START
from ...rere.pattern_lib import STRING_END as _STRING_END
from ...rere.find import first_named_capture as _first_named_capture
from ...rere.find import all_captures as _all_captures

_NONWHITESPACES_NONGREEDY = _one_or_more(_NONWHITESPACE, greedy=False)
_INCHI_SUBLAYER_END = _one_of_these([_escape('/'), _STRING_END])
_STEREO_UNKNOWN_VAL = 'u'
_STEREO_UNDEFINED_VAL = '?'
_STEREO_MINUS_VAL = '-'
_STEREO_PLUS_VAL = '+'


def _key_layer(key):
    assert key in _ascii_lowercase

    class _KEYxLAYER():
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = key + _named_capturing(_NONWHITESPACES_NONGREEDY,
                                        name=CONTENT_KEY)
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    _KEYxLAYER.__name__ = 'KEYxLAYER(\'{:s}\')'.format(key)
    return _KEYxLAYER


class PARSE():
    """ InChI format specifications """

    class PREFIX():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _STRING_START
        _LAYER = (_escape('InChI=') +
                  _named_capturing(_NONWHITESPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    class FORMULA():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = (_not_followed_by(_LOWERCASE_LETTER) +
                  _named_capturing(_NONWHITESPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    KEY_LAYER = _key_layer

    class ATOMxSTEREO(_key_layer('t')):
        """ _ """

        class TERM():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = _UNSIGNED_INTEGER
            _VAL = _one_of_these(list(map(_escape, VALS)))
            PATTERN = (_named_capturing(_KEY, name=KEY_KEY) +
                       _named_capturing(_VAL, name=VAL_KEY))

    class BONDxSTEREO(_key_layer('b')):
        """ _ """

        class TERM():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = _UNSIGNED_INTEGER + _escape('-') + _UNSIGNED_INTEGER
            _VAL = _one_of_these(list(map(_escape, VALS)))
            PATTERN = (_named_capturing(_KEY, name=KEY_KEY) +
                       _named_capturing(_VAL, name=VAL_KEY))


def smiles(ich):
    """ SMILES string from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    smi = _rdm_to_smiles(rdm)
    return smi


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdm_from_inchi(ich)
    ich = _rdm_to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def is_closed(ich):
    """ regenerating the InChI string yields the same thing
    """
    return recalculate(ich) == ich


def prefix(ich):
    """ InChI prefix
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    pfx = cap_dct[PARSE.PREFIX.LAYER_KEY]
    return pfx


def version(ich):
    """ InChI version
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    ver = cap_dct[PARSE.PREFIX.CONTENT_KEY]
    return ver


def formula_layer(ich):
    """ InChI formula
    """
    cap_dct = _first_named_capture(PARSE.FORMULA.PATTERN, ich)
    assert cap_dct
    fml = cap_dct[PARSE.FORMULA.LAYER_KEY]
    return fml


def key_layer(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = PARSE.KEY_LAYER(key)
    cap_dct = _first_named_capture(key_layer_parser.PATTERN, ich)
    return cap_dct[key_layer_parser.LAYER_KEY] if cap_dct else None


def key_layer_content(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = PARSE.KEY_LAYER(key)
    cap_dct = _first_named_capture(key_layer_parser.PATTERN, ich)
    return cap_dct[key_layer_parser.CONTENT_KEY] if cap_dct else None


def core_parent(ich):
    """ get the InChI string of the core parent structure
    """
    lyrs = [prefix(ich), formula_layer(ich)]
    for key in ('c', 'h'):
        lyr = key_layer(ich, key=key)
        if lyr is not None:
            lyrs.append(lyr)
    return '/'.join(lyrs)


def atom_stereo_elements(ich):
    """ atom stereo keys and values
    """
    cap_dct = _first_named_capture(PARSE.ATOMxSTEREO.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[PARSE.ATOMxSTEREO.LAYER_KEY]
        ret = _all_captures(PARSE.ATOMxSTEREO.TERM.PATTERN, lyr)
    return ret


def bond_stereo_elements(ich):
    """ bond stereo keys and values
    """
    cap_dct = _first_named_capture(PARSE.BONDxSTEREO.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[PARSE.BONDxSTEREO.LAYER_KEY]
        ret = _all_captures(PARSE.BONDxSTEREO.TERM.PATTERN, lyr)
    return ret


def has_unknown_stereo_elements(ich):
    """ does this InChI string have unknown stereo elements?
    """
    ich_ste = recalculate(ich, force_stereo=True)
    return (atom_stereo_elements(ich_ste) !=
            _known_atom_stereo_elements(ich_ste) or
            bond_stereo_elements(ich_ste) !=
            _known_bond_stereo_elements(ich_ste))


def compatible_stereoisomers(ich):
    """ expand InChI string to its compatible stereoisomers
    """
    atm_terms_lst = []
    for num, (key, val) in enumerate(
            atom_stereo_elements(recalculate(ich, force_stereo=True))):
        terms = ([key + val] if val not in PARSE.ATOMxSTEREO.TERM.UNKNOWN_VALS
                 else [key + PARSE.ATOMxSTEREO.TERM.MINUS_VAL] if num == 0
                 else [key + PARSE.ATOMxSTEREO.TERM.MINUS_VAL,
                       key + PARSE.ATOMxSTEREO.TERM.PLUS_VAL])
        atm_terms_lst.append(terms)

    bnd_terms_lst = []
    for key, val in bond_stereo_elements(recalculate(ich, force_stereo=True)):
        terms = ([key + val] if val not in PARSE.BONDxSTEREO.TERM.UNKNOWN_VALS
                 else [key + PARSE.BONDxSTEREO.TERM.MINUS_VAL,
                       key + PARSE.BONDxSTEREO.TERM.PLUS_VAL])
        bnd_terms_lst.append(terms)

    ich_cp = core_parent(ich)
    ich_lst = [ich_cp]
    if bnd_terms_lst:
        lyr_lst = ['b' + ','.join(terms)
                   for terms in itertools.product(*bnd_terms_lst)]
        ich_lst = [ich_start + '/' + lyr
                   for ich_start, lyr in itertools.product(ich_lst, lyr_lst)]

    if atm_terms_lst:
        lyr_lst = ['t' + ','.join(terms)
                   for terms in itertools.product(*atm_terms_lst)]
        ich_lst = [ich_start + '/' + lyr
                   for ich_start, lyr in itertools.product(ich_lst, lyr_lst)]

    ich_lst = tuple(map(recalculate, ich_lst))
    return ich_lst


def inchi_key(ich):
    """ computes InChIKey from an InChI string
    """
    return _inchi_to_inchi_key(ich)


def connectivity_graph(ich):
    """ connectivity graph from an InChI string
    """
    rdm = _rdm_from_inchi(ich)

    # make sure the InChI string was valid
    ich_, _ = _rdm_to_inchi(rdm, with_aux_info=True)
    assert core_parent(ich) == core_parent(ich_)

    cgr = _rdm_to_connectivity_graph(rdm)
    cgr_ich = _inchi_from_graph(cgr)
    assert _has_same_connectivity(ich, cgr_ich)
    return cgr


def stereo_graph(ich):
    """ stereo graph from an InChI string
    """
    def _int_minus_one(int_str):
        return int(int_str) - 1

    def _atom_key(ich_atm_key):
        return _int_minus_one(ich_atm_key)

    def _bond_key(ich_bnd_key):
        return frozenset(map(_int_minus_one, str.split(ich_bnd_key, '-')))

    def _value(ich_ste_val):
        assert ich_ste_val in ('-', '+')
        return ich_ste_val == '+'

    atms, cnns = connectivity_graph(ich)
    assert not has_unknown_stereo_elements(ich)
    atm_ste_dct = {_atom_key(key): _value(val)
                   for key, val in atom_stereo_elements(ich)}
    bnd_ste_dct = {_bond_key(key): (1, _value(val))
                   for key, val in bond_stereo_elements(ich)}
    assert set(atm_ste_dct.keys()) <= set(range(len(atms)))
    assert set(bnd_ste_dct.keys()) <= set(cnns.keys())
    atms = {atm_key: ((sym, hcnt, atm_ste_dct[atm_key])
                      if atm_key in atm_ste_dct else (sym, hcnt, None))
            for atm_key, (sym, hcnt, _) in atms.items()}
    bnds = cnns.copy()
    bnds.update(bnd_ste_dct)
    sgr = (atms, bnds)
    sgr_ich = _inchi_from_stereo_graph(sgr)
    assert _has_same_connectivity(ich, sgr_ich)
    assert _has_compatible_stereo(ich, sgr_ich)
    return sgr


def geometry(ich):
    """ cartesian geometry from an InChI string
    """
    try:
        rdm = _rdm_from_inchi(ich)
        geo = _rdm_to_geometry(rdm)
        geo_ich = _inchi_from_geometry(geo)
        assert _has_same_connectivity(ich, geo_ich)
        assert _has_compatible_stereo(ich, geo_ich)
    except (AssertionError, RuntimeError):
        pbm = _pbm_from_inchi(ich)
        geo = _pbm_to_geometry(pbm)
        geo_ich = _inchi_from_geometry(geo)
        assert _has_same_connectivity(ich, geo_ich)
        assert _has_compatible_stereo(ich, geo_ich)
    return geo


def _has_same_connectivity(ich, other_ich):
    """ do these InChI strings have the same connectivity?
    """
    return (key_layer(ich, 'c') == key_layer(other_ich, 'c') and
            key_layer(ich, 'h') == key_layer(other_ich, 'h'))


def _has_compatible_stereo(ich, other_ich):
    """ is `other_ich` compatible with `ich`?
    """
    return (set(_known_atom_stereo_elements(ich)) <=
            set(_known_atom_stereo_elements(other_ich)) and
            set(_known_bond_stereo_elements(ich)) <=
            set(_known_bond_stereo_elements(other_ich)))


def _known_atom_stereo_elements(ich):
    return tuple((key, val) for key, val in atom_stereo_elements(ich)
                 if val not in PARSE.ATOMxSTEREO.TERM.UNKNOWN_VALS)


def _known_bond_stereo_elements(ich):
    return tuple((key, val) for key, val in bond_stereo_elements(ich)
                 if val not in PARSE.BONDxSTEREO.TERM.UNKNOWN_VALS)
