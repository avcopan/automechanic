""" CHEMKIN parsing
"""
from itertools import chain
from more_itertools import lstrip
from more_itertools import split_before
from numpy import multiply as _scale
from autoparse.pattern import maybe
from autoparse.pattern import escape
from autoparse.pattern import capturing
from autoparse.pattern import one_or_more
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_of_these
from autoparse.pattern import not_followed_by
from autoparse.pattern import STRING_START
from autoparse.pattern import LINE_START
from autoparse.pattern import LINE_END
from autoparse.pattern import WILDCARD
from autoparse.pattern import NONSPACE
from autoparse.pattern import LINESPACES
from autoparse.pattern import NEWLINE
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import DIGIT
from autoparse.pattern import UNSIGNED_FLOAT
from autoparse.pattern import PLUS
from autoparse.pattern import INTEGER
from autoparse.pattern import FLOAT
from autoparse.pattern import EXPONENTIAL_INTEGER
from autoparse.pattern import EXPONENTIAL_FLOAT
from autoparse.find import has_match as find_if_has_match
from autoparse.find import split as find_split
from autoparse.find import remove as find_remove
from autoparse.find import matcher as find_matcher
from autoparse.find import all_captures as find_captures
from autoparse.find import split_words as find_split_words
from autoparse.find import split_lines as find_split_lines
from autoparse.find import strip_spaces as find_strip_spaces
from autoparse.find import first_capture as find_first_capture

CHEMKIN_ARROW = maybe(escape('<')) + escape('=') + maybe(escape('>'))
CHEMKIN_PLUS_EM = PLUS + 'M'
CHEMKIN_PAREN_PLUS_EM = escape('(') + PLUS + 'M' + escape(')')
A_UNITS = (
    ('MOLECULES', 6.02214076e23),
    ('MOLES', 1.)
)
E_UNITS = (
    ('KCAL/MOLE', 1.e3),
    ('CAL/MOLE', 1.),
    ('KJOULES/MOLE', 0.239006 * 1.e3),
    ('JOULES/MOLE', 0.239006),
    ('KELVINS', 0.001987191686485529 * 1.e3)
)


def species_names(mech_str):
    """ find all species
    """
    block_str = species_block(mech_str)
    spcs = find_split_words(block_str)
    return spcs


def reaction_data(mech_str):
    """ find all reaction data
    """
    a_key, e_key = reaction_unit_names(mech_str)
    rxn_dstr_lst = reaction_data_strings(mech_str)

    rxn_lst = list(map(reaction_data_reaction_name, rxn_dstr_lst))
    arrh_lst = list(map(reaction_data_high_p_coeffs, rxn_dstr_lst))
    arrh_lst = _convert_units(arrh_lst, a_key=a_key, e_key=e_key)
    rxn_dat_lst = tuple(zip(rxn_lst, arrh_lst))
    return rxn_dat_lst


def _convert_units(arrh_lst, a_key, e_key):
    a_unit_dct = dict(A_UNITS)
    e_unit_dct = dict(E_UNITS)
    a_lst, b_lst, e_lst = zip(*arrh_lst)
    if a_key is not None:
        assert a_key in a_unit_dct
        a_lst = _scale(a_lst, a_unit_dct[a_key])
    if e_key is not None:
        assert e_key in e_unit_dct
        e_lst = _scale(e_lst, e_unit_dct[e_key])
    return list(zip(a_lst, b_lst, e_lst))


def reaction_data_strings(mech_str):
    """ find all reaction data strings
    """
    block_str = remove_blanks(reactions_block(mech_str))
    headline_pattern = CHEMKIN_ARROW
    rxn_dat_lst = _headlined_sections(headline_pattern, block_str)
    return rxn_dat_lst


def reaction_data_reaction_name(rxn_dstr):
    """ get the reaction name from a reaction data string
    """
    headline = find_split_lines(rxn_dstr)[0]
    number = one_of_these(
        [EXPONENTIAL_FLOAT, EXPONENTIAL_INTEGER, FLOAT, INTEGER])
    pattern = (LINESPACES + number + LINESPACES + number + LINESPACES + number)
    rxn = find_remove(pattern, headline)
    return rxn


def split_reaction_name(rxn):
    """ split a CHEMKIN reaction name into reactants and products
    """
    rct_str, prd_str = find_split(CHEMKIN_ARROW, rxn)
    rcts = _split_reagent_string(rct_str)
    prds = _split_reagent_string(prd_str)
    return rcts, prds


def _split_reagent_string(rgt_str):

    def _interpret_reagent_count(rgt_cnt_str):
        _pattern = (STRING_START + capturing(maybe(DIGIT)) +
                    capturing(one_or_more(NONSPACE)))
        cnt, rgt = find_first_capture(_pattern, rgt_cnt_str)
        cnt = int(cnt) if cnt else 1
        rgts = (rgt,) * cnt
        return rgts

    rgt_str = find_remove(LINESPACES, rgt_str)
    rgt_str = find_remove(CHEMKIN_PAREN_PLUS_EM, rgt_str)
    rgt_str = find_remove(CHEMKIN_PLUS_EM, rgt_str)
    pattern = PLUS + not_followed_by(PLUS)
    rgt_cnt_strs = find_split(pattern, rgt_str)
    rgts = tuple(chain(*map(_interpret_reagent_count, rgt_cnt_strs)))
    return rgts


def reaction_data_high_p_coeffs(rxn_dstr):
    """ get the high-pressure Arrhenius coefficients from a reaction data string
    """
    headline = find_split_lines(rxn_dstr)[0]
    number = one_of_these(
        [EXPONENTIAL_FLOAT, EXPONENTIAL_INTEGER, FLOAT, INTEGER])
    pattern = (LINESPACES + capturing(number) + LINESPACES +
               capturing(number) + LINESPACES + capturing(number))
    captures = find_first_capture(pattern, headline)
    assert captures
    cfts = tuple(map(float, captures))
    return cfts


def reaction_data_is_duplicate(rxn_dstr):
    """ is this a duplicate reaction?
    """
    pattern = 'DUP'
    return find_if_has_match(pattern, rxn_dstr)


def thermo_data(mech_str):
    """ find all thermo data
    """
    thm_dstr_lst = thermo_data_strings(mech_str)
    thm_dat_lst = tuple(zip(
        map(thermo_data_species_name, thm_dstr_lst),
        map(thermo_data_lo_coefficients, thm_dstr_lst),
        map(thermo_data_hi_coefficients, thm_dstr_lst),
        map(thermo_data_temperatures, thm_dstr_lst)))
    return thm_dat_lst


def thermo_data_strings(mech_str):
    """ find all thermo data strings
    """
    block_str = remove_blanks(thermo_block(mech_str))
    start_pattern = LINE_START + not_followed_by(
        one_of_these([DIGIT, PLUS, escape('=')]))
    end_pattern = '1' + LINE_END
    headline_pattern = start_pattern + one_or_more(NONNEWLINE) + end_pattern
    thm_dstr_lst = _headlined_sections(headline_pattern, block_str)
    assert all(len(find_split_lines(thm_dstr)) == 4
               for thm_dstr in thm_dstr_lst)
    return thm_dstr_lst


def thermo_data_species_name(thm_dstr):
    """ get the species name from a thermo data string
    """
    pattern = STRING_START + capturing(one_or_more(NONSPACE))
    spc = find_first_capture(pattern, thm_dstr)
    return spc


def thermo_data_temperatures(thm_dstr):
    """ get the common temperature from a thermo data string
    """
    headline = find_split_lines(thm_dstr)[0]
    pattern = (LINESPACES + capturing(UNSIGNED_FLOAT) +
               LINESPACES + capturing(UNSIGNED_FLOAT) +
               LINESPACES + capturing(UNSIGNED_FLOAT))
    captures = find_first_capture(pattern, headline)
    assert captures
    tmps = tuple(map(float, captures))
    return tmps


def thermo_data_lo_coefficients(thm_dstr):
    """ get the low temperature thermo coefficients
    """
    capture_lst = find_captures(EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[7:14]))
    return cfts


def thermo_data_hi_coefficients(thm_dstr):
    """ get the low temperature thermo coefficients
    """
    capture_lst = find_captures(EXPONENTIAL_FLOAT, thm_dstr)
    assert len(capture_lst) in (14, 15)
    cfts = tuple(map(float, capture_lst[:7]))
    return cfts


def reaction_unit_names(mech_str):
    """ units specified in the reaction block
    """
    block_str = remove_blanks(reactions_block(mech_str))
    a_unit_names, _ = zip(*A_UNITS)
    e_unit_names, _ = zip(*E_UNITS)
    a_pattern = (STRING_START +
                 maybe(one_of_these(e_unit_names) + LINESPACES) +
                 capturing(one_of_these(a_unit_names)))
    e_pattern = (STRING_START +
                 maybe(one_of_these(a_unit_names) + LINESPACES) +
                 capturing(one_of_these(e_unit_names)))
    a_unit_name = find_first_capture(a_pattern, block_str)
    e_unit_name = find_first_capture(e_pattern, block_str)
    return a_unit_name, e_unit_name


def thermo_t_common_default(mech_str):
    """ temperature defaults from the thermo block
    """
    block_str = remove_blanks(thermo_block(mech_str))
    pattern = (STRING_START +
               UNSIGNED_FLOAT + LINESPACES +
               capturing(UNSIGNED_FLOAT) + LINESPACES +
               UNSIGNED_FLOAT)
    capture = find_first_capture(pattern, block_str)
    assert capture
    tmp_com_def = float(capture)
    return tmp_com_def


def species_block(mech_str):
    """ find the species block
    """
    return _block(mech_str, block_keys=['SPECIES', 'SPEC'])


def reactions_block(mech_str):
    """ find the reactions block
    """
    return _block(mech_str, block_keys=['REACTIONS', 'REAC'])


def thermo_block(mech_str):
    """ find the thermodynamics block
    """
    return _block(mech_str, block_keys=['THERMO ALL', 'THERM ALL', 'THER ALL',
                                        'THERMO', 'THERM', 'THER'])


def remove_comments(mech_str):
    """ remove comments
    """
    pattern = escape('!') + zero_or_more(NONNEWLINE)
    clean_mech_str = find_remove(pattern, mech_str)
    assert clean_mech_str is not None
    return clean_mech_str


def remove_blanks(mech_str):
    """ remove blank lines as well as leading and trailing blanks
    """
    blank_line = LINE_START + maybe(LINESPACES) + NEWLINE
    trailing_blanks = LINESPACES + LINE_END
    leading_blanks = LINE_START + LINESPACES
    pattern = one_of_these([blank_line, trailing_blanks, leading_blanks])
    return find_remove(pattern, mech_str)


def _block(mech_str, block_keys):
    start_key = one_of_these(block_keys)
    contents = capturing(one_or_more(WILDCARD, greedy=False))
    pattern = start_key + contents + 'END'
    mech_str = remove_comments(mech_str)
    block = find_first_capture(pattern, mech_str)
    return find_strip_spaces(block) if block else None


def _headlined_sections(pattern, string):
    """ return sections with headlines matching a pattern
    """
    lines = string.splitlines()
    join_lines = '\n'.join
    pattern_matcher = find_matcher(pattern)
    lines = lstrip(lines, pred=lambda line: not pattern_matcher(line))
    sections = list(map(join_lines, split_before(lines, pattern_matcher)))
    return sections
