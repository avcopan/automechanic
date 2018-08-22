""" functions and regexes for parsing CHEMKIN files
"""
import re
from functools import partial
from .parse import WHITESPACES
from .parse import LINE
from .parse import PLUS
from .parse import OPEN_PAREN
from .parse import CLOSE_PAREN
from .parse import escape
from .parse import maybe
from .parse import one_of_these
from .parse import zero_or_more
from .parse import repeat_range
from .parse import capture
from .parse import named_capture


ARROW = maybe(escape('<')) + escape('=') + maybe(escape('>'))
PADDED_PLUS = maybe(WHITESPACES) + PLUS + maybe(WHITESPACES)
PADDED_ARROW = maybe(WHITESPACES) + ARROW + maybe(WHITESPACES)


def remove_comments(mech_str):
    """ remove comments from a CHEMKIN filestring

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: commentr-free CHEMKIN filestring
    :rtype: str
    """
    comment = '!.*'
    clean_mech_str = re.sub(comment, '', mech_str)
    return clean_mech_str


def find_species(mech_str):
    """ find species in a CHEMKIN file
    """
    species_block = _find_chemkin_block(mech_str, 'SPECIES')
    species = species_block.split()
    return species


def find_reactions(mech_str):
    """ find reactions in a CHEMKIN file
    """
    species = _long_to_short(map(escape, find_species(mech_str)))
    pi_rxn_pattern = _pressure_independent_reaction_pattern(species)
    lp_rxn_pattern = _low_pressure_reaction_pattern(species)
    fo_rxn_pattern = _falloff_reaction_pattern(species)
    split = partial(_split_reaction, species=species)

    reactions_block = _find_chemkin_block(mech_str, 'REACTIONS')
    pi_rxns = list(map(split, re.findall(pi_rxn_pattern, reactions_block)))
    lp_rxns = list(map(split, re.findall(lp_rxn_pattern, reactions_block)))
    fo_rxns = list(map(split, re.findall(fo_rxn_pattern, reactions_block)))
    return pi_rxns, lp_rxns, fo_rxns


def reaction_to_smirks(rxn, smi_dct):
    """ convert a CHEMKIN reaction to the SMIRKS format
    """
    rcts, prds = rxn
    rsmis = [smi_dct[spc] for spc in rcts]
    psmis = [smi_dct[spc] for spc in prds]
    smrk = '.'.join(rsmis) + '>>' + '.'.join(psmis)
    return smrk


def _split_reaction(rxn, species):
    reagent_capture = capture(one_of_these(species))
    rct_str, prd_str = re.split(PADDED_ARROW, rxn)
    rcts = re.findall(reagent_capture, rct_str)
    prds = re.findall(reagent_capture, prd_str)
    return rcts, prds


def _long_to_short(iterable):
    return list(reversed(sorted(iterable, key=len)))


def _pressure_independent_reaction_pattern(species):
    reagent = one_of_these(species)
    reagents = reagent + repeat_range(PADDED_PLUS + reagent, 0, 2)
    reaction = reagents + PADDED_ARROW + reagents
    return reaction


def _low_pressure_reaction_pattern(species):
    reagent = one_of_these(species)
    reagents = reagent + maybe(PADDED_PLUS + reagent) + PADDED_PLUS + 'M'
    reaction = reagents + PADDED_ARROW + reagents
    return reaction


def _falloff_reaction_pattern(species):
    reagent = one_of_these(species)
    reagents = (reagent + maybe(PADDED_PLUS + reagent)
                + OPEN_PAREN + PADDED_PLUS + 'M' + CLOSE_PAREN)
    reaction = reagents + PADDED_ARROW + reagents
    return reaction


def _find_chemkin_block(mech_str, name):
    lines_nongreedy = zero_or_more(LINE, greedy=False)
    head = r'^\s*{:s}\s*$'.format(name.upper())
    body = named_capture(lines_nongreedy, name='block')
    foot = r'^\s*END\s*$'
    block = head + body + foot

    clean_mech_str = remove_comments(mech_str)
    match = re.search(block, clean_mech_str, re.MULTILINE)
    assert match
    gdct = match.groupdict()
    block = gdct['block']
    return block