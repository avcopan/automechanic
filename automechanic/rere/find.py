""" re finders
"""
import re
from functools import partial
from more_itertools import split_before
from more_itertools import lstrip
from .pattern import one_or_more
from .pattern_lib import STRING_START
from .pattern_lib import STRING_END
from .pattern_lib import NEWLINE
from .pattern_lib import WHITESPACE

WHITESPACES = one_or_more(WHITESPACE)


def has_match(pattern, string):
    """ does this string have a pattern match?
    """
    match = re.search(pattern, string, flags=re.MULTILINE)
    return bool(match)


def ends_with(pattern, string):
    """ does the string end with this pattern
    """
    end_pattern = pattern + STRING_END
    return has_match(end_pattern, string)


def matcher(pattern):
    """ return a boolean matching function
    """
    return partial(has_match, pattern)


def all_captures(pattern, string):
    """ capture(s) for all matches of a capturing pattern
    """
    return tuple(re.findall(pattern, string, flags=re.MULTILINE))


def first_capture(pattern, string):
    """ capture(s) from first match for a capturing pattern
    """
    match = re.search(pattern, string, flags=re.MULTILINE)
    return (match.group(1) if match and len(match.groups()) == 1 else
            match.groups() if match else None)


def last_capture(pattern, string):
    """ capture(s) from first match for a capturing pattern
    """
    caps_lst = all_captures(pattern, string)
    return caps_lst[-1] if caps_lst else None


def first_named_capture(pattern, string):
    """ capture dictionary from first match for a pattern with named captures
    """
    match = re.search(pattern, string, flags=re.MULTILINE)
    return match.groupdict() if match and match.groupdict() else None


def split(pattern, string):
    """ split string at matches
    """
    return tuple(re.split(pattern, string, maxsplit=0, flags=re.MULTILINE))


def split_words(string):
    """ split string at whitespaces
    """
    return split(WHITESPACES, strip_spaces(string))


def split_lines(string):
    """ split string at newlines
    """
    return split(NEWLINE, string)


def remove(pattern, string):
    """ remove pattern matches
    """
    return replace(pattern, '', string)


def strip_spaces(string):
    """ strip spaces from the string ends
    """
    lspaces = STRING_START + WHITESPACES
    rspaces = WHITESPACES + STRING_END
    return remove(lspaces, remove(rspaces, string))


def replace(pattern, repl, string):
    """ replace pattern matches
    """
    return re.sub(pattern, repl, string, count=0, flags=re.MULTILINE)


def headlined_sections(pattern, string):
    """ return sections with headlines matching a pattern
    """
    lines = string.splitlines()
    join_lines = '\n'.join
    pattern_matcher = matcher(pattern)
    lines = lstrip(lines, pred=lambda line: not pattern_matcher(line))
    sections = list(map(join_lines, split_before(lines, pattern_matcher)))
    return sections
