""" CLI commands
"""
import os
from . import arglib as al
from .arg import specifier_from_kernel as specifier
from .clihelp import call_subcommand
from .clihelp import call_task
from .. import task
from .. import params as par

RXN_CSV_CHAR = 'r'
SPC_CSV_CHAR = 's'

HOME_DIR = os.path.expanduser("~")
FILESYSTEM_NAME = 'automech_fs'
FILESYSTEM_PREFIX_DEF = os.path.join(HOME_DIR, FILESYSTEM_NAME)
FILESYSTEM_PREFIX_CHAR = 'f'

STEREO_MODE_CHAR = 'm'


def automech(argt):
    """ automech command
    """
    call_subcommand(
        argt,
        subcmds=(
            (par.CHK_ENG, chemkin),
            (par.SPC_ENG, species),
            (par.RXN_ENG, reactions),
        )
    )


def chemkin(argt):
    """ chemkin sub-command
    """
    call_subcommand(
        argt,
        subcmds=(
            (par.PRS_ENG, chemkin__parse),
        )
    )


def chemkin__parse(argt):
    """ parse CHEMKIN mechanism
    """
    call_task(
        argt,
        task.chemkin.parse,
        specs=(
            specifier(
                al.MECHANISM_TXT, inp=True,
                extra_kwargs=(('nargs', '+'),),
                extra_helps=('followed by a CHEMKIN thermo file, if needed',),
                # value_map=(lambda x: tuple(map(os.path.abspath, x)))
                # ^ not necessary after deleting the PREFIX option
            ),
            specifier(
                al.REACTIONS_CSV, out=True, opt_char=RXN_CSV_CHAR.upper(),
                extra_kwargs=(
                    ('default', task.chemkin.DEFS.PARSE.RXN_CSV_DEF),),
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(
                    ('default', task.chemkin.DEFS.PARSE.SPC_CSV_DEF),),
            ),
        )
    )


def species(argt):
    """ species sub-command
    """
    call_subcommand(
        argt,
        subcmds=(
            (par.ICH_ENG, species__inchi),
            (par.FLS_ENG, species__filesystem),
        )
    )


def species__inchi(argt):
    """ expand species into their possible stereoisomers
    """
    call_task(
        argt,
        task.species.inchi,
        specs=(
            specifier(
                al.SPECIES_ID,
                allowed_values=task.species.VALS.INCHI.INP_ID_KEY,
            ),
            specifier(
                al.SPECIES_CSV, inp=True,
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(
                    ('default', task.species.DEFS.INCHI.SPC_CSV),),
            ),
            specifier(
                al.STEREO_MODE, opt_char=STEREO_MODE_CHAR,
                allowed_values=task.species.VALS.INCHI.STEREO_MODE,
                extra_kwargs=(
                    ('default', task.species.DEFS.INCHI.STEREO_MODE),)
            ),
        )
    )


def species__filesystem(argt):
    """ create the species filesystem
    """
    call_task(
        argt,
        task.species.filesystem,
        specs=(
            specifier(
                al.SPECIES_CSV, inp=True,
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(
                    ('default', task.species.DEFS.FILESYSTEM.SPC_CSV),),
            ),
            specifier(
                al.FILESYSTEM_PREFIX, out=True,
                opt_char=FILESYSTEM_PREFIX_CHAR.upper(),
                extra_kwargs=(('default', FILESYSTEM_PREFIX_DEF),),
            ),
        )
    )


def reactions(argt):
    """ reactions sub-command
    """
    call_subcommand(
        argt,
        subcmds=(
            ('filesystem', reactions__filesystem),
        )
    )


def reactions__filesystem(argt):
    """ create the reactions filesystem
    """
    raise NotImplementedError(argt)
