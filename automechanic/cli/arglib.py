""" library of command-line arguments for automech
"""

# common subcommand arg(s)
SUBCMD = (
    'subcommand',
    (
        ('type', str),
    )
)

LOG_NAME = (
    'log_name',
    (
        ('type', str),
        ('help', "log file for automech output")
    )
)

LOG_LEVEL = (
    'log_level',
    (
        ('type', int),
        ('default', 20),
        ('help', "log verbosity: 20=INFO, 10=DEBUG")
    )
)

PRINT_OUT = (
    'print_out',
    (
        ('action', 'store_true'),
        ('help', "print log output to screen")
    )
)

FILESYSTEM_PREFIX = (
    'filesystem_prefix',
    (
        ('type', str),
        ('help', "Path to automech filesystem"),
    )
)

REACTIONS_CSV = (
    'reactions_csv',
    (
        ('type', str),
        ('help', "CSV with reaction data"),
    )
)

SPECIES_CSV = (
    'species_csv',
    (
        ('type', str),
        ('help', "CSV with species data")
    )
)

# chemkin arg(s)
MECHANISM_TXT = (
    'mechanism_txt',
    (
        ('type', str),
        ('help', "CHEMKIN mechanism file")
    )
)

# species arg(s)
SPECIES_ID = (
    'species_id',
    (
        ('type', str),
        ('help', "Species identifier type")
    )
)

STEREO_MODE = (
    'stereo_mode',
    (
        ('type', str),
        ('help', "How to handle stereoisomers")
    )
)
