""" automechanic run parameters
"""
from . import tab

# some english names for things
CHK_ENG = 'chemkin'
SPC_ENG = 'species'
RXN_ENG = 'reactions'

PRS_ENG = 'parse'
ICH_ENG = 'inchi'
FLS_ENG = 'filesystem'

# extensions
CSV_EXT = 'csv'

# common schema parameters
_NAME_KEY = 'name'
_NAME_TYP = tab.dt_(str)
_INP_ID_ICH_KEY = ICH_ENG
_INP_ID_SMI_KEY = 'smiles'
_INP_ID_TYP = tab.dt_(str)
_ICH_KEY = ICH_ENG
_ICH_TYP = tab.dt_(str)
_MULT_KEY = 'mult'
_MULT_TYP = tab.dt_(int)
_FILESYSTEM_PATH_KEY = 'path'
_FILESYSTEM_PATH_TYP = tab.dt_(str)
_STRING_PLACEHOLDER = '{:s}'


class SPC():
    """ species parameters
    """
    ENGLISH_NAME = SPC_ENG

    INP_ID_ICH_KEY = _INP_ID_ICH_KEY
    INP_ID_SMI_KEY = _INP_ID_SMI_KEY
    MULT_KEY = _MULT_KEY

    ICH_KEY = _ICH_KEY

    FILESYSTEM_DIR_NAME = 'SPC'

    PICK_STEREO = 'pick'
    EXPAND_STEREO = 'expand'

    class IO():
        """ file and directory names for CLI input/output
        """
        _CSV_WITH_KEYWORD = '{:s}_{:s}.{:s}'.format(
            SPC_ENG, _STRING_PLACEHOLDER, CSV_EXT)

        # name only
        CHEMKIN_CSV = _CSV_WITH_KEYWORD.format(CHK_ENG)
        # name and inchi
        INCHI_CSV = _CSV_WITH_KEYWORD.format(ICH_ENG)
        # name and inchi and filesystem path
        FILESYSTEM_CSV = _CSV_WITH_KEYWORD.format(FLS_ENG)

    class TAB():
        """ species table parameters
        """
        NAME_KEY = _NAME_KEY
        NAME_TYP = _NAME_TYP

        INP_ID_TYP = _INP_ID_TYP
        ICH_TYP = _ICH_TYP

        MULT_TYP = _MULT_TYP

        FILESYSTEM_PATH_KEY = _FILESYSTEM_PATH_KEY
        FILESYSTEM_PATH_TYP = _FILESYSTEM_PATH_TYP

        NASA_C_TYP = tab.dt_(float)
        NASA_C_LO_KEYS = ('nasa_lo_1', 'nasa_lo_2', 'nasa_lo_3', 'nasa_lo_4',
                          'nasa_lo_5', 'nasa_lo_6', 'nasa_lo_7')
        NASA_C_HI_KEYS = ('nasa_hi_1', 'nasa_hi_2', 'nasa_hi_3', 'nasa_hi_4',
                          'nasa_hi_5', 'nasa_hi_6', 'nasa_hi_7')
        NASA_T_TYP = tab.dt_(float)
        NASA_T_KEYS = ('t_lo', 't_hi', 't_c')


class RXN():
    """ reaction parameters
    """
    ENGLISH_NAME = RXN_ENG

    MULT_KEY = _MULT_KEY

    FILESYSTEM_DIR_NAME = 'RXN'

    class IO():
        """ file and directory names for CLI input/output
        """
        _CSV_WITH_KEYWORD = '{:s}_{:s}.{:s}'.format(
            RXN_ENG, _STRING_PLACEHOLDER, CSV_EXT)

        # name only
        CHEMKIN_CSV = _CSV_WITH_KEYWORD.format(CHK_ENG)

    class TAB():
        """ species table parameters
        """
        NAME_KEY = _NAME_KEY
        NAME_TYP = _NAME_TYP

        INP_ID_TYP = _INP_ID_TYP

        MULT_TYP = _MULT_TYP

        FILESYSTEM_PATH_KEY = _FILESYSTEM_PATH_KEY
        FILESYSTEM_PATH_TYP = _FILESYSTEM_PATH_TYP

        ARRH_TYP = tab.dt_(float)
        ARRH_KEYS = ('arrh_a', 'arrh_b', 'arrh_e')
