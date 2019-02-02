""" tasks that operate on chemkin-format files
"""
from .. import params as par
from .. import tab
from ..iohelp import read_string
from ..iohelp import timestamp_if_exists
from ..parse.chemkin import species_names
from ..parse.chemkin import thermo_data
from ..parse.chemkin import reaction_data


class DEFS():
    """ function argument defaults"""
    class PARSE():
        """_"""
        SPC_CSV_DEF = par.SPC.IO.CHEMKIN_CSV
        RXN_CSV_DEF = par.RXN.IO.CHEMKIN_CSV


def parse(mech_txt_lst, rxn_csv_out, spc_csv_out, logger):
    """ parse chemkin information to CSV
    """
    mech_str = _read_mechanism_string(mech_txt_lst, logger)

    spc_eng = par.SPC.ENGLISH_NAME
    rxn_eng = par.RXN.ENGLISH_NAME

    spc_tbl = _create_species_table(mech_str, logger)
    _write_csv_file(spc_eng, spc_csv_out, spc_tbl, logger)

    rxn_tbl = _create_reactions_table(mech_str, logger)
    _write_csv_file(rxn_eng, rxn_csv_out, rxn_tbl, logger)


def _read_mechanism_string(mech_txt_lst, logger):
    logger.info("Reading in chemkin mechanism file(s)")

    mech_str = '\n'.join(map(read_string, mech_txt_lst))
    return mech_str


def _write_csv_file(eng, csv, tbl, logger):
    logger.info("Writing {:s} data to {:s}".format(eng, csv))

    timestamp_if_exists(csv)
    tab.write_csv(csv, tbl, float_format='%.6e')


def _create_species_table(mech_str, logger):
    eng = par.SPC.ENGLISH_NAME
    logger.info("Finding {:s} data in chemkin mechanism string".format(eng))

    spcs = species_names(mech_str)
    thm_dat_lst = thermo_data(mech_str)
    assert len(thm_dat_lst) == len(spcs)
    keys = (par.SPC.TAB.NAME_KEY, par.SPC.TAB.NASA_C_LO_KEYS,
            par.SPC.TAB.NASA_C_HI_KEYS, par.SPC.TAB.NASA_T_KEYS)
    typs = (par.SPC.TAB.NAME_TYP, par.SPC.TAB.NASA_C_TYP,
            par.SPC.TAB.NASA_C_TYP, par.SPC.TAB.NASA_T_TYP)
    spc_tbl = tab.from_records(vals=thm_dat_lst, keys=keys, typs=typs)
    return spc_tbl


def _create_reactions_table(mech_str, logger):
    eng = par.RXN.ENGLISH_NAME
    logger.info("Finding {:s} data in chemkin mechanism string".format(eng))

    rxn_dat_lst = reaction_data(mech_str)
    keys = (par.RXN.TAB.NAME_KEY, par.RXN.TAB.ARRH_KEYS)
    typs = (par.RXN.TAB.NAME_TYP, par.RXN.TAB.ARRH_TYP)
    rxn_tbl = tab.from_records(vals=rxn_dat_lst, keys=keys, typs=typs)
    return rxn_tbl
