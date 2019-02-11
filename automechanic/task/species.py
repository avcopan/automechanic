""" tasks operating on species
"""
import automol as mol
from .. import params as par
from .. import tab
from .. import fslib
from .. import fs
from ._util import read_csv as _read_csv
from ._util import write_csv as _write_csv


class VALS():
    """ function argument values """
    class INCHI():
        """_"""
        INP_ID_KEY = (par.SPC.INP_ID_SMI_KEY, par.SPC.INP_ID_ICH_KEY)
        STEREO_MODE = (par.SPC.EXPAND_STEREO, par.SPC.PICK_STEREO)

    class FILESYSTEM():
        """_"""


class DEFS():
    """ function argument defaults"""
    class INCHI():
        """_"""
        SPC_CSV = par.SPC.IO.INCHI_CSV
        STEREO_MODE = par.SPC.PICK_STEREO

    class FILESYSTEM():
        """_"""
        SPC_CSV = par.SPC.IO.FILESYSTEM_CSV
        # RXN_CSV = par.RXN.IO.FILESYSTEM_CSV


def inchi(inp_id_key, spc_csv, spc_csv_out, stereo_mode, logger):
    """ convert species identifiers to InChI
    """
    assert inp_id_key in VALS.INCHI.INP_ID_KEY
    assert stereo_mode in VALS.INCHI.STEREO_MODE

    tbl = _read_csv(spc_csv, logger)
    tbl = _inchi_table(tbl, inp_id_key, logger)
    tbl = _assign_stereo(tbl, mode=stereo_mode, logger=logger)
    _write_csv(spc_csv_out, tbl, logger)


def filesystem(spc_csv, spc_csv_out, filesystem_prefix, logger):
    """ chart the species filesystem structure
    """

    tbl = _read_csv(spc_csv, logger)
    tbl = _create_filesystem(tbl, fs_root_pth=filesystem_prefix, logger=logger)
    _write_csv(spc_csv_out, tbl, logger)


def _inchi_table(tbl, inp_id_key, logger):
    assert inp_id_key in (par.SPC.INP_ID_SMI_KEY, par.SPC.INP_ID_ICH_KEY)
    ich_key = par.SPC.ICH_KEY
    action = ("Calculating" if inp_id_key == ich_key else "Recalculating")
    logger.info("{:s} '{:s}' from '{:s}'".format(action, ich_key, inp_id_key))

    tbl = tab.enforce_schema(tbl,
                             keys=(inp_id_key,),
                             typs=(par.SPC.TAB.INP_ID_TYP,))

    conv_ = (mol.inchi.recalculate if inp_id_key == par.SPC.INP_ID_ICH_KEY else
             mol.smiles.inchi)

    sids = tbl[inp_id_key]
    ichs = list(map(conv_, sids))
    tbl = tbl[[key for key in tab.keys_(tbl) if key != inp_id_key]]
    tbl[par.SPC.INP_ID_ICH_KEY] = ichs
    return tbl


def _assign_stereo(tbl, mode, logger):
    assert mode in (par.SPC.EXPAND_STEREO, par.SPC.PICK_STEREO)
    logger.info("Assigning stereo in mode '{:s}'".format(mode))

    tbl = (_assign_stereo_by_expanding(tbl) if mode == par.SPC.EXPAND_STEREO
           else _assign_stereo_by_picking(tbl))
    return tbl


def _assign_stereo_by_picking(tbl):
    tbl = tab.enforce_schema(tbl,
                             keys=(par.SPC.INP_ID_ICH_KEY,),
                             typs=(par.SPC.TAB.INP_ID_TYP,))
    tbl = tbl.copy()
    # use coordinates to get stereo assignments
    ichs = tbl[par.SPC.INP_ID_ICH_KEY]
    ichs = list(map(mol.geom.stereo_inchi, map(mol.inchi.geometry, ichs)))
    assert not any(map(mol.inchi.has_unknown_stereo_elements, ichs))
    tbl[par.SPC.INP_ID_ICH_KEY] = ichs
    return tbl


def _assign_stereo_by_expanding(tbl):
    tbl = tab.enforce_schema(tbl,
                             keys=(par.SPC.INP_ID_ICH_KEY,),
                             typs=(par.SPC.TAB.INP_ID_TYP,))

    ichs = tbl[par.SPC.INP_ID_ICH_KEY]
    ichsts_lst = list(map(mol.inchi.compatible_stereoisomers, ichs))

    idx_save_key = tab.next_index_save_key(tbl)
    tbl = tab.save_index(tbl)

    tbl_idxs = tab.idxs_(tbl)
    vals = [[idx, ichst]
            for idx, ichsts in zip(tbl_idxs, ichsts_lst) for ichst in ichsts]
    ste_tbl = tab.from_records(vals,
                               keys=(idx_save_key, par.SPC.INP_ID_ICH_KEY),
                               typs=(tab.IDX_TYP, par.SPC.TAB.INP_ID_TYP))

    keys = [key for key in tab.keys_(tbl) if key != par.SPC.INP_ID_ICH_KEY]
    tbl = tab.left_join(tbl[keys], ste_tbl, key=idx_save_key)
    return tbl


def _create_filesystem(tbl, fs_root_pth, logger):
    logger.info("Creating filesystem at '{:s}'".format(fs_root_pth))

    id_keys = (par.SPC.TAB.ICH_KEY, par.SPC.TAB.MULT_KEY)
    id_typs = (par.SPC.TAB.ICH_TYP, par.SPC.TAB.MULT_TYP)
    tbl = tab.enforce_schema(tbl, keys=id_keys, typs=id_typs)

    def __create_branch(ich, mult):
        sgms = fslib.species.branch_segments(ich, mult)
        return fs.branch.create(sgms)

    with fs.enter(fs_root_pth):
        pth_tbl = tab.from_starmap(tbl, __create_branch, id_keys,
                                   keys=(par.SPC.TAB.FILESYSTEM_PATH_KEY,),
                                   typs=(par.SPC.TAB.FILESYSTEM_PATH_TYP,))
        tbl = tab.left_join(tbl, pth_tbl)

    return tbl


# used elswhere -- decide if this is the best place
def grouped_identifier_lists(tbl, logger):
    """ species identifier lists, grouped by name
    """
    logger.info("Determining species identifiers, by name")

    keys = (par.SPC.TAB.NAME_KEY, par.SPC.TAB.ICH_KEY, par.SPC.TAB.MULT_KEY)
    typs = (par.SPC.TAB.NAME_TYP, par.SPC.TAB.ICH_TYP, par.SPC.TAB.MULT_TYP)
    tbl = tab.enforce_schema(tbl, keys=keys, typs=typs)

    key1 = par.SPC.TAB.NAME_KEY
    key2 = (par.SPC.TAB.ICH_KEY, par.SPC.TAB.MULT_KEY)
    ids_dct = tab.group_dictionary(tbl, key1, key2)
    return ids_dct
