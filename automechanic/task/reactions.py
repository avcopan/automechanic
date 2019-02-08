""" tasks operation on reactions
"""
# from ..parse.chemkin import split_reaction_name
from .species import grouped_identifier_lists as _grouped_identifier_lists
from ._util import read_csv as _read_csv
# from ._util import write_csv as _write_csv


def classify(rxn_csv, spc_csv, logger):
    """ classify reactions by type
    """
    spc_tbl = _read_csv(spc_csv, logger)
    rxn_tbl = _read_csv(rxn_csv, logger)

    logger.info(rxn_tbl)
    spc_ids_dct = _grouped_identifier_lists(spc_tbl, logger)
    print(spc_ids_dct)
