""" some common task utilities
"""
from .. import tab
from ..iohelp import timestamp_if_exists


def read_csv(csv, logger):
    """ read a table from a CSV file
    """
    logger.info("Reading in {:s}".format(csv))

    tbl = tab.read_csv(csv)
    return tbl


def write_csv(csv, tbl, logger):
    """ write a table to a CSV file
    """
    logger.info("Writing to {:s}".format(csv))

    timestamp_if_exists(csv)
    tab.write_csv(csv, tbl)
