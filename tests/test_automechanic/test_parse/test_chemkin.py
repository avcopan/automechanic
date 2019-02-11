""" test the automechanic.parse.chemkin module
"""
from __future__ import unicode_literals
from builtins import open
import os
from automechanic.parse import chemkin

PATH = os.path.dirname(os.path.realpath(__file__))
NATGAS_PATH = os.path.join(PATH, '../../../examples/natgas')
HEPTANE_PATH = os.path.join(PATH, '../../../examples/heptane')


def test__reaction_data():
    """ test chemkin.reaction_data_strings
    """
    mech_txt = os.path.join(HEPTANE_PATH, 'mechanism.txt')
    mech_str = open(mech_txt, encoding='utf8', errors='ignore').read()
    rxn_dat_lst = chemkin.reaction_data(mech_str)
    assert len(rxn_dat_lst) == 5336


def test__thermo_data():
    """ test chemkin.thermo_data
    """
    ther_txt = os.path.join(HEPTANE_PATH, 'thermo_data.txt')
    ther_str = open(ther_txt, encoding='utf8', errors='ignore').read()
    thm_dat_lst = chemkin.thermo_data(ther_str)
    assert len(thm_dat_lst) == 1268


def test__reaction_unit_names():
    """ test chemkin.reaction_unit_names()
    """
    mech_txt = os.path.join(HEPTANE_PATH, 'mechanism.txt')
    mech_str = open(mech_txt, encoding='utf8', errors='ignore').read()
    unts = chemkin.reaction_unit_names(mech_str)
    assert unts == (None, None)

    mech_txt = os.path.join(NATGAS_PATH, 'mechanism.txt')
    mech_str = open(mech_txt, encoding='utf8', errors='ignore').read()
    unts = chemkin.reaction_unit_names(mech_str)
    assert unts == ('MOLES', 'KCAL/MOLE')


def test__thermo_t_common_default():
    """ test chemkin.thermo_t_common_default()
    """
    ther_txt = os.path.join(HEPTANE_PATH, 'thermo_data.txt')
    ther_str = open(ther_txt, encoding='utf8', errors='ignore').read()
    tmp_com_def = chemkin.thermo_t_common_default(ther_str)
    assert tmp_com_def == 1000.

    mech_txt = os.path.join(NATGAS_PATH, 'mechanism.txt')
    mech_str = open(mech_txt, encoding='utf8', errors='ignore').read()
    tmp_com_def = chemkin.thermo_t_common_default(mech_str)
    assert tmp_com_def == 1000.


if __name__ == '__main__':
    # test__thermo_t_common_default()
    test__thermo_data()
    # test__reaction_unit_names()
    # test__reaction_data()
