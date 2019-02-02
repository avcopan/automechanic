""" Test the automech CLI
"""
import os
import tempfile
import subprocess
from automechanic import fs

AUTOMECH_CMD = 'automech'
PATH = os.path.dirname(os.path.realpath(__file__))
HEPTANE_PATH = os.path.join(PATH, '../examples/heptane')


def test__help():
    """ test `automech -h`
    """
    subprocess.check_call([AUTOMECH_CMD, '-h'])


def test__chemkin__help():
    """ test `automech chemkin -h`
    """
    subprocess.check_call([AUTOMECH_CMD, 'chemkin', '-h'])


def test__chemkin__parse():
    """ test `automech chemkin parse`
    """
    subprocess.check_call([AUTOMECH_CMD, 'chemkin', 'parse', '-h'])

    tmp_dir = tempfile.mkdtemp()
    print(tmp_dir)

    with fs.enter(tmp_dir):
        mech_txt = os.path.join(HEPTANE_PATH, 'mechanism.txt')
        ther_txt = os.path.join(HEPTANE_PATH, 'thermo_data.txt')
        subprocess.check_call([AUTOMECH_CMD, 'chemkin', 'parse',
                               mech_txt, ther_txt, '-p'])


def test__species__help():
    """ test `automech species -h`
    """
    subprocess.check_call([AUTOMECH_CMD, 'species', '-h'])


def test__species__inchi():
    """ test `automech species inchi`
    """
    subprocess.check_call([AUTOMECH_CMD, 'species', 'inchi', '-h'])

    tmp_dir = tempfile.mkdtemp()
    print(tmp_dir)

    with fs.enter(tmp_dir):
        smi_csv = os.path.join(HEPTANE_PATH, 'smiles.csv')
        subprocess.check_call([AUTOMECH_CMD, 'species', 'inchi',
                               'smiles', smi_csv, '-p'])


def test__species__filesystem():
    """ test `automech species filesystem`
    """
    subprocess.check_call([AUTOMECH_CMD, 'species', 'filesystem', '-h'])

    tmp_dir = tempfile.mkdtemp()
    print(tmp_dir)

    with fs.enter(tmp_dir):
        smi_csv = os.path.join(HEPTANE_PATH, 'smiles.csv')
        subprocess.check_call([AUTOMECH_CMD, 'species', 'inchi',
                               'smiles', smi_csv, '-m', 'expand', '-p'])
        subprocess.check_call([AUTOMECH_CMD, 'species', 'filesystem',
                               'species_inchi.csv', '-F', 'fs', '-p'])


if __name__ == '__main__':
    # test__help()
    # test__chemkin__help()
    # test__chemkin__parse()
    # test__species__help()
    # test__species__inchi()
    test__species__filesystem()
