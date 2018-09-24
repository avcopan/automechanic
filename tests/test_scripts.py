""" Test scripts
"""
import os
import subprocess
import tempfile

PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(PATH, '../examples/templates')
SYNGAS_PATH = os.path.join(PATH, '../examples/syngas')


def test__automech__syngas_from_rmg():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    print tmp_path

    mech_json = os.path.join(SYNGAS_PATH, 'mechanism.json')
    spc_json = os.path.join(SYNGAS_PATH, 'species.json')
    subprocess.check_call(['automech', 'init_from_rmg',
                           mech_json,
                           spc_json,
                           '-P', tmp_path,
                           '-p'])
    subprocess.check_call(['automech', 'abstractions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p'])
    subprocess.check_call(['automech', 'migrations', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           '-y', 'nodes:d',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'migrations', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])

    # test divided runs
    subprocess.check_call(['automech', 'divide',
                           'rad-rad', 'rad-rad', 'rad-mol',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'divide',
                           'high-spin', 'high-spin', 'low-spin',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad'),
                           '-p'])
    subprocess.check_call(['automech', 'abstractions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol', 'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'abstractions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])


def test__automech__syngas_from_chemkin():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    print tmp_path

    subprocess.check_call(['automech', 'init',
                           os.path.join(SYNGAS_PATH, 'mechanism.txt'),
                           os.path.join(SYNGAS_PATH, 'species.csv'),
                           '-P', tmp_path])
    subprocess.check_call(['automech', 'abstractions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p'])
    subprocess.check_call(['automech', 'migrations', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p'])
    subprocess.check_call(['automech', 'additions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'additions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           '-y', 'nodes:d',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'migrations', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'migrations',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])

    # test divided runs
    subprocess.check_call(['automech', 'divide',
                           'rad-rad', 'rad-rad', 'rad-mol',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech', 'divide',
                           'high-spin', 'high-spin', 'low-spin',
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad'),
                           '-p'])
    subprocess.check_call(['automech', 'abstractions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol', 'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol', 'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-mol'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech', 'abstractions', 'run_batch',
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-s', os.path.join(tmp_path, 'species.csv'),
                           '-r', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin',
                                              'reactions.csv'),
                           '-b', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin',
                                              'reactions.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions',
                                              'rad-rad', 'high-spin'),
                           '-y', 'nodes:d',
                           '-p',
                           'cmd', 'ls'])


if __name__ == '__main__':
    # test__automech__syngas_from_rmg()
    test__automech__syngas_from_chemkin()
