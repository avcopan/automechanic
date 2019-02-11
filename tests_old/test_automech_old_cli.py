""" Test scripts
"""
import os
import sys
import subprocess
import tempfile

PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(PATH, '../examples/old/templates')
SYNGAS_PATH = os.path.join(PATH, '../examples/old/syngas')


def test__automech__syngas_from_rmg():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    sys.stdout.write("{:s}\n".format(tmp_path))

    mech_json = os.path.join(SYNGAS_PATH, 'mechanism.json')
    spc_json = os.path.join(SYNGAS_PATH, 'species.json')
    subprocess.check_call(['automech_old', 'init_from_rmg',
                           mech_json,
                           spc_json,
                           '-P', tmp_path,
                           '-p'])
    _main_routine(tmp_path)


def test__automech__syngas_from_chemkin():
    """ test the automech script
    """
    tmp_path = tempfile.mkdtemp()

    sys.stdout.write("{:s}\n".format(tmp_path))

    subprocess.check_call(['automech_old', 'init',
                           os.path.join(SYNGAS_PATH, 'mechanism.txt'),
                           os.path.join(SYNGAS_PATH, 'species.csv'),
                           '-P', tmp_path])
    _main_routine(tmp_path)


def _main_routine(tmp_path):
    subprocess.check_call(['automech_old', 'abstractions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p'])
    subprocess.check_call(['automech_old', 'additions', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p'])
    subprocess.check_call(['automech_old', 'migrations', 'init',
                           os.path.join(tmp_path, 'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p'])
    subprocess.check_call(['automech_old', 'additions', 'setup_run',
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech_old', 'additions', 'run',
                           os.path.join(tmp_path, 'additions',
                                        'reactions.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'addition.txt'),
                           '-n', 'b444', 'b445', 'b447', 'b448', 'b449',
                           '-x', '1-2', '3', '4-20',
                           '-P', os.path.join(tmp_path, 'additions'),
                           '-p',
                           'cmd', 'sleep', '1'])
    subprocess.check_call(['automech_old', 'migrations', 'setup_run',
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech_old', 'migrations', 'run',
                           os.path.join(tmp_path, 'migrations',
                                        'reactions.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'migration.txt'),
                           '-n', 'b444', 'b445', 'b447', 'b448', 'b449',
                           '-P', os.path.join(tmp_path, 'migrations'),
                           '-p',
                           'cmd', 'sleep', '1'])
    subprocess.check_call(['automech_old', 'abstractions', 'setup_run',
                           os.path.join(tmp_path, 'abstractions',
                                        'reactions.csv'),
                           os.path.join(tmp_path, 'species.csv'),
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p',
                           'cmd', 'ls'])
    subprocess.check_call(['automech_old', 'abstractions', 'run',
                           os.path.join(tmp_path, 'abstractions',
                                        'reactions.csv'),
                           '-t', os.path.join(TEMPLATE_PATH,
                                              'abstraction.txt'),
                           '-n', 'b444', 'b445', 'b447', 'b448', 'b449',
                           '-P', os.path.join(tmp_path, 'abstractions'),
                           '-p',
                           'cmd', 'sleep', '1'])


if __name__ == '__main__':
    test__automech__syngas_from_rmg()
    test__automech__syngas_from_chemkin()
