""" test the automechanic.form module
"""
from automechanic_old import form


def test__subtract():
    """ test form.subtract()
    """
    assert form.subtract({'H': 3, 'C': 2}, {'H': 4, 'C': 2}) == {'H': -1}
