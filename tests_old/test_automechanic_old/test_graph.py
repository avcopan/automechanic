""" test the automechanic.graph module
"""
from itertools import permutations
import numpy.random
from automechanic_old import graph


def test__atoms():
    """ test graph.atoms()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert graph.atoms(mgrph) == ('C', 'C', 'H', 'H', 'H', 'H')


def test__bonds():
    """ test graph.bonds()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.bonds(mgrph1)
            == frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                          (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                          (frozenset([0, 1]), 2)]))
    assert graph.bonds(mgrph2) == frozenset()


def test__bond_keys():
    """ test graph.bond_keys()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (set(graph.bond_keys(mgrph1))
            == set([frozenset([0, 3]), frozenset([1, 4]), frozenset([0, 2]),
                    frozenset([0, 1]), frozenset([1, 5])]))
    assert graph.bond_keys(mgrph2) == ()


def test__bond_orders():
    """ test graph.bond_orders()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.bond_orders(mgrph1)
            == {frozenset([1, 5]): 1, frozenset([1, 4]): 1,
                frozenset([0, 2]): 1, frozenset([0, 3]): 1,
                frozenset([0, 1]): 2})
    assert graph.bond_orders(mgrph2) == {}


def test__subgraph():
    """ test graph.subgraph()
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H'),
             frozenset([(frozenset([5, 6]), 1), (frozenset([4, 7]), 1),
                        (frozenset([8, 4]), 1), (frozenset([0, 5]), 1),
                        (frozenset([1, 5]), 1), (frozenset([2, 5]), 1),
                        (frozenset([2, 3]), 1), (frozenset([2, 4]), 2)]))
    ref_mgrph = (('C', 'C', 'C'),
                 frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 2)]))
    assert graph.subgraph(mgrph, (2, 4, 5)) == ref_mgrph
    assert graph.subgraph(mgrph, (0,)) == (('H',), frozenset([]))
    assert graph.subgraph(mgrph, ()) == ((), frozenset([]))


def test__heavy_atom_subgraph():
    """ test graph.heavy_atom_subgraph()
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H'),
             frozenset([(frozenset([5, 6]), 1), (frozenset([4, 7]), 1),
                        (frozenset([8, 4]), 1), (frozenset([0, 5]), 1),
                        (frozenset([1, 5]), 1), (frozenset([2, 5]), 1),
                        (frozenset([2, 3]), 1), (frozenset([2, 4]), 2)]))
    ref_mgrph = (('C', 'C', 'C'),
                 frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 2)]))
    assert graph.heavy_atom_subgraph(mgrph) == ref_mgrph


def test__skeleton_graph():
    """ test graph.skeleton_graph()
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([10, 11]), 1), (frozenset([5, 6]), 1),
                        (frozenset([4, 7]), 1), (frozenset([8, 4]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 5]), 1),
                        (frozenset([2, 5]), 1), (frozenset([2, 3]), 1),
                        (frozenset([2, 4]), 2)]))
    ref_mgrph = ((('C', 1), ('C', 2), ('C', 3), ('H', 0), ('H', 1)),
                 frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 2)]))
    assert graph.skeleton_graph(mgrph) == ref_mgrph


def test__radical_sites():
    """ test graph.radical_sites()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    mgrph3 = (('C', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([8, 3]), 1),
                         (frozenset([4, 6]), 1), (frozenset([4, 5]), 1),
                         (frozenset([0, 2]), 1), (frozenset([3, 4]), 1),
                         (frozenset([3, 7]), 1), (frozenset([0, 1]), 1)]))
    assert graph.radical_sites(mgrph1) == ()
    assert graph.radical_sites(mgrph2) == (0,)
    assert graph.radical_sites(mgrph3) == (0, 4)


def test__heavy_atom_indices():
    """ test graph.heavy_atom_indices
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H'),
             frozenset([(frozenset([5, 6]), 1), (frozenset([4, 7]), 1),
                        (frozenset([8, 4]), 1), (frozenset([0, 5]), 1),
                        (frozenset([1, 5]), 1), (frozenset([2, 5]), 1),
                        (frozenset([2, 3]), 1), (frozenset([2, 4]), 2)]))
    assert graph.heavy_atom_indices(mgrph) == (2, 4, 5)


def test__backbone_indices():
    """ test graph.backbone_indices
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([10, 11]), 1), (frozenset([5, 6]), 1),
                        (frozenset([4, 7]), 1), (frozenset([8, 4]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 5]), 1),
                        (frozenset([2, 5]), 1), (frozenset([2, 3]), 1),
                        (frozenset([2, 4]), 2)]))
    assert graph.backbone_indices(mgrph) == (2, 4, 5, 9, 11)


def test__atom_bonds():
    """ test graph.atom_bonds()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.atom_bonds(mgrph1, 0)
            == frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                          (frozenset([0, 1]), 2)]))
    assert graph.atom_bonds(mgrph2, 0) == frozenset()


def test__delete_atom():
    """ test graph.delete_atom()
    """
    mgrph = (('C', 'O', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                        (frozenset([0, 1]), 2)]))
    assert graph.delete_atom(mgrph, 0) == (('O', 'H', 'H'), frozenset([]))
    assert graph.delete_atom(mgrph, 1) == (('C', 'H', 'H'),
                                           frozenset([(frozenset([0, 2]), 1),
                                                      (frozenset([0, 1]), 1)]))
    assert graph.delete_atom(mgrph, 2) == (('C', 'O', 'H'),
                                           frozenset([(frozenset([0, 2]), 1),
                                                      (frozenset([0, 1]), 2)]))


def test__bind_atom():
    """ test graph.bind_atom()
    """
    mgrph = (('O', 'H'), frozenset([(frozenset([0, 1]), 1)]))
    assert (graph.bind_atom(mgrph, 0, 'H', order=1)
            == (('O', 'H', 'H'), frozenset([(frozenset([0, 2]), 1),
                                            (frozenset([0, 1]), 1)])))


def test__atom_neighborhood_indices():
    """ test graph.atom_neighborhood_indices()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert graph.atom_neighborhood_indices(mgrph1, 0) == frozenset([1, 2, 3])
    assert graph.atom_neighborhood_indices(mgrph2, 0) == frozenset()


def test__atom_neighborhood():
    """ test graph.atom_neighborhood()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    assert graph.atom_neighborhood(mgrph1, 0) == ('C', 'H', 'H')


def test__atom_hydrogen_indices():
    """ test graph.atom_hydrogen_indices()
    """
    mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 5]), 1),
                        (frozenset([1, 2]), 1), (frozenset([0, 4]), 1),
                        (frozenset([8, 2]), 1), (frozenset([2, 6]), 1),
                        (frozenset([0, 1]), 2), (frozenset([2, 7]), 1)]))
    assert graph.atom_hydrogen_indices(mgrph, 0) == (3, 4)
    assert graph.atom_hydrogen_indices(mgrph, 1) == (5,)
    assert graph.atom_hydrogen_indices(mgrph, 2) == (8, 6, 7)
    assert graph.atom_hydrogen_indices(mgrph, 3) == ()


def test__atom_hydrogen_count():
    """ test graph.atom_hydrogen_count()
    """
    mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 5]), 1),
                        (frozenset([1, 2]), 1), (frozenset([0, 4]), 1),
                        (frozenset([8, 2]), 1), (frozenset([2, 6]), 1),
                        (frozenset([0, 1]), 2), (frozenset([2, 7]), 1)]))
    assert graph.atom_hydrogen_count(mgrph, 0) == 2
    assert graph.atom_hydrogen_count(mgrph, 1) == 1
    assert graph.atom_hydrogen_count(mgrph, 2) == 3
    assert graph.atom_hydrogen_count(mgrph, 3) == 0


def test__atom_bond_count():
    """ test graph.atom_bond_count()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert graph.atom_bond_count(mgrph1, 0) == 4
    assert graph.atom_bond_count(mgrph2, 0) == 0


def test__atom_free_electrons():
    """ test graph.atom_free_electrons()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    mgrph3 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    assert graph.atom_free_electrons(mgrph1, 0) == 0
    assert graph.atom_free_electrons(mgrph2, 0) == 1
    assert graph.atom_free_electrons(mgrph3, 0) == 2


def test__union():
    """ test graph.union()
    """
    mgrph1 = (('O',), frozenset())
    mgrph2 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    mgrph3 = (('H',), frozenset())
    assert (graph.union(mgrph1, mgrph2)
            == (('O', 'C', 'H', 'H'),
                frozenset([(frozenset([1, 3]), 1), (frozenset([1, 2]), 1)])))
    assert (graph.union(mgrph2, mgrph3)
            == (('C', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)])))


def test__bind():
    """ test graph.bind()
    """
    mgrph1 = (('O',), frozenset())
    mgrph2 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    mgrph3 = (('H',), frozenset())
    assert (graph.bind(mgrph1, mgrph2, 0, 0, order=2)
            == (('O', 'C', 'H', 'H'),
                frozenset([(frozenset([1, 3]), 1), (frozenset([1, 2]), 1),
                           (frozenset([0, 1]), 2)])))
    assert (graph.bind(mgrph2, mgrph3, 0, 0, order=1)
            == (('C', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                           (frozenset([0, 1]), 1)])))


def test__permute_atoms():
    """ test graph.permute_atoms()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert (graph.permute_atoms(mgrph, (5, 2, 1, 4, 3, 0)) ==
            (('H', 'H', 'C', 'H', 'H', 'C'),
             frozenset([(frozenset([2, 5]), 2), (frozenset([2, 3]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([4, 5]), 1)])))


def test__skeleton_sort():
    """ test graph.skeleton_sort
    """
    mgrph = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([10, 11]), 1), (frozenset([5, 6]), 1),
                        (frozenset([4, 7]), 1), (frozenset([8, 4]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 5]), 1),
                        (frozenset([2, 5]), 1), (frozenset([2, 3]), 1),
                        (frozenset([2, 4]), 2)]))
    ref_mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
                 frozenset([(frozenset([8, 2]), 1), (frozenset([11, 4]), 1),
                            (frozenset([1, 6]), 1), (frozenset([1, 7]), 1),
                            (frozenset([0, 2]), 1), (frozenset([0, 5]), 1),
                            (frozenset([0, 1]), 2), (frozenset([2, 10]), 1),
                            (frozenset([9, 2]), 1)]))
    assert graph.skeleton_sort(mgrph) == ref_mgrph


def test__isomorphic():
    """ test graph.isomorphic()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert all(graph.isomorphic(mgrph, graph.permute_atoms(mgrph, p))
               for p in permutations(range(6)))


def test__isomorphism():
    """ test graph.isomorphism()
    """
    mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([2, 4]), 1), (frozenset([0, 6]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 2]), 1),
                        (frozenset([0, 1]), 2), (frozenset([2, 3]), 1),
                        (frozenset([1, 7]), 1)]))

    natms = len(graph.atoms(mgrph))
    for _ in range(10):
        pmt = tuple(numpy.random.permutation(natms))
        pmt_mgrph = graph.permute_atoms(mgrph, pmt)
        iso = graph.isomorphism(mgrph, pmt_mgrph)
        assert graph.permute_atoms(mgrph, iso) == pmt_mgrph


def test__multibond_opening_resonances():
    """ test graph.multibond_opening_resonances
    """
    mgrph1 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
               'H'),
              frozenset([(frozenset([0, 6]), 1), (frozenset([8, 2]), 1),
                         (frozenset([1, 2]), 1), (frozenset([1, 7]), 1),
                         (frozenset([3, 4]), 1), (frozenset([9, 3]), 1),
                         (frozenset([2, 3]), 2), (frozenset([0, 5]), 1),
                         (frozenset([0, 1]), 2), (frozenset([10, 4]), 1),
                         (frozenset([4, 12]), 1), (frozenset([11, 4]), 1)]))
    mgrph2 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
               'H'),
              frozenset([(frozenset([0, 6]), 1), (frozenset([8, 2]), 1),
                         (frozenset([1, 2]), 1), (frozenset([1, 7]), 1),
                         (frozenset([3, 4]), 1), (frozenset([9, 3]), 1),
                         (frozenset([2, 3]), 1), (frozenset([0, 5]), 1),
                         (frozenset([0, 1]), 2), (frozenset([10, 4]), 1),
                         (frozenset([4, 12]), 1), (frozenset([11, 4]), 1)]))
    mgrph3 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
               'H'),
              frozenset([(frozenset([0, 6]), 1), (frozenset([8, 2]), 1),
                         (frozenset([1, 2]), 1), (frozenset([1, 7]), 1),
                         (frozenset([3, 4]), 1), (frozenset([9, 3]), 1),
                         (frozenset([2, 3]), 2), (frozenset([0, 5]), 1),
                         (frozenset([0, 1]), 1), (frozenset([10, 4]), 1),
                         (frozenset([4, 12]), 1), (frozenset([11, 4]), 1)]))

    assert (set(graph.multibond_opening_resonances(mgrph1))
            == {mgrph1, mgrph2, mgrph3})


def test__multibond_forming_resonances():
    """ test graph.multibond_forming_resonances()
    """
    mgrph1 = (('O', 'C', 'O'),
              frozenset([(frozenset([1, 2]), 1), (frozenset([0, 1]), 1)]))
    mgrph2 = (('O', 'C', 'O'),
              frozenset([(frozenset([1, 2]), 1), (frozenset([0, 1]), 2)]))
    mgrph3 = (('O', 'C', 'O'),
              frozenset([(frozenset([1, 2]), 2), (frozenset([0, 1]), 2)]))
    mgrph4 = (('O', 'C', 'O'),
              frozenset([(frozenset([1, 2]), 2), (frozenset([0, 1]), 1)]))
    mgrph5 = (('O', 'C', 'O'),
              frozenset([(frozenset([1, 2]), 2), (frozenset([0, 1]), 2)]))

    assert (set(graph.multibond_forming_resonances(mgrph1))
            == {mgrph1, mgrph2, mgrph3, mgrph4, mgrph5})


def test__forward_abstraction_indices():
    """ test graph.forward_abstraction_indices()
    """
    qh_mgrph = (('C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 6]), 1), (frozenset([1, 2]), 1),
                           (frozenset([0, 5]), 1), (frozenset([8, 1]), 1),
                           (frozenset([1, 7]), 1), (frozenset([2, 3]), 2),
                           (frozenset([10, 3]), 1), (frozenset([0, 4]), 1),
                           (frozenset([3, 11]), 1), (frozenset([9, 2]), 1),
                           (frozenset([0, 1]), 1)]))
    q_mgrph = (('C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([0, 6]), 1), (frozenset([3, 5]), 1),
                          (frozenset([0, 7]), 1), (frozenset([9, 1]), 1),
                          (frozenset([3, 4]), 1), (frozenset([2, 3]), 2),
                          (frozenset([8, 0]), 1), (frozenset([0, 1]), 1),
                          (frozenset([10, 2]), 1), (frozenset([1, 2]), 1)]))
    qh_idx, q_idx = graph.forward_abstraction_indices(qh_mgrph, q_mgrph)
    assert (qh_idx, q_idx) == (8, 1)
    assert graph.isomorphism(q_mgrph, graph.delete_atom(qh_mgrph, qh_idx))
    assert graph.isomorphism(qh_mgrph, graph.bind_atom(q_mgrph, q_idx, 'H'))


def test__addition_indices():
    """ test graph.addition_indices()
    """
    x_mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([0, 3]), 1), (frozenset([1, 5]), 1),
                          (frozenset([1, 2]), 1), (frozenset([0, 4]), 1),
                          (frozenset([8, 2]), 1), (frozenset([2, 6]), 1),
                          (frozenset([0, 1]), 2), (frozenset([2, 7]), 1)]))
    y_mgrph = (('O', 'H'), frozenset([(frozenset([0, 1]), 1)]))
    xy_mgrph = (('C', 'H', 'H', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H'),
                frozenset([(frozenset([10, 5]), 1), (frozenset([3, 5]), 1),
                           (frozenset([3, 6]), 1), (frozenset([4, 7]), 1),
                           (frozenset([0, 2]), 1), (frozenset([3, 4]), 1),
                           (frozenset([0, 3]), 1), (frozenset([8, 5]), 1),
                           (frozenset([9, 5]), 1), (frozenset([0, 1]), 1)]))
    idxs = graph.addition_indices(x_mgrph, y_mgrph, xy_mgrph)
    assert idxs == (1, 0, 3, 4)


def test__migration_indices():
    """ test graph.migration_indices()
    """
    r_mgrph = (('O', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([3, 5]), 1), (frozenset([1, 2]), 1),
                          (frozenset([1, 6]), 1), (frozenset([1, 7]), 1),
                          (frozenset([3, 4]), 1), (frozenset([8, 2]), 1),
                          (frozenset([2, 3]), 1), (frozenset([9, 2]), 1),
                          (frozenset([0, 1]), 1)]))
    p_mgrph = (('O', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([1, 5]), 1), (frozenset([1, 2]), 1),
                          (frozenset([8, 3]), 1), (frozenset([2, 7]), 1),
                          (frozenset([9, 3]), 1), (frozenset([0, 4]), 1),
                          (frozenset([1, 6]), 1), (frozenset([2, 3]), 1),
                          (frozenset([0, 1]), 1)]))
    idxs = graph.migration_indices(r_mgrph, p_mgrph)
    assert idxs == (9, 0, 4, 2)

    r_idx_h, r_idx_a, p_idx_h, p_idx_a = idxs
    r_mgrph_ = graph.delete_atom(graph.bind_atom(p_mgrph, p_idx_a, 'H'),
                                 p_idx_h)
    p_mgrph_ = graph.delete_atom(graph.bind_atom(r_mgrph, r_idx_a, 'H'),
                                 r_idx_h)
    assert graph.isomorphism(r_mgrph, r_mgrph_)
    assert graph.isomorphism(p_mgrph, p_mgrph_)


if __name__ == '__main__':
    test__heavy_atom_indices()
    test__atom_hydrogen_indices()
    test__atom_hydrogen_count()
    test__heavy_atom_subgraph()
    test__backbone_indices()
    test__skeleton_sort()
