""" test the automechanc.mol.graph module
"""
import numpy
from automechanic.mol import graph


C8H13O_CGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
C8H13O_RGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (2, None), frozenset({3, 5}): (2, None),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
C8H13O_SGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
C8H13O_ICH = ('InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'
              '/b5-3-,6-4-/t8-/m0/s1')

CH2FH2H_CGR_IMP = (
    {1: ('F', 0, None), 3: ('C', 2, None), 4: ('H', 1, None),
     6: ('H', 0, None)},
    {frozenset({1, 3}): (1, None)})
CH2FH2H_CGR_EXP = (
    {0: ('H', 0, None), 1: ('F', 0, None), 2: ('H', 0, None),
     3: ('C', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None)},
    {frozenset({1, 3}): (1, None), frozenset({2, 3}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({4, 5}): (1, None)})
C2H2CL2F2_MM_ICH = 'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2-/m0/s1'
C2H2CL2F2_MM_SGR = (
    {0: ('C', 1, False), 1: ('C', 1, False),
     2: ('F', 0, None), 3: ('Cl', 0, None),
     4: ('F', 0, None), 5: ('Cl', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None)})
C2H2CL2F2_MP_ICH = 'InChI=1S/C2H2Cl2F2/c3-1(5)2(4)6/h1-2H/t1-,2+'
C2H2CL2F2_MP_SGR = (
    {0: ('C', 1, False), 1: ('C', 1, True),
     2: ('F', 0, None), 3: ('Cl', 0, None),
     4: ('F', 0, None), 5: ('Cl', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None)})

C2H2F2_P_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C4H8O_M_ICH = 'InChI=1S/C4H8O/c1-3-4(2)5/h3,5H,1-2H3/b4-3-'
C2H2F2_P_SGR = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('F', 0, None),
                 3: ('F', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
                {frozenset({0, 1}): (1, True), frozenset({0, 2}): (1, None),
                 frozenset({1, 3}): (1, None), frozenset({0, 4}): (1, None),
                 frozenset({1, 5}): (1, None)})
C4H8O_M_SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('C', 0, None), 4: ('O', 1, None), 5: ('H', 0, None)},
               {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
                frozenset({2, 3}): (1, False), frozenset({3, 4}): (1, None),
                frozenset({2, 5}): (1, None)})


# test constructors and value getters
def test__from_data():
    """ test graph.from_data

    also tests a bunch of accessors
    """
    assert graph.from_data(
        graph.atom_symbols(CH2FH2H_CGR_EXP),
        graph.bond_keys(CH2FH2H_CGR_EXP)
    ) == CH2FH2H_CGR_EXP
    assert graph.from_data(
        graph.atom_symbols(C8H13O_CGR),
        graph.bond_keys(C8H13O_CGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR)
    ) == C8H13O_CGR
    assert graph.from_data(
        graph.atom_symbols(C8H13O_RGR),
        graph.bond_keys(C8H13O_RGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGR)
    ) == C8H13O_RGR
    assert graph.from_data(
        graph.atom_symbols(C8H13O_SGR),
        graph.bond_keys(C8H13O_SGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_SGR),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGR),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGR)
    ) == C8H13O_SGR


def test__atom_stereo_keys():
    """ test graph.atom_stereo_keys
    """
    assert graph.atom_stereo_keys(C8H13O_SGR) == (7,)


def test__bond_stereo_keys():
    """ test graph.bond_stereo_keys
    """
    assert (graph.bond_stereo_keys(C8H13O_SGR)
            == (frozenset({2, 4}), frozenset({3, 5})))


# test value setters
def test__set_atom_implicit_hydrogen_valences():
    """ test graph.set_atom_implicit_hydrogen_valences
    """
    assert graph.set_atom_implicit_hydrogen_valences(
        CH2FH2H_CGR_IMP, {3: 1, 4: 0, 6: 1}
    ) == ({1: ('F', 0, None), 3: ('C', 1, None), 4: ('H', 0, None),
           6: ('H', 1, None)},
          {frozenset({1, 3}): (1, None)})


def test__set_atom_stereo_parities():
    """ test graph.set_atom_stereo_parities
    """
    assert graph.atom_stereo_parities(
        graph.set_atom_stereo_parities(C8H13O_CGR, {7: False})
    ) == graph.atom_stereo_parities(C8H13O_SGR)


def test__set_bond_orders():
    """ test graph.set_bond_orders
    """
    assert graph.set_bond_orders(
        C8H13O_CGR, {frozenset({2, 4}): 2, frozenset({3, 5}): 2},
    ) == C8H13O_RGR


def test__set_bond_stereo_parities():
    """ test graph.set_bond_stereo_parities
    """
    assert graph.bond_stereo_parities(
        graph.set_bond_stereo_parities(
            C8H13O_CGR, {frozenset({2, 4}): False, frozenset({3, 5}): False},
        )
    ) == graph.bond_stereo_parities(C8H13O_SGR)


def test__increment_bond_orders():
    """ test graph.increment_bond_orders
    """
    assert graph.increment_bond_orders(
        C8H13O_CGR, {frozenset({2, 4}): 1, frozenset({3, 5}): 1}
    ) == C8H13O_RGR


# test derived values
def test__is_chiral():
    """ test graph.is_chiral
    """
    assert graph.is_chiral(C8H13O_SGR) is True
    assert graph.is_chiral(C2H2CL2F2_MM_SGR) is True
    assert graph.is_chiral(C2H2CL2F2_MP_SGR) is False


def test__maximum_spin_multiplicity():
    """ test graph.maximum_spin_multiplicity
    """
    catm_cgr = ({0: ('C', 0, None)}, {})
    ch0f1_cgr = ({0: ('C', 0, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch1f1_cgr = ({0: ('C', 1, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch2f1_cgr = ({0: ('C', 2, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch2f2_cgr = ({0: ('C', 2, None), 1: ('F', 0, None), 2: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None), frozenset([0, 2]): (1, None)})
    o2_cgr = ({0: ('O', 0, None), 1: ('O', 0, None)},
              {frozenset([0, 1]): (1, None)})
    assert graph.maximum_spin_multiplicity(catm_cgr) == 5
    assert graph.maximum_spin_multiplicity(ch0f1_cgr) == 4
    assert graph.maximum_spin_multiplicity(ch1f1_cgr) == 3
    assert graph.maximum_spin_multiplicity(ch2f1_cgr) == 2
    assert graph.maximum_spin_multiplicity(ch2f2_cgr) == 1
    assert graph.maximum_spin_multiplicity(o2_cgr) == 3


def test__possible_spin_multiplicities():
    """ test graph.possible_spin_multiplicities
    """
    catm_cgr = ({0: ('C', 0, None)}, {})
    ch0f1_cgr = ({0: ('C', 0, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch1f1_cgr = ({0: ('C', 1, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch2f1_cgr = ({0: ('C', 2, None), 1: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None)})
    ch2f2_cgr = ({0: ('C', 2, None), 1: ('F', 0, None), 2: ('F', 0, None)},
                 {frozenset([0, 1]): (1, None), frozenset([0, 2]): (1, None)})
    o2_cgr = ({0: ('O', 0, None), 1: ('O', 0, None)},
              {frozenset([0, 1]): (1, None)})
    assert graph.possible_spin_multiplicities(catm_cgr) == (1, 3, 5)
    assert graph.possible_spin_multiplicities(ch0f1_cgr) == (2, 4)
    assert graph.possible_spin_multiplicities(ch1f1_cgr) == (1, 3)
    assert graph.possible_spin_multiplicities(ch2f1_cgr) == (2,)
    assert graph.possible_spin_multiplicities(ch2f2_cgr) == (1,)
    assert graph.possible_spin_multiplicities(o2_cgr) == (1, 3)


def test__ring_keys_list():
    """ test graph.ring_keys_list
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
            6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
            9: ('N', 1, None), 10: ('O', 1, None)},
           {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
            frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
            frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    assert graph.ring_keys_list(cgr) == ((0, 1, 3, 6, 7), (1, 2, 3, 4, 8, 9))


def test__backbone_keys():
    """ test graph.backbone_keys
    """
    assert graph.backbone_keys(CH2FH2H_CGR_EXP) == (1, 3, 4, 6)


def test__explicit_hydrogen_keys():
    """ test graph.explicit_hydrogen_keys
    """
    assert graph.explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == (0, 2, 5)


def test__atom_nuclear_charges():
    """ test graph.atom_nuclear_charges
    """
    assert (graph.atom_nuclear_charges(C8H13O_CGR) ==
            {0: 6, 1: 6, 2: 6, 3: 6, 4: 6, 5: 6, 6: 6, 7: 6, 8: 8})


def test__atom_total_valences():
    """ test graph.atom_total_valences
    """
    assert (graph.atom_total_valences(C8H13O_CGR) ==
            {0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 2})


def test__atom_bond_valences():
    """ test graph.atom_bond_valences
    """
    assert (graph.atom_bond_valences(C8H13O_CGR) ==
            {0: 4, 1: 4, 2: 3, 3: 3, 4: 3, 5: 3, 6: 4, 7: 4, 8: 1})


def test__atom_radical_valences():
    """ test graph.atom_radical_valences
    """
    assert (graph.atom_radical_valences(C8H13O_CGR) ==
            {0: 0, 1: 0, 2: 1, 3: 1, 4: 1, 5: 1, 6: 0, 7: 0, 8: 1})


def test__atom_neighbor_keys():
    """ test graph.atom_neighbor_keys
    """
    assert (graph.atom_neighbor_keys(C8H13O_CGR) ==
            {0: (2,), 1: (3,), 2: (0, 4), 3: (1, 5), 4: (2, 6), 5: (3, 7),
             6: (4, 7), 7: (5, 6, 8), 8: (7,)})


def test__atom_explicit_hydrogen_keys():
    """ test graph.atom_explicit_hydrogen_keys
    """
    assert (graph.atom_explicit_hydrogen_keys(CH2FH2H_CGR_EXP) ==
            {0: (), 1: (), 2: (), 3: (0, 2), 4: (5,), 5: (), 6: ()})


def test__atom_bond_keys():
    """ test graph.atom_bond_keys
    """
    assert (graph.atom_bond_keys(C8H13O_CGR) ==
            {0: (frozenset({0, 2}),),
             1: (frozenset({1, 3}),),
             2: (frozenset({0, 2}), frozenset({2, 4})),
             3: (frozenset({1, 3}), frozenset({3, 5})),
             4: (frozenset({2, 4}), frozenset({4, 6})),
             5: (frozenset({3, 5}), frozenset({5, 7})),
             6: (frozenset({4, 6}), frozenset({6, 7})),
             7: (frozenset({5, 7}), frozenset({6, 7}), frozenset({8, 7})),
             8: (frozenset({8, 7}),)})


def test__atom_neighborhoods():
    """ test graph.atom_neighborhoods
    """
    assert (graph.atom_neighborhoods(C8H13O_CGR) == {
        0: ({0: ('C', 3, None), 2: ('C', 1, None)},
            {frozenset({0, 2}): (1, None)}),
        1: ({1: ('C', 3, None), 3: ('C', 1, None)},
            {frozenset({1, 3}): (1, None)}),
        2: ({0: ('C', 3, None), 2: ('C', 1, None), 4: ('C', 1, None)},
            {frozenset({0, 2}): (1, None), frozenset({2, 4}): (1, None)}),
        3: ({1: ('C', 3, None), 3: ('C', 1, None), 5: ('C', 1, None)},
            {frozenset({1, 3}): (1, None), frozenset({3, 5}): (1, None)}),
        4: ({2: ('C', 1, None), 4: ('C', 1, None), 6: ('C', 2, None)},
            {frozenset({2, 4}): (1, None), frozenset({4, 6}): (1, None)}),
        5: ({3: ('C', 1, None), 5: ('C', 1, None), 7: ('C', 1, None)},
            {frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)}),
        6: ({4: ('C', 1, None), 6: ('C', 2, None), 7: ('C', 1, None)},
            {frozenset({4, 6}): (1, None), frozenset({6, 7}): (1, None)}),
        7: ({8: ('O', 0, None), 5: ('C', 1, None), 6: ('C', 2, None),
             7: ('C', 1, None)},
            {frozenset({5, 7}): (1, None), frozenset({6, 7}): (1, None),
             frozenset({8, 7}): (1, None)}),
        8: ({8: ('O', 0, None), 7: ('C', 1, None)},
            {frozenset({8, 7}): (1, None)})})


def test__atom_inchi_numbers():
    """ test graph.atom_inchi_numbers
    """
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        inv_pmt_dct = dict(map(reversed, pmt_dct.items()))
        assert graph.atom_inchi_numbers(cgr_pmt) == inv_pmt_dct

    ch_cgr = ({5: ('C', 1, None)}, {})
    assert graph.atom_inchi_numbers(ch_cgr) == {5: 0}

    ch_cgr = ({5: ('C', 0, None), 2: ('H', 0, None)},
              {frozenset({5, 2}): (1, None)})
    assert graph.atom_inchi_numbers(ch_cgr) == {5: 0, 2: -1}

    cf_cgr = ({5: ('C', 0, None), 2: ('F', 0, None)},
              {frozenset({5, 2}): (1, None)})
    assert graph.atom_inchi_numbers(cf_cgr) == {5: 0, 2: 1}

    ccl_cgr = ({5: ('C', 0, None), 2: ('F', 0, None)},
               {frozenset({5, 2}): (1, None)})
    assert graph.atom_inchi_numbers(ccl_cgr) == {5: 0, 2: 1}


def test__inchi():
    """ test graph.inchi
    """
    co_cgr = ({0: ('C', 0, None), 1: ('O', 0, None)},
              {frozenset({0, 1}): (1, None)})
    assert graph.inchi(co_cgr) == 'InChI=1S/CO/c1-2'

    assert graph.inchi(C8H13O_SGR) == (
        'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')

    c_cgr = ({5: ('C', 0, None)}, {})
    assert graph.inchi(c_cgr) == 'InChI=1S/C'

    n_cgr = ({5: ('N', 0, None)}, {})
    assert graph.inchi(n_cgr) == 'InChI=1S/N'

    ch_cgr = ({5: ('C', 1, None)}, {})
    assert graph.inchi(ch_cgr) == 'InChI=1S/CH/h1H'

    ch_cgr = ({5: ('C', 0, None), 2: ('H', 0, None)},
              {frozenset({5, 2}): (1, None)})
    assert graph.inchi(ch_cgr) == 'InChI=1S/CH/h1H'

    cf_cgr = ({5: ('C', 0, None), 2: ('F', 0, None)},
              {frozenset({5, 2}): (1, None)})
    assert graph.inchi(cf_cgr) == 'InChI=1S/CF/c1-2'

    ccl_cgr = ({5: ('C', 0, None), 2: ('Cl', 0, None)},
               {frozenset({5, 2}): (1, None)})
    assert graph.inchi(ccl_cgr) == 'InChI=1S/CCl/c1-2'

    nh_cgr = ({5: ('N', 1, None)}, {})
    assert graph.inchi(nh_cgr) == 'InChI=1S/HN/h1H'

    ch2_cgr = ({5: ('C', 2, None)}, {})
    assert graph.inchi(ch2_cgr) == 'InChI=1S/CH2/h1H2'


def test__stereo_inchi():
    """ test graph.stereo_inchi
    """
    assert graph.stereo_inchi(C8H13O_SGR) == C8H13O_ICH
    assert graph.stereo_inchi(C2H2CL2F2_MM_SGR) == C2H2CL2F2_MM_ICH
    assert graph.stereo_inchi(C2H2CL2F2_MP_SGR) == C2H2CL2F2_MP_ICH
    assert graph.stereo_inchi(C2H2F2_P_SGR) == C2H2F2_P_ICH
    assert graph.stereo_inchi(C4H8O_M_SGR) == C4H8O_M_ICH


# test transformations
def test__implicit():
    """ test graph.implicit
    """
    assert graph.implicit(CH2FH2H_CGR_EXP) == CH2FH2H_CGR_IMP
    assert graph.implicit(CH2FH2H_CGR_EXP, (1, 3, 4, 6)) == CH2FH2H_CGR_IMP


def test__explicit():
    """ test graph.explicit
    """
    ch2fh2h_cgr_exp = graph.explicit(CH2FH2H_CGR_IMP)
    assert graph.backbone_isomorphic(ch2fh2h_cgr_exp, CH2FH2H_CGR_EXP)
    assert (graph.atom_explicit_hydrogen_keys(ch2fh2h_cgr_exp) ==
            {1: (), 3: (7, 8), 4: (9,), 6: (), 7: (), 8: (), 9: ()})


def test__explicit_stereo_sites():
    """ test graph.explicit_stereo_sites
    """
    assert graph.explicit_stereo_sites(C8H13O_CGR) == C8H13O_CGR
    assert (graph.explicit_stereo_sites(C8H13O_SGR)
            == ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                 3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
                 6: ('C', 2, None), 7: ('C', 0, False), 8: ('O', 0, None),
                 9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
                 12: ('H', 0, None), 13: ('H', 0, None)},
                {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
                 frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
                 frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
                 frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
                 frozenset({9, 7}): (1, None), frozenset({2, 10}): (1, None),
                 frozenset({3, 11}): (1, None), frozenset({4, 12}): (1, None),
                 frozenset({5, 13}): (1, None)}))


def test__delete_atoms():
    """ test graph.delete_atoms
    """
    assert (graph.delete_atoms(CH2FH2H_CGR_EXP, (0, 2, 5)) ==
            ({1: ('F', 0, None), 3: ('C', 0, None), 4: ('H', 0, None),
              6: ('H', 0, None)},
             {frozenset({1, 3}): (1, None)}))


def test__add_explicit_hydrogens():
    """ test graph.add_explicit_hydrogens
    """
    assert graph.add_explicit_hydrogens(
        CH2FH2H_CGR_IMP, {3: 2, 4: 1}
    ) == ({1: ('F', 0, None), 3: ('C', 2, None), 4: ('H', 1, None),
           6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
           9: ('H', 0, None)},
          {frozenset({1, 3}): (1, None), frozenset({3, 7}): (1, None),
           frozenset({8, 3}): (1, None), frozenset({9, 4}): (1, None)})


def test__subgraph():
    """ test graph.subgraph
    """
    assert (graph.subgraph(CH2FH2H_CGR_EXP, (1, 3, 4, 6)) ==
            ({1: ('F', 0, None), 3: ('C', 0, None), 4: ('H', 0, None),
              6: ('H', 0, None)},
             {frozenset({1, 3}): (1, None)}))


def test__subgraph_by_bonds():
    """ test graph.subgraph_by_bonds
    """
    assert (graph.subgraph_by_bonds(C8H13O_CGR,
                                    {frozenset({1, 3}), frozenset({3, 5}),
                                     frozenset({5, 7}), frozenset({8, 7})}) ==
            ({1: ('C', 3, None), 3: ('C', 1, None), 5: ('C', 1, None),
              7: ('C', 1, None), 8: ('O', 0, None)},
             {frozenset({1, 3}): (1, None), frozenset({3, 5}): (1, None),
              frozenset({5, 7}): (1, None), frozenset({8, 7}): (1, None)}))


def test__relabel():
    """ test graph.relabel
    """
    assert graph.relabel(
        CH2FH2H_CGR_IMP, {1: 0, 3: 1, 4: 2, 6: 3}
    ) == ({0: ('F', 0, None), 1: ('C', 2, None), 2: ('H', 1, None),
           3: ('H', 0, None)},
          {frozenset({0, 1}): (1, None)})


def test__subresonances():
    """ test graph.subresonances
    """
    c2_cgr = ({0: ('C', 0, None), 1: ('C', 0, None)},
              {frozenset({0, 1}): (1, None)})
    assert graph.subresonances(c2_cgr) == (
        ({0: ('C', 0, None), 1: ('C', 0, None)},
         {frozenset({0, 1}): (1, None)}),
        ({0: ('C', 0, None), 1: ('C', 0, None)},
         {frozenset({0, 1}): (2, None)}),
        ({0: ('C', 0, None), 1: ('C', 0, None)},
         {frozenset({0, 1}): (3, None)}),
        ({0: ('C', 0, None), 1: ('C', 0, None)},
         {frozenset({0, 1}): (4, None)}),
    )

    c3h3_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
                {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
                 frozenset({2, 0}): (1, None)})

    assert graph.subresonances(c3h3_cgr) == (
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
         {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
          frozenset({0, 2}): (1, None)}),
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
         {frozenset({0, 1}): (1, None), frozenset({1, 2}): (2, None),
          frozenset({0, 2}): (1, None)}),
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
         {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
          frozenset({0, 2}): (2, None)}),
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
         {frozenset({0, 1}): (2, None), frozenset({1, 2}): (1, None),
          frozenset({0, 2}): (1, None)}),
    )


def test__lowspin_resonance():
    """ test graph.lowspin_resonance
    """
    c6h6_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                 3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
                {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
                 frozenset({2, 3}): (1, None), frozenset({3, 4}): (1, None),
                 frozenset({4, 5}): (1, None), frozenset({5, 0}): (1, None)})
    assert graph.lowspin_resonance(c6h6_cgr) in [
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
          3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
         {frozenset({0, 1}): (2, None), frozenset({1, 2}): (1, None),
          frozenset({2, 3}): (2, None), frozenset({3, 4}): (1, None),
          frozenset({4, 5}): (2, None), frozenset({5, 0}): (1, None)}),
        ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
          3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
         {frozenset({0, 1}): (1, None), frozenset({1, 2}): (2, None),
          frozenset({2, 3}): (1, None), frozenset({3, 4}): (2, None),
          frozenset({4, 5}): (1, None), frozenset({5, 0}): (2, None)})
    ]


def test__reflection():
    """ test graph.reflection
    """
    assert (graph.reflection(C8H13O_SGR) ==
            graph.set_atom_stereo_parities(C8H13O_SGR, {7: True}))


# test comparisons
def test__backbone_isomorphic():
    """ test graph.backbone_isomorphic
    """
    assert graph.backbone_isomorphic(CH2FH2H_CGR_EXP, CH2FH2H_CGR_IMP)
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphic(cgr, cgr_pmt)


def test__backbone_isomorphism():
    """ test graph.backbone_isomorphism
    """
    assert (graph.backbone_isomorphism(CH2FH2H_CGR_EXP, CH2FH2H_CGR_IMP) ==
            {1: 1, 3: 3, 4: 4, 6: 6})
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphism(cgr, cgr_pmt) == pmt_dct


if __name__ == '__main__':
    # test constructors and value getters
    test__from_data()
    test__atom_stereo_keys()
    test__bond_stereo_keys()
    # test value setters
    test__set_atom_implicit_hydrogen_valences()
    test__set_atom_stereo_parities()
    test__set_bond_orders()
    test__set_bond_stereo_parities()
    test__increment_bond_orders()
    # test derived values
    test__is_chiral()
    test__maximum_spin_multiplicity()
    test__possible_spin_multiplicities()
    test__ring_keys_list()
    test__backbone_keys()
    test__explicit_hydrogen_keys()
    test__atom_nuclear_charges()
    test__atom_total_valences()
    test__atom_bond_valences()
    test__atom_radical_valences()
    test__atom_neighbor_keys()
    test__atom_explicit_hydrogen_keys()
    test__atom_bond_keys()
    test__atom_neighborhoods()
    test__atom_inchi_numbers()
    test__inchi()
    test__stereo_inchi()
    # test transformations
    test__implicit()
    test__explicit()
    test__explicit_stereo_sites()
    test__delete_atoms()
    test__add_explicit_hydrogens()
    test__subgraph()
    test__subgraph_by_bonds()
    test__relabel()
    test__subresonances()
    test__lowspin_resonance()
    test__reflection()
    # test comparisons
    test__backbone_isomorphic()
    test__backbone_isomorphism()
