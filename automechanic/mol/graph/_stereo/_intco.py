""" coordinate library
"""
from itertools import chain as _chain
from functools import partial as _partial
from ._intco_linalg import unit_direction
from ._intco_linalg import aligning_rotation_matrix
from ._intco_linalg import local_coordinate_interpreter
from .._conn import atom_inchi_numbers
from .._base import atom_keys
from .._base import atom_stereo_keys
from .._base import bond_stereo_keys
from .._base import atom_stereo_parities
from .._base import bond_stereo_parities
from .._base import atom_neighbor_keys
from .._base import explicit_stereo_sites
from .._base import ring_keys_list


def atom_stereo_coordinates(sgr):
    """ determine stereo-specific coordinates for this molecular graph
    """

    if ring_keys_list(sgr):
        raise NotImplementedError("Not yet implemented for rings")

    assert sgr == explicit_stereo_sites(sgr)

    atm_ich_num_dct = atom_inchi_numbers(sgr)
    key_sorter = _partial(sorted, key=atm_ich_num_dct.__getitem__)

    atm_ngb_keys_dct = atom_neighbor_keys(sgr)

    atm_ste_keys = atom_stereo_keys(sgr)
    bnd_ste_keys = bond_stereo_keys(sgr)
    atm_par_dct = atom_stereo_parities(sgr)
    bnd_par_dct = bond_stereo_parities(sgr)

    def _recurse_coordinates(atm_key, anchor_key, xyz_dct):
        boundary_edges = ()
        bnd_key = next(
            filter(lambda bnd_key: atm_key in bnd_key, bnd_ste_keys), None)

        assert bnd_key is None or atm_key not in atm_ste_keys

        if atm_key in atm_ste_keys:
            atm_par = atm_par_dct[atm_key]
            xyz_dct, boundary_edges = _atom_stereo_coordinates(
                anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct, key_sorter,
                atm_par)
        elif bnd_key is not None:
            bnd_par = bnd_par_dct[bnd_key]
            xyz_dct, boundary_edges = _bond_stereo_coordinates(
                anchor_key, bnd_key, atm_ngb_keys_dct, xyz_dct, key_sorter,
                bnd_par)
        else:
            xyz_dct, boundary_edges = _nonstereo_coordinates(
                anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct)

        for anchor_key_, atm_key_ in boundary_edges:
            xyz_dct = _recurse_coordinates(atm_key_, anchor_key_, xyz_dct)

        return xyz_dct

    xyz_dct = {}

    atm_keys = key_sorter(atom_keys(sgr))

    if atm_keys:
        atm_key = next(iter(atm_keys))
        xyz_dct[atm_key] = (0, 0, 0)

        atm_ngb_keys = atm_ngb_keys_dct[atm_key]

        if atm_ngb_keys:
            atm_ngb_key = next(iter(atm_ngb_keys))
            xyz_dct[atm_ngb_key] = (1, 0, 0)

            xyz_dct = _recurse_coordinates(atm_ngb_key, atm_key, xyz_dct)
            xyz_dct = _recurse_coordinates(atm_key, atm_ngb_key, xyz_dct)

    return xyz_dct


def _nonstereo_coordinates(anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct):
    """ assign non-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),    # atm 1
                    (-1, 0, 0),   # atm 1 ngb 0
                    (1, 0, 0),    # atm 1 ngb 1
                    (0, 1, 0),    # atm 1 ngb 2
                    (0, -1, 0))   # atm 1 ngb 3

    assert atm_key in atm_ngb_keys_dct
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    stencil_keys = list(_chain([atm_key], atm_ngb_keys))

    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple((atm_key, ngb_key) for ngb_key in atm_ngb_keys
                           if ngb_key != anchor_key)

    return xyz_dct, boundary_edges


def _atom_stereo_coordinates(anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct,
                             key_sorter, parity):
    """ assign atom-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),                     # atm 1
                    (-1, 0, 0),                    # atm 1 ngb 0
                    (0, 1, 0),                     # atm 1 ngb 1
                    (0, 0, (-1)**(not parity)),    # atm 1 ngb 2
                    (0, -1, 0))                    # atm 1 ngb 3

    assert atm_key in atm_ngb_keys_dct
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    stencil_keys = list(_chain([atm_key], key_sorter(atm_ngb_keys)))

    assert len(stencil_keys) == len(stencil_xyzs)
    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple((atm_key, ngb_key) for ngb_key in atm_ngb_keys
                           if ngb_key != anchor_key)

    return xyz_dct, boundary_edges


def _bond_stereo_coordinates(anchor_key, bnd_key, atm_ngb_keys_dct, xyz_dct,
                             key_sorter, parity):
    """ assign bond-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),                     # atm 1
                    (0, 1, 0),                     # atm 2
                    (-1, 0, 0),                    # atm 1 ngb 0 (anchor?)
                    (1, 0, 0),                     # atm 1 ngb 1 (anchor?)
                    ((-1) ** (not parity), 1, 0),  # atm 2 ngb 0
                    ((-1) ** parity, 1, 0))        # atm 2 ngb 1

    atm1_key, atm2_key = bnd_key
    assert atm1_key in atm_ngb_keys_dct
    assert atm2_key in atm_ngb_keys_dct
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key]
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]

    if anchor_key in atm2_ngb_keys:
        atm1_key, atm2_key = atm2_key, atm1_key
        atm1_ngb_keys, atm2_ngb_keys = atm2_ngb_keys, atm1_ngb_keys

    stencil_keys = list(_chain(
        [atm1_key, atm2_key],
        key_sorter(filter(lambda x: x != atm2_key, atm1_ngb_keys)),
        key_sorter(filter(lambda x: x != atm1_key, atm2_ngb_keys))))

    assert len(stencil_keys) == len(stencil_xyzs)
    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm1_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple(_chain(
        ((atm1_key, ngb_key) for ngb_key in atm1_ngb_keys
         if ngb_key not in (anchor_key, atm2_key)),
        ((atm2_key, ngb_key) for ngb_key in atm2_ngb_keys
         if ngb_key != atm1_key)))

    return xyz_dct, boundary_edges


def _from_stencil(atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs):
    assert atm_key in xyz_dct
    assert anchor_key in xyz_dct

    true_atm_xyz = xyz_dct[atm_key]
    true_anchor_xyz = xyz_dct[anchor_key]
    true_bond_xyz = unit_direction(true_atm_xyz, true_anchor_xyz)

    assert atm_key in stencil_keys
    assert anchor_key in stencil_keys
    stencil_atm_xyz = stencil_xyzs[stencil_keys.index(atm_key)]
    stencil_anchor_xyz = stencil_xyzs[stencil_keys.index(anchor_key)]
    stencil_bond_xyz = unit_direction(stencil_atm_xyz, stencil_anchor_xyz)

    interp_stencil_ = local_coordinate_interpreter(
        trans=true_atm_xyz,
        rot=aligning_rotation_matrix(stencil_bond_xyz, true_bond_xyz)
    )

    assert len(stencil_keys) <= len(stencil_xyzs)
    ret_xyz_dct = dict(zip(stencil_keys, map(interp_stencil_, stencil_xyzs)))

    # check that the transformation worked
    assert ret_xyz_dct[atm_key] == true_atm_xyz
    assert ret_xyz_dct[anchor_key] == true_anchor_xyz

    return ret_xyz_dct
