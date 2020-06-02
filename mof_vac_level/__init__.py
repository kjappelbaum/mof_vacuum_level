# -*- coding: utf-8 -*-
"""Calculate the vacuum level using the lowest variance estimate from a cube file of a porous MOF"""
from __future__ import absolute_import

import math

from . import macrodensity as md
import numpy as np
from .macrodensity.cp2k_tools import cell_to_cellpar, read_cube_density, test_point
import multiprocess as mp
import itertools
from functools import partial
from tqdm import tqdm

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def _nested_search(
    res,
    vector_a,
    vector_b,
    vector_c,
    coord,
    params,
    num_atoms,
    threshold,
    cube_size,
    grid_pot,
    ngx,
    ngy,
    ngz,
):
    pool = mp.Pool()

    c1 = res / vector_a
    c2 = res / vector_b
    c3 = res / vector_c

    partial_point_test = partial(
        _test_one_point,
        coord=coord,
        params=params,
        num_atoms=num_atoms,
        threshold=threshold,
        cube_size=cube_size,
        grid_pot=grid_pot,
        ngx=ngx,
        ngy=ngy,
        ngz=ngz,
    )

    dimension_1 = np.arange(0.0, 1.0, c1)
    dimension_2 = np.arange(0.0, 1.0, c2)
    dimension_3 = np.arange(0.0, 1.0, c3)

    var = np.zeros([len(dimension_1), len(dimension_2), len(dimension_3)])
    potential = np.zeros([len(dimension_1), len(dimension_2), len(dimension_3)])

    for ii, i in tqdm(enumerate(dimension_1)):
        for jj, j in enumerate(dimension_2):
            result = pool.map(partial_point_test, [(i, j, k) for k in dimension_3])

            potentials, variances = [], []
            for k, resultpoint in enumerate(result):

                var[ii][jj][k] = resultpoint[1]
                potential[ii][jj][k] = resultpoint[0]

    minimum_variance = np.argmin(var)
    return var[minimum_variance], potential[minimum_variance], minimum_variance


def _test_one_point(
    cube_origin, coord, params, num_atoms, threshold, cube_size, grid_pot, ngx, ngy, ngz
):
    travelled = [0, 0, 0]
    logical = test_point(cube_origin, coord, params, num_atoms, threshold)
    if logical == 1:

        cube_potential, cube_var = md.cube_potential(
            cube_origin, travelled, cube_size, grid_pot, ngx, ngy, ngz,
        )
        return cube_potential, cube_var

    return np.nan, np.nan


class MOFVacLevel:
    """Vacuum level from the cube file of a porous MOF"""

    def __init__(self, input_cube: str):
        self.input_cube = input_cube

        (
            self.pot,
            self.ngx,
            self.ngy,
            self.ngz,
            self.lattice,
            self.num_atoms,
            self.coord,
            _,
        ) = read_cube_density(self.input_cube)

        (
            self.vector_a,
            self.vector_b,
            self.vector_c,
            self.av,
            self.bv,
            self.cv,
        ) = md.matrix_2_abc(self.lattice)

        self.grid_pot = md.density_grid_cube(self.pot, self.ngx, self.ngy, self.ngz)
        self.params = cell_to_cellpar(self.lattice, radians=False)

        self.minimum_variance_coords = None
        self.cube_potential = None
        self.cube_variance = None

    def _find_lowest_var_point(
        self, threshold: float = 4.0, res: float = 0.2, cube_size: list = [30, 30, 30]
    ):
        lowest_variance, potential, minimum_variance_coords = _nested_search(
            res,
            self.vector_a,
            self.vector_b,
            self.vector_c,
            self.coord,
            self.params,
            self.num_atoms,
            threshold,
            cube_size,
            self.grid_pot,
            self.ngx,
            self.ngy,
            self.ngz,
        )

        self.lowest_variance = lowest_variance
        self.cube_potential = potential
        self.minimum_variance_coords = minimum_variance_coords

    def get_vacuum_potential(
        self, threshold: float = 2.0, res: float = 0.2, cube_size: list = [30, 30, 30]
    ):
        self._find_lowest_var_point(threshold, res, cube_size)
        return self.cube_potential

    @property
    def vacuum_potential(self):
        return self.cube_potential

    @property
    def minimum_variance(self):
        return self.minimum_variance

    @property
    def minimum_variance_indices(self):
        return self.minimum_variance_coords
