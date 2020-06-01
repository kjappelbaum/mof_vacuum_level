# -*- coding: utf-8 -*-
"""Calculate the vacuum level using the lowest variance estimate from a cube file of a porous MOF"""
from __future__ import absolute_import

import math

import macrodensity as md
import numpy as np
from macrodensity.cp2k_tools import cell_to_cellpar, read_cube_density, test_point

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


class MOFVacLevel:
    """Vacuum level from the cube file of a porous MOF"""

    def __init__(self, input_cube: str, cube_size: list = [35, 35, 35]):
        self.input_cube = input_cube
        self.cube_size = cube_size
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

        self.d_l = None
        self.f_l = None
        self.g_l = None
        self.cube_potential = None
        self.cube_variance = None

    def _find_lowest_var_point(self, threshold: float = 2.0):
        travelled = [0, 0, 0]
        cube_var_l = 1000
        cube_potential_l = 0
        d_l = 0
        f_l = 0
        g_l = 0
        logical = 0

        c1 = 0.2 / self.vector_a
        c2 = 0.2 / self.vector_b
        c3 = 0.2 / self.vector_c

        for d in np.arange(0.0, 1.0, c1):
            for f in np.arange(0.0, 1.0, c2):
                for g in np.arange(0.0, 1.0, c3):
                    cube_origin = [d, f, g]
                    logical = test_point(
                        cube_origin, self.coord, self.params, self.num_atoms, threshold
                    )
                    if logical == 1:
                        cube_potential, cube_var = md.cube_potential(
                            cube_origin,
                            travelled,
                            self.cube_size,
                            self.grid_pot,
                            self.ngx,
                            self.ngy,
                            self.ngz,
                        )
                        if cube_var < cube_var_l:
                            cube_var_l = cube_var
                            cube_potential_l = cube_potential
                            d_l = d
                            f_l = f
                            g_l = g

        self.d_l = d_l
        self.f_l = f_l
        self.g_l = g_l
        self.cube_potential = cube_potential_l
        self.cube_variance = cube_var_l

    def get_vacuum_potential(self):
        self._find_lowest_var_point()
        return self.cube_potential

    @property
    def vacuum_potential(self):
        return self.cube_potential
