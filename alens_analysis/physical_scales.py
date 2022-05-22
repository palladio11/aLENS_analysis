#!/usr/bin/env python

"""@package docstring
File: physical_scales.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import numpy as np


def get_drag_coeff(bead_rad, viscosity):
    ''' Dimensions of Length/(Time*Force) = Length^2/(Time*Energy)'''
    return 6. * np.pi * viscosity * bead_rad


def get_char_time(bead_diam, viscosity, kT=.00411):
    ''' Characteristic time i.e. the time it takes for bead to diffuse its diameter '''
    return 4. * np.pi * viscosity * np.power(bead_diam, 3) / kT


def get_poly_diffuse_dist(time, n_beads, bead_diam,
                          viscosity, frac_dim=.5, kT=.00411):
    return np.sqrt(kT * time / (np.pi * viscosity *
                                bead_diam * np.power(n_beads, frac_dim)))


def get_rouse_time(n_beads, bead_diam, viscosity, kT=.00411):
    return 2. * viscosity * \
        np.power(bead_diam, 3) * n_beads * n_beads / (3 * np.pi * kT)


def get_link_relax_time(bead_rad, viscosity, spring_const):
    return (6. * np.pi * viscosity * bead_rad) / spring_const


def get_pfract(n_chrom, rad_chrom, n_crowd, rad_crowd, rad_sys):
    return (n_chrom * (rad_chrom**3) + n_crowd * (rad_crowd**3)) / (rad_sys**3)


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
