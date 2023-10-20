#!/usr/bin/env python

"""@package docstring
File: physical_scales.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import numpy as np
from scipy.special import erf
import yaml
import h5py


def get_drag_coeff(bead_rad, viscosity):
    ''' Drag coefficient of a sphere of radius bead_rad in a viscous fluid using Stokes-Einstein relations. Dimensions of Length/(Time*Force) = Length^2/(Time*Energy)'''
    return 6. * np.pi * viscosity * bead_rad


def get_char_time(bead_diam, viscosity, kT=.00411):
    ''' Characteristic time i.e. the time it takes for bead to diffuse its diameter '''
    return 4. * np.pi * viscosity * np.power(bead_diam, 3) / kT


def get_poly_diffuse_dist(time, n_beads, bead_diam,
                          viscosity, frac_dim=.5, kT=.00411):
    return np.sqrt(kT * time / (np.pi * viscosity *
                                bead_diam * np.power(n_beads, frac_dim)))


def get_rouse_time(n_beads, bead_diam, viscosity, kT=.00411):
    return viscosity * np.power(bead_diam, 3) * n_beads * n_beads / (2 * np.pi * kT)


def get_link_relax_time(bead_rad, viscosity, spring_const):
    return (6. * np.pi * viscosity * bead_rad) / spring_const


def get_pfract(n_chrom, rad_chrom, n_crowd, rad_crowd, rad_sys):
    """Packing fraction calculation. 

    Parameters
    ----------
    n_chrom : int
        Number of filaments
    rad_chrom : float
        Radius of the beads in the filament
    n_crowd : int
        Number of crowding particles
    rad_crowd : float
        Radius of the crowding particles
    rad_sys : float
        Radius of the confining sphere


    Returns
    -------
    float
        Packing fraction of the beads to total space available.
    """
    return (n_chrom * (rad_chrom**3) + n_crowd * (rad_crowd**3)) / (rad_sys**3)

def calc_sticky_search_volume(l, ks, kbT=4.11e-3):
    a = .5 * ks / kbT
    return 4*np.pi*((1. + (2*a*l*l))*np.sqrt(np.pi)*(1. + erf(np.sqrt(a)*l))) /(4.*np.power(a, 1.5)) + (3. * np.exp(-a*l*l)*l)/(2*a)
    



def get_fundamental_consts(h5_path):
    with h5py.File(h5_path, 'r+') as h5_data:
        rc_dict = yaml.safe_load(h5_data.attrs['RunConfig'])
        # p_dict = yaml.safe_load(h5_data.attrs['ProteinConfig'])
        ks = float(rc_dict['linkKappa'])
        visc = float(rc_dict['viscosity'])
        kbT = float(rc_dict['KBT'])
        diam = float(rc_dict['sylinderDiameter'])
        nbeads = h5_data['raw_data']['sylinders'].shape[0]
        const_dict = {'drag_coeff': get_drag_coeff(diam * .5, visc),
                      'bead_diff_time': get_char_time(diam, visc, kbT),
                      'rouse_time': get_rouse_time(nbeads, diam, visc, kbT),
                      'link_relax_time': get_link_relax_time(diam * .5, visc, ks)}
        return const_dict


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
