#!/usr/bin/env python

"""@package docstring
File: chrom_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
import re
import time
import yaml
from pprint import pprint
from pathlib import Path
import h5py

# Data manipulation
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
import scipy.stats as stats


def get_rouse_modes_at_t(pos_arr, nmodes=20):
    """TODO: Docstring for get_rouse_modes.

    @param sphere_dat TODO
    @param nmodes TODO
    @return: TODO

    """
    modes = []
    mode_0 = pos_arr[0]
    nbeads = pos_arr.shape[0]

    for k in range(nmodes):
        modes += [np.zeros(3)]
        for n in range(nbeads - 1):
            modes[-1] += ((pos_arr[n] - mode_0) *
                          np.cos(np.pi * (n + .5) * k / nbeads))

    return np.asarray(modes) / (nbeads)


def get_rouse_modes(pos_mat, nmodes=20):
    """TODO: Docstring for get_rouse_modes.

    @param pos_mat TODO
    @param nmodes TODO
    @return: TODO

    """
    nsteps = pos_mat.shape[-1]
    mode_arr = np.zeros((nmodes, 3, nsteps))
    for i in range(nsteps):
        mode_arr[:, :, i] = get_rouse_modes_at_t(pos_mat[:, :, i], nmodes)

    return mode_arr


def get_rouse_mode_corr(mode_mat):
    """Get the autocorrelation function of rouse modes.

    @param mode_mat TODO
    @return: TODO

    """
    nsteps = mode_mat.shape[-1]
    nmodes = mode_mat.shape[0]
    mode_corr = np.zeros((nmodes, nsteps))
    print(mode_corr.shape)

    for t in range(nsteps):
        for j in range(nsteps - t):
            mode_corr[:, t] += np.einsum('ij,ij->i', mode_mat[:, :, t + j],
                                         mode_mat[:, :, j])
        mode_corr[:, t] /= (nsteps - t)

    return mode_corr


def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i


def get_rouse_mode_corr_fast(mode_mat):
    """Get the autocorrelation function of rouse modes using fftw.

    @param mode_mat TODO
    @return: TODO

    """
    nsteps = mode_mat.shape[-1]
    n = next_pow_two(nsteps)

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fftn(mode_mat, s=[2 * n], axes=[-1])
    mode_corr = np.fft.ifftn(np.einsum('ijk,ijk->ik', f,
                                       np.conjugate(f)),
                             axes=[-1])[:nsteps].real

    mode_corr /= 4 * n
    return mode_corr


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
