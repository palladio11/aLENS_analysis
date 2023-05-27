#!/usr/bin/env python

"""@package docstring
File: chrom_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

# Data manipulation
import numpy as np
import torch


def get_rouse_modes_at_t(pos_arr, n_modes=20):
    """TODO: Docstring for get_rouse_modes.

    @param sphere_dat TODO
    @param nmodes TODO
    @return: TODO

    """
    # modes = []
    # mode_0 = pos_arr[0]
    com_arr = pos_arr - pos_arr[0]
    nbeads = pos_arr.shape[0]
    mode_coeff_arr = np.zeros((n_modes, 3))
    kn_mat = np.einsum('n,k->nk',
                       .5+np.arange(0, nbeads),
                       (np.pi/nbeads)*np.arange(0, n_modes))
    modes = np.einsum('nk,ni->ki', np.cos(kn_mat), com_arr)

    # for k in range(nmodes):
    #     modes += [np.zeros(3)]
    #     for n in range(nbeads - 1):
    #         modes[-1] += ((pos_arr[n] - mode_0) *
    #                       np.cos(np.pi * (n + .5) * k / nbeads))
    # return np.asarray(modes) / (nbeads)
    return modes / nbeads


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


# def get_rouse_mode_corr_fast(mode_mat):
#     """Get the autocorrelation function of rouse modes using fftw.

#     @param mode_mat TODO
#     @return: TODO

#     """
#     nsteps = mode_mat.shape[-1]
#     n = next_pow_two(nsteps)

#     # Compute the FFT and then (from that) the auto-correlation function
#     f = np.fft.fftn(mode_mat, s=[2 * n], axes=[-1])
#     mode_corr = np.fft.ifftn(np.einsum('ijk,ijk->ik', f,
#                                        np.conjugate(f)),
#                              axes=[-1])[:nsteps].real

#     mode_corr /= 4 * n
#     return mode_corr

def get_rouse_mode_corr_fast(mode_mat, device='cpu'):
    """Get the autocorrelation function of rouse modes using fftw.

    @param mode_mat TODO
    @return: TODO

    """
    tmode_mat = torch.from_numpy(mode_mat).to(device)
    nsteps = tmode_mat.shape[-1]
    # n = next_pow_two(nsteps)

    # Compute the FFT and then (from that) the auto-correlation function
    f = torch.fft.fftn(tmode_mat, dim=[-1], norm='forward')
    power_spec = torch.einsum('ijk,ijk->ik', f, torch.conj(f))
    n_pos_vals = int(power_spec.size(-1)/2)
    mode_corr = torch.fft.ifftn(
        power_spec, norm='forward', dim=[-1])[:, :nsteps].real

    return mode_corr[:, :n_pos_vals]


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
