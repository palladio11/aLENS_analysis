#!/usr/bin/env python

"""@package docstring
File: chrom_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
import yaml
from copy import deepcopy
import gc

# Data manipulation
import numpy as np
import scipy.stats as stats
from scipy.signal import savgol_filter
from scipy import fftpack
import torch

import numpy as np

import alens_analysis as aa
from alens_analysis.helpers import gen_id


def avg_dist_from_poly_com(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(dim=0).to(device)
    tcom_dist = torch.norm(tcom_arr-pol_com, dim=1)
    tcom_dist_avg = torch.mean(tcom_dist, dim=1)
    return tcom_dist_avg


def poly_bead_msd(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    tcom_arr -= tcom_arr.mean(axis=0)
    n = com_arr.shape[0]
    Ttot = com_arr.shape[-1]
    msd = torch.zeros(Ttot, device=device)
    for i in range(1, Ttot):
        tdiff_mat = tcom_arr[:, :, i:] - tcom_arr[:, :, :-i]
        msd[i] = torch.einsum('ijk,ijk->', tdiff_mat, tdiff_mat)/((Ttot-i)*n)

    return msd


def poly_autocorr(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    Ttot = com_arr.shape[-1]
    pol_com = tcom_arr.mean(axis=0).to(device)
    autocorr = torch.zeros(Ttot, device=device)
    autocorr[0] = torch.einsum('ijk,ijk->ik',
                               (tcom_arr[:, :, :] - pol_com[:, :]),
                               (tcom_arr[:, :, :] - pol_com[:, :])).mean()
    for i in list(range(1, Ttot)):
        autocorr[i] = torch.einsum('ijk,ijk->ik',
                                   (tcom_arr[:, :, i:] - pol_com[:, i:]),
                                   (tcom_arr[:, :, :-i] - pol_com[:, :-i])).mean()

    return autocorr


def poly_autocorr_fast(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(dim=0).to(device)
    nsteps = tcom_arr.shape[-1]

    # Compute the FFT and then (from that) the auto-correlation function
    f = torch.fft.fftn(tcom_arr-pol_com, dim=[-1], norm='forward')
    power_spec = torch.einsum('ijk,ijk->ik', f, torch.conj(f))
    autocorr = torch.fft.ifftn(
        power_spec, norm='forward', dim=[-1])[:, :nsteps].real
    return autocorr.mean(dim=0)


def poly_dist_autocorr_fast(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(dim=0).to(device)
    nsteps = tcom_arr.shape[-1]
    tcom_dist = torch.norm(tcom_arr-pol_com, dim=1)

    # Compute the FFT and then (from that) the auto-correlation function
    f = torch.fft.fftn(tcom_dist, dim=[-1], norm='forward')
    power_spec = torch.einsum('ik,ik->ik', f, torch.conj(f))
    autocorr = torch.fft.ifftn(
        power_spec, norm='forward', dim=[-1])[:, :nsteps].real
    return autocorr.mean(dim=0)


def sep_autocorr(com_arr, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    tsep_mat = (tcom_arr[:, None, :, :] -
                tcom_arr[None, :, :, :]).norm(dim=2).to(device)
    n = tsep_mat.shape[0]
    Ttot = tsep_mat.shape[-1]
    tcorr_d = torch.zeros(tsep_mat.shape[-1], device=device)
    avg_tsep_mat = tsep_mat.mean(dim=(0, 1)).to(device)
    avg_sep = avg_tsep_mat.mean()
    for i in range(1, tsep_mat.shape[-1]):
        tcorr_d[i] = ((tsep_mat[:, :, i:] - avg_tsep_mat[i:]) * (tsep_mat[:,
                                                                          :, :-i]-avg_tsep_mat[:-i])).sum()/(n*n*avg_sep*avg_sep*(Ttot-i))

    return tcorr_d


def sep_autocorr_fast(com_arr, device='cpu'):
    # Create necessary matrices for analysis
    tcom_arr = torch.from_numpy(com_arr).to(device)
    nsteps = tcom_arr.shape[-1]
    tsep_mat = (tcom_arr[:, None, :, :] -
                tcom_arr[None, :, :, :]).norm(dim=2).to(device)

    # Transform into frequency-space
    f = torch.fft.fftn(tsep_mat, dim=[-1], norm='forward')

    # Clean up large arrays to prevent running out of memory
    del tcom_arr
    del tsep_mat
    torch.cuda.empty_cache() if device == 'cuda' else gc.collect()

    # Power spectrum is the product separation matrs Fourier transform
    power_spec = torch.einsum('ijk,ijk->ijk', f, torch.conj(f))

    # Transform back into time-space
    autocorr = torch.fft.ifftn(
        power_spec, norm='forward', dim=[-1])[:, :, :nsteps].real

    # Clean up power spectrum to prevent running out of memory
    del power_spec
    torch.cuda.empty_cache() if device == 'cuda' else gc.collect()
    return autocorr


def power_spec(com_arr, dt, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    n = tcom_arr.size(-1)
    pol_com = tcom_arr.mean(axis=0).to(device)

    # Compute the FFT and then (from that) the power spectrum
    f = torch.fft.fftn(tcom_arr-pol_com, dim=[-1], norm='ortho')
    power_spec = dt*torch.einsum('ijk,ijk->ik',
                                 f, torch.conj(f)).mean(dim=0)

    # Compute the frequencies
    n_modes = power_spec.size(dim=0)  # Includes negative modes
    n_ps_pos_vals = int(n_modes/2)
    freq = torch.fft.fftfreq(n_modes, dt)[:n_ps_pos_vals]

    return power_spec[:n_ps_pos_vals], freq


def poly_dist_power_spec(com_arr, dt, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(dim=0).to(device)
    tcom_dist = torch.norm(tcom_arr-pol_com, dim=1)

    # Compute the FFT and then (from that) the power spectrum
    # 0 dim: bead dimension
    # 1 dim: xyz
    # 2 dim: mode dimension
    f = torch.fft.fftn(tcom_dist, dim=[-1], norm='ortho')
    power_spec = dt*torch.einsum('ik,ik->ik',
                                 f, torch.conj(f)).mean(dim=0)

    # Compute the frequencies
    n_modes = power_spec.size(dim=0)  # Includes negative modes
    n_ps_pos_vals = int(n_modes/2)
    freq = torch.fft.fftfreq(n_modes, dt)[:n_ps_pos_vals]

    return power_spec[:n_ps_pos_vals], freq


def poly_ang_power_spec(com_arr, dt, device='cpu'):
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(dim=0).to(device)
    tdir_arr = tcom_arr-pol_com
    tdir_arr /= torch.norm(tdir_arr, dim=1)[:, None, :]

    # Compute the FFT and then (from that) the power spectrum
    # 0 dim: bead dimension
    # 1 dim: xyz
    # 2 dim: mode dimension
    f = torch.fft.fftn(tdir_arr, dim=[-1], norm='forward')
    power_spec = (1./dt)*torch.einsum('ijk,ijk->ik',
                                      f, torch.conj(f)).mean(dim=0)

    # Compute the frequencies
    n_modes = power_spec.size(dim=0)  # Includes negative modes
    n_ps_pos_vals = int(n_modes/2)
    freq = torch.fft.fftfreq(n_modes, dt)[:n_ps_pos_vals]  # Remove neg modes

    return power_spec[:n_ps_pos_vals], freq


def imag_poly_response_func(com_arr, dt, kT=.0041, device='cpu'):
    """Refer to  F. Gittes, et al. PRL 1997 
    https://doi.org/10.1103/PhysRevLett.79.3286

    Parameters
    ----------
    com_arr : N X 3 X T
        _description_
    dt : int
        _description_
    beta : float, optional
        _description_, by default .0041
    device : str, optional
        _description_, by default 'cpu'

    Returns
    -------
    _type_
        _description_
    """
    beta = 1./kT
    tcom_arr = torch.from_numpy(com_arr).to(device)
    pol_com = tcom_arr.mean(axis=0).to(device)
    nsteps = tcom_arr.shape[-1]

    # Compute the FFT and then (from that) the power spectrum
    f = torch.fft.fftn(tcom_arr-pol_com, dim=[-1], norm='forward')
    power_spec = dt * torch.einsum('ijk,ijk->ik', f, torch.conj(f))
    freq_arr = torch.fft.fftfreq(nsteps, dt).to(device)
    iresp = .5*beta*freq_arr*power_spec.mean(dim=0)
    return iresp, freq_arr


def real_poly_response_func(iresp_arr):
    """Refer to  F. Gittes, et al. PRL 1997 
    https://doi.org/10.1103/PhysRevLett.79.3286

    Parameters
    ----------
    com_arr : N X 3 X T
        _description_
    dt : int
        _description_
    beta : float, optional
        _description_, by default .0041
    device : str, optional
        _description_, by default 'cpu'

    Returns
    -------
    _type_
        _description_
    """
    # TODO add check for if this is torch or numpy array
    # Discrete cosine transform
    dct_arr = fftpack.dct(iresp_arr.numpy(), norm='ortho')
    # Discrete sine transform
    dst_arr = fftpack.dst(dct_arr, norm='ortho')

    return (2/np.pi) * dst_arr
