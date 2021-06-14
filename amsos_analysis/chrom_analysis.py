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


def gauss_weighted_contact(sep_mat, sigma=.020):
    return np.exp(-np.power(sep_mat, 2) / (2. * (sigma * sigma)))


def log_gauss_weighted_contact(sep_mat, sigma=.020):
    return -np.power(sep_mat, 2) / (2. * (sigma * sigma)) / np.log(10)


def get_link_energy_arrays(h5_data, write=False):
    """ Get the mean, standard deviation, and expected energy of all links in
    a bead-spring chain

    @param h5_data HDF5 data file to analyze with all raw data about filaments
    @param write If true, will write data directly to the analysis group in
                 the h5_data file.
    @return: TODO

    """
    sy_dat = h5_data['raw_data']['sylinders'][...]
    params = yaml.safe_load(h5_data.attrs['RunConfig'])
    k_spring = params['linkKappa']
    kbt = params['KBT']

    rest_length = params['linkGap'] + sy_dat[1:, 1, :] + sy_dat[:-1, 1, :]
    sep_vec = sy_dat[1:, 2:5, :] - sy_dat[:-1, 5:8, :]

    sep_mag = np.linalg.norm(sep_vec, axis=1)

    energy_arr = .5 * k_spring * np.power(sep_mag - rest_length, 2)
    mean_energy = np.mean(energy_arr, axis=0)
    sem_energy = stats.sem(energy_arr, axis=0)
    expt_energy = kbt * \
        (.5 - 1. /
         (1. + (k_spring * rest_length[0, 0] * rest_length[0, 0] / kbt)))
    if write:
        energy_dset = h5_data['analysis'].create_dataset(
            'link_energy', data=np.stack(mean_energy, sem_energy))
        energy_dset.attrs['nsylinders'] = energy_arr.shape[0]
    return mean_energy, sem_energy, kbt


def get_sep_hist(h5_data, nbins=100, ss_ind=0, write=False):
    """Returns a 2D histogram of bead separations vs time

    @param h5_data TODO
    @return: TODO

    """
    params = yaml.safe_load(h5_data.attrs['RunConfig'])
    hist_min = params['sylinderDiameter'] * .8
    hist_max = params['sylinderDiameter'] * 1.2

    dist_hist = []
    dist_mat = get_sep_dist_mat(h5_data, ss_ind)

    for i in range(dist_mat.shape[-1]):
        hist, bin_edges = np.histogram(
            dist_mat[:, :, i].flatten(), nbins, range=(hist_min, hist_max))
        dist_hist += [hist * .5]

    return dist_hist, bin_edges


def get_sep_dist_mat(h5_data, ss_ind=0, write=False):
    """Returns a NxNxM matrix of NXN filaments distances over M time points
    starting at ss_ind time point.

    @param h5_data TODO
    @return: TODO

    """
    sy_dat = h5_data['raw_data']['sylinders'][...]

    com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    dist_mat = np.linalg.norm((com_arr[:, np.newaxis, :, ss_ind:] -
                               com_arr[np.newaxis, :, :, ss_ind:]),
                              axis=2)

    return dist_mat


def get_overlap_arrs(dist_mat, sy_diam):
    """Returns a NxNxM matrix of NXN filaments distances over M time points
    starting at ss_ind time point.

    @param h5_data TODO
    @return: TODO

    """
    is_overlap_mat = (dist_mat < sy_diam).astype(int)
    num_overlap = .5 * (is_overlap_mat.sum(axis=(0, 1))
                        - dist_mat.shape[0])  # remove self-overlap
    overlap_dist_mat = np.einsum('ijk, ijk -> ijk', dist_mat, is_overlap_mat)
    avg_overlap_arr = (.5 * overlap_dist_mat.sum(axis=(0, 1))) / num_overlap
    min_overlap_arr = np.ma.masked_values(overlap_dist_mat, 0).min(axis=(0, 1))

    return num_overlap, avg_overlap_arr, min_overlap_arr


def autocorr_bead_pos(com_arr, ignore_id=None):
    """Find the autocorrelation function for bead positions

    @param com_arr TODO
    @param ignore_id TODO
    @return: TODO

    """

    #com_rel_arr = com_arr - com_arr.mean(axis=-1)[:, :, np.newaxis]
    com_rel_arr = com_arr[...]

    if ignore_id is not None:
        com_rel_arr = np.delete(com_rel_arr, ignore_id, axis=0)
    nsteps = com_rel_arr.shape[-1]
    nbeads = com_rel_arr.shape[0]
    auto_corr = np.zeros((nbeads, nsteps))
    for t in range(nsteps):
        for j in range(nsteps - t):
            auto_corr[:, t] += np.einsum('ij,ij->i', com_rel_arr[:, :, t + j],
                                         com_rel_arr[:, :, j])
        auto_corr[:, t] /= (nsteps - t)

    return auto_corr


def distr_hists(pos_mat, free_frac_chain=.5,
                rel_ind=0, nbins=100, hist_max=1.):
    """TODO: Docstring for radial_distr.
    @param pos_mat TODO
    @param free_frac_chain TODO
    @param rel_ind TODO
    @param nbins TODO
    @return: TODO
    """
    nbeads = pos_mat.shape[0]
    ind = int(nbeads * free_frac_chain)

    rel_vec_arr = pos_mat[ind, :, :] - pos_mat[rel_ind, :, :]
    dist_arr = np.linalg.norm(rel_vec_arr, axis=0)

    dist_hist, dist_bin_edges = np.histogram(
        dist_arr, nbins, range=[0, hist_max], density=True)
    z_rho_hist, rho_bin_edges, z_bin_edges = np.histogram2d(
        np.linalg.norm(rel_vec_arr[:-1, :], axis=0), rel_vec_arr[-1, :],
        int(nbins / 2), range=[[0, hist_max], [-hist_max, hist_max]], density=True)

    return ((dist_hist, dist_bin_edges),
            (z_rho_hist, rho_bin_edges, z_bin_edges))


def total_distr_hists(pos_mat, rel_ind=0, nbins=100, hist_max=1):
    """TODO: Docstring for radial_distr.
    @param pos_mat TODO
    @param free_frac_chain TODO
    @param rel_ind TODO
    @param nbins TODO
    @return: TODO
    """
    rel_vec_arr = pos_mat - (pos_mat[rel_ind])[np.newaxis, :, :]
    dist_arr = np.linalg.norm(rel_vec_arr, axis=1).flatten()

    dist_hist, dist_bin_edges = np.histogram(
        dist_arr, nbins, range=[0, hist_max], density=True)
    z_rho_hist, rho_bin_edges, z_bin_edges = np.histogram2d(
        np.linalg.norm(rel_vec_arr[:, :-1, :], axis=1).flatten(),
        rel_vec_arr[:, -1, :].flatten(), int(nbins / 2),
        range=[[0, hist_max], [-hist_max, hist_max]], density=True)

    return ((dist_hist, dist_bin_edges),
            (z_rho_hist, rho_bin_edges, z_bin_edges))


def cart_distr_hists(pos_mat, rel_pos, e0_ind, e1_ind, nbins=100, hist_max=1.):
    """TODO: Docstring for radial_distr.
    @param pos_mat TODO
    @param rel_ind TODO
    @param nbins TODO
    @return: TODO
    """
    rel_vec_arr = pos_mat - (rel_pos)[np.newaxis, :, :]
    e0_e1_hist, e0_edges, e1_edges = np.histogram2d(
        rel_vec_arr[:, e0_ind, :].flatten(),
        rel_vec_arr[:, e1_ind, :].flatten(),
        int(nbins / 2),
        range=[[-hist_max, hist_max], [-hist_max, hist_max]], density=True)

    return (e0_e1_hist, e0_edges, e1_edges)


def cylin_distr_hists(pos_mat, zero_pos, z_uvec, nbins=100, hist_max=1.):
    """TODO: Docstring for cylindrical histogram.
    @param pos_mat TODO
    @param rel_ind TODO
    @param nbins TODO
    @return: TODO
    """
    rel_vec_arr = pos_mat - (zero_pos)[np.newaxis, :, :]
    z_proj_arr = np.einsum('ijk,jk->ik', rel_vec_arr, z_uvec)
    rho_proj_arr = np.linalg.norm(
        rel_vec_arr - np.einsum('jk,ik->ijk', z_uvec, z_proj_arr), axis=1)
    rho_z_hist, rho_bin_edges, z_bin_edges = np.histogram2d(
        rho_proj_arr.flatten(), z_proj_arr.flatten(), int(nbins / 2),
        range=[[0, hist_max], [-hist_max, hist_max]], density=True)

    return (rho_z_hist, rho_bin_edges, z_bin_edges)


def rad_distr_hists(pos_mat, zero_pos, nbins=100, hist_max=1.):
    """TODO: Docstring for cylindrical histogram.
    @param pos_mat TODO
    @param nbins TODO
    @return: TODO
    """
    rel_vec_arr = pos_mat - (zero_pos)[np.newaxis, :, :]
    rad_arr = np.linalg.norm(rel_vec_arr, axis=1).flatten()

    rad_hist, rad_bin_edges = np.histogram(
        rad_arr, nbins, range=[
            0, hist_max], density=True)
    return (rad_hist, rad_bin_edges)


def get_all_rog_stats(pos_mat, rel_ind=0):
    rel_vec_arr = pos_mat - (pos_mat[rel_ind])[np.newaxis, :, :]
    pos_avg_arr = rel_vec_arr.mean(axis=2)
    pos_std_arr = rel_vec_arr.std(axis=2)
    rad_pos_arr = np.linalg.norm(pos_avg_arr, axis=1)
    rog_arr = np.linalg.norm(pos_std_arr, axis=1)

    return(pos_avg_arr, pos_std_arr, rad_pos_arr, rog_arr)


def get_time_avg_contact_mat(com_arr, sigma=.02, avg_block_step=1):
    reduc_com_arr = com_arr[::avg_block_step, :, :]  # simple downsampling
    sep_mat = np.linalg.norm(
        reduc_com_arr[:, np.newaxis, :, :] - reduc_com_arr[np.newaxis, :, :, :], axis=2)
    log_contact_map = log_gauss_weighted_contact(sep_mat, sigma)
    return log_contact_map.mean(axis=-1)


def get_end_end_distance(com_arr):
    return np.linalg.norm(com_arr[0, :, :] - com_arr[-1, :, :], axis=0)


def calc_rad_of_gyration(pos_mat):
    """Calculate the radius of gyration of filament

    @param pos_mat TODO
    @return: TODO

    """
    n_beads = float(pos_mat.shape[0])
    rel_pos_arr = pos_mat - np.mean(pos_mat, axis=0)

    rog_sqr_arr = np.einsum('ijk,ijk->k', rel_pos_arr, rel_pos_arr) / n_beads
    return np.sqrt(rog_sqr_arr)


def find_neighbors(com_arr, diam, time_ind=0):
    """Find beads that are in close proximity with one another at any given time.

    """
    neighbor_mat = (np.linalg.norm((com_arr[:, np.newaxis, :, time_ind] -
                                    com_arr[np.newaxis, :, :, time_ind]),
                                   axis=2) < diam * 1.2).astype(int)
    return neighbor_mat


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
