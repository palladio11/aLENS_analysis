#!/usr/bin/env python

"""@package docstring
File: jn_scratch.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import h5py
import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def get_energy_arrays(h5_data, write=False):
    """TODO: Docstring for get_mean_energy_array.

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
    if write:
        energy_dset = h5_data['analysis'].create_dataset(
            'energy', data=np.stack(mean_energy, sem_energy))
        energy_dset.attrs['nsylinders'] = energy_arr.shape[0]
    return mean_energy, sem_energy, kbt


fig, ax = plt.subplots(1, 1, figsize=(10, 6))


def plot_energy_vs_time(data_path):
    """TODO: Docstring for plot_energy_vs_time.

    @param data_path TODO
    @return: TODO

    """
    with h5py.File(next(data_path.glob('*.h5')), 'r+') as h5_data:
        mean_energy, std_energy, kbt = get_energy_arrays(h5_data)
        time = h5_data['time'][...]

        # print(sy_dat)
        # end_end_sep = np.linalg.norm(sy_dat[-1,5:8,:] - sy_dat[0,5:8,:],axis=0)
        # print(end_end_sep)
        ax.plot(time, mean_energy, label="Self-avoidng")
        energy_mean = mean_energy[int(mean_energy.size * .75):].mean()
        print(energy_mean)
        ax.axhline(energy_mean, color='r')
        ax.set_ylabel(r"Mean spring energy (pN$\cdot$nm)")
        ax.set_xlabel("Time (sec)")


def plot_diff_hist_vs_time(data_path):
    """TODO: Docstring for plot_diff_hist_vs_time.

    @param data_path TODO
    @return: TODO

    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    with h5py.File(next(data_path.glob('*.h5')), 'r+') as h5_data:
        diff_hist, bin_edges = get_sep_hist(h5_data)
        time = h5_data['time'][...]

        c = ax.pcolorfast(bin_edges, time, diff_hist)
        fig.colorbar(c, ax=ax, label=r'Pairs')
        ax.set_xlabel(r"Bead COM separation ($\mu$m)")
        ax.set_ylabel("Time (sec)")


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


def get_sep_dist_mat(h5_data, ss_ind=0):
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
    avg_overlap_arr = overlap_dist_mat.sum(axis=(0, 1)) / num_overlap
    min_overlap_arr = overlap_dist_mat.min(axis=(0, 1))

    return num_overlap, avg_overlap_arr, min_overlap_arr


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

    @param sphere_dat TODO
    @param nmodes TODO
    @return: TODO

    """
    nsteps = pos_mat.shape[-1]
    mode_arr = np.zeros((nmodes, 3, nsteps))
    for i in range(nsteps):
        mode_arr[:, :, i] = get_rouse_modes_at_t(pos_mat[:, :, i], nmodes)

    return mode_arr


def get_rouse_mode_corr(mode_mat):
    """TODO: Docstring for get_rouse_modes.

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


def distr_hists(pos_mat, free_frac_chain=.5, rel_ind=0, nbins=100,):
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
    dist_hist, dist_bin_edges = np.histogram(dist_arr, nbins)
    z_rho_hist, rho_bin_edges, z_bin_edges = np.histogram2d(
        np.linalg.norm(rel_vec_arr[:-1, :], axis=0), rel_vec_arr[-1, :], nbins)

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
        dist_arr, nbins, range=[
            0, hist_max], density=True)
    z_rho_hist, rho_bin_edges, z_bin_edges = np.histogram2d(
        np.linalg.norm(rel_vec_arr[:, :-1, :], axis=1).flatten(),
        rel_vec_arr[:, -1, :].flatten(), int(nbins / 2),
        range=[[0, hist_max], [-hist_max, hist_max]], density=True)

    return ((dist_hist, dist_bin_edges),
            (z_rho_hist, rho_bin_edges, z_bin_edges))


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
