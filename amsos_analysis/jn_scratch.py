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


def plot_energy_vs_time(data_path):
    """TODO: Docstring for plot_energy_vs_time.

    @param data_path TODO
    @return: TODO

    """
    with h5py.File(next(data_path.glob('*.h5')), 'r+') as h5_data:
        mean_energy, std_energy, kbt = get_energy_arrays(h5_data)
        time = h5_data['time'][...]

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


def find_neighbors(com_arr, diam, time_ind=0):
    """TODO: Docstring for find_neighbors.

    """
    neighbor_mat = (np.linalg.norm((com_arr[:, np.newaxis, :, time_ind] -
                                    com_arr[np.newaxis, :, :, time_ind]),
                                   axis=2) < diam * 1.2).astype(int)
    return neighbor_mat


def make_segment_distr_graphs(
        com_arr, sep_ids, rel_ids, e0_ind, e1_ind, hist_max=1.):
    """TODO: Docstring for make_segment_distr_graphs.

    @param com_arr Center-of-mass array
    @param sep_ids IDs to define the boundaries of the different segments
    @param rel_ids List or tuple of gids of segements to define the zero
                   position from the average of both
    @param z_uvec The z-direction of the cartesian and cylindrical graphs
    @param hist_max absoulte maximum range of histograms
    @return: TODO

    """
    n_rows = len(sep_ids) + 1

    fig, axarr = plt.subplots(n_rows, 3, figsize=(18 + 4, n_rows * 8))

    zero_pos = .5 * (com_arr[rel_ids[0]] + com_arr[rel_ids[1]])

    for i in range(n_rows):
        if i == 0:
            seg_com_arr = com_arr[:sep_ids[i]]
            seg_zpos_arr = zero_pos[:sep_ids[i]]
        elif i == n_rows - 1:
            seg_com_arr = com_arr[sep_ids[i - 1]:]
            seg_zpos_arr = zero_pos[sep_ids[i - 1]:]
        else:
            seg_com_arr = com_arr[sep_ids[i - 1]:sep_ids[i]]
            seg_zpos_arr = zero_pos[sep_ids[i - 1]:sep_ids[i]]

        z_uvec = np.zeros((3, seg_zpos_arr.shape[-1]))
        z_uvec[e1_ind] = 1.
        e0_e1_hist, e0_edges, e1_edges = cart_distr_hists(
            seg_com_arr, zero_pos, e0_ind, e1_ind, hist_max=hist_max)
        X, Y = np.meshgrid(e0_edges, e1_edges)
        axarr[i, 0].pcolorfast(X, Y, e0_e1_hist.T)
        axarr[i, 0].set_xlabel(r'$x-x_0$ ($\mu$m)')

        rho_z_hist, rho_edges, z_edges = cylin_distr_hists(
            seg_com_arr, zero_pos, z_uvec, hist_max=hist_max)
        Rho, Z = np.meshgrid(rho_edges, z_edges)
        axarr[i, 1].pcolorfast(Rho, Z, rho_z_hist.T)
        axarr[i, 1].set_xlabel(r'$\rho - \rho_0$ ($\mu$m)')

        rad_hist, rad_edges = rad_distr_hists(
            seg_com_arr, zero_pos, hist_max=hist_max)
        rad_mean = np.sum(rad_hist * (rad_edges[:-1] + np.diff(rad_edges)))
        axarr[i, 2].bar(rad_edges[:-1], rad_hist,
                        width=np.diff(rad_edges), align='edge')
        axarr[i, 2].set_ylabel(r'Probability density ($\mu$m$^{-1}$)')
        axarr[i, 2].set_xlabel(r'$|{\bf r} - {\bf r}_0|$ ($\mu$m)')
        axarr[i, 2].axvline(rad_mean, color='r',
                            label=f'mean = {rad_mean:.4} $\mu$m')
        axarr[i, 2].legend()

        for ax in axarr[i, :-1]:
            ax.set_aspect('equal')
            ax.set_ylabel(r'$z-z_0$ ($\mu$m)')

    plt.tight_layout()


def cart_distr_hists(pos_mat, rel_pos, e0_ind=0,
                     e1_ind=0, nbins=100, hist_max=1.):
    """TODO: Docstring for radial_distr.
    @param pos_mat TODO
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


def autocorr_bead_pos(com_arr, ignore_id=None):
    """Find the autocorrelation function for bead positions

    @param com_arr TODO
    @param ignore_id TODO
    @return: TODO

    """

    com_rel_arr = com_arr - com_arr.mean(axis=-1)[:, :, np.newaxis]

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


def make_rouse_analysis_plots(
        time_arr, mode_mat, nbeads, nmodes=20, max_ind=-1):
    fig, ax = plt.subplots(1, 2, figsize=(14, 8))

    # nbeads = com_arr.shape[0]

    mode_corr = get_rouse_mode_corr(mode_mat)
    mode_arr = np.arange(nmodes)

    ax[0].loglog(mode_arr / nbeads, np.power(nbeads * mode_arr, -2.) * 1e6)
    ax[0].loglog(mode_arr / nbeads, nbeads * mode_corr[:, 0])

    for k in range(nmodes):
        ax[1].loglog(np.power(k / nbeads, 2) * time_arr[:max_ind],
                     -np.log(mode_corr[k, :max_ind] / mode_corr[k, 0]), '+')

    ax[1].loglog(time_arr[:max_ind] / nbeads, time_arr[:max_ind] / nbeads)

    ax[0].set_xlabel(r'Mode index ($k/N_b$)')
    ax[0].set_ylabel('Mode amplitude correlation\n' + r'($N_bX_{kk}(0)$)')

    ax[1].set_xlabel(r'Time (($k/N_b)^2t$)')
    ax[1].set_ylabel(r'$-\ln(X_{kk}(t)/X_{kk}(0))$)')

    fig.tight_layout()


def calc_rad_of_gyration(pos_mat):
    """Calculate the radius of gyration of filament

    @param pos_mat TODO
    @return: TODO

    """
    n_beads = float(pos_mat.shape[0])
    rel_pos_arr = pos_mat - np.mean(pos_mat, axis=0)

    rog_sqr_arr = np.einsum('ijk,ijk->k', rel_pos_arr, rel_pos_arr) / n_beads
    return np.sqrt(rog_sqr_arr)


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
