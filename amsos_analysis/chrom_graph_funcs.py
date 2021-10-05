#!/usr/bin/env python

"""@package docstring
File: chrom_graph_funcs.py
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

# Visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import (Circle,
                                RegularPolygon,
                                FancyArrowPatch,
                                ArrowStyle)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, NullFormatter)
import matplotlib.colors as colors

from .chrom_analysis import (get_link_energy_arrays, total_distr_hists,
                             get_all_rog_stats, cart_distr_hists,
                             cylin_distr_hists, rad_distr_hists,
                             calc_rad_of_gyration)


def graph_link_energy_vs_time(ax, h5_data, ss_frac=.75):
    """TODO: Docstring for plot_energy_vs_time.

    @param h5_data TODO
    @return: TODO

    """
    mean_energy, sem_energy, expt_energy = get_link_energy_arrays(h5_data)
    time = h5_data['time'][...]

    ax.plot(time, mean_energy)
    ax.fill_between(time,
                    (mean_energy - sem_energy),
                    (mean_energy + sem_energy),
                    color='blue',
                    alpha=0.1)

    energy_mean = mean_energy[int(mean_energy.size * ss_frac):].mean()

    ax.axhline(energy_mean, color='r')
    ax.set_ylabel(r"Mean spring energy (pN$\cdot \mu$m)")
    ax.set_xlabel("Time (sec)")


def make_total_distr_plots(com_arr, log_contact_avg=None, hist_max=1.,
                           rel_ind=0, vmin=-50):
    fig, axarr = plt.subplots(2, 3, figsize=(26, 16))
    ax = axarr.flatten()
    ax[1].set_aspect('equal')
    dist_hist_dat, z_rho_hist_dat = total_distr_hists(
        com_arr, hist_max=hist_max, rel_ind=rel_ind)
    ax[0].bar(dist_hist_dat[1][:-1], dist_hist_dat[0],
              width=np.diff(dist_hist_dat[1]), align='edge')

    X, Y = np.meshgrid(z_rho_hist_dat[1], z_rho_hist_dat[2])
    ax[1].pcolorfast(X, Y, z_rho_hist_dat[0].T)
    pos_avg_arr, pos_std_arr, rad_pos_arr, rog_arr = get_all_rog_stats(
        com_arr, rel_ind=rel_ind)
    ax[2].scatter(np.arange(rog_arr.size), rad_pos_arr)
    ax[3].scatter(np.arange(rog_arr.size), rog_arr)

    if isinstance(log_contact_avg, np.ndarray):
        nbeads = com_arr.shape[0]
        x = np.arange(nbeads + 1)[::int(nbeads / log_contact_avg.shape[0])]
        X, Y = np.meshgrid(x, x)
        c = ax[4].pcolorfast(X, Y, log_contact_avg, vmin=vmin)
        fig.colorbar(c, label="Log contact probability")

    ax[0].set_xlabel(r'$|{\bf r} - {\bf r}_0|$ ($\mu$m)')
    ax[0].set_ylabel(r'Probability density ($\mu$m$^{-1}$)')
    ax[1].set_xlabel(r'$\rho - \rho_0$ ($\mu$m)')
    ax[1].set_ylabel(r'$z-z_0$ ($\mu$m)')
    ax[2].set_xlabel(r'Bead index')
    ax[2].set_ylabel(r'$|{\bf r} - {\bf r}_0|$ ($\mu$m)')
    ax[3].set_xlabel(r'Bead index')
    ax[3].set_ylabel(r'$\sigma_{|{\bf r} - {\bf r}_0|}$ ($\mu$m)')

    fig.tight_layout()

    return fig, ax


def make_min_distr_plots(com_arr, log_contact_avg=None,
                         hist_max=1., rel_ind=0, vmin=-30):
    fig, axarr = plt.subplots(1, 3, figsize=(26, 8))
    ax = axarr.flatten()
    pos_avg_arr, pos_std_arr, rad_pos_arr, rog_arr = get_all_rog_stats(
        com_arr, rel_ind=rel_ind)
    ax[0].scatter(np.arange(rog_arr.size), rad_pos_arr)
    ax[1].scatter(np.arange(rog_arr.size), rog_arr)

    if isinstance(log_contact_avg, np.ndarray):
        nbeads = com_arr.shape[0]
        x = np.arange(nbeads + 1)[::int((nbeads) / log_contact_avg.shape[0])]
        X, Y = np.meshgrid(x, x)
        c = ax[2].pcolorfast(X, Y, log_contact_avg, vmin=vmin)
        ax[2].set_aspect('equal')
        fig.colorbar(c, label="Log contact probability")

    ax[0].set_xlabel(r'Bead index')
    ax[0].set_ylabel(r'$\langle|{\bf r}_i - {\bf r}_0|^2\rangle$ ($\mu$m)')
    ax[1].set_xlabel(r'Bead index')
    ax[1].set_ylabel(
        r'$\langle \sigma_{|{\bf r}_i - {\bf r}_0|}\rangle$ ($\mu$m)')
    ax[2].set_xlabel(r'Bead index')
    ax[2].set_ylabel(r'Bead index')
    plt.tight_layout()

    return fig, ax


def make_hic_plot(com_arr, log_contact_avg, vmin=-30):
    fig, ax = plt.subplots(figsize=(10, 8))

    nbeads = com_arr.shape[0]
    x = np.arange(nbeads + 1)[::int((nbeads) / log_contact_avg.shape[0])]
    X, Y = np.meshgrid(x, x)
    c = ax.pcolorfast(X, Y, log_contact_avg, vmin=vmin)
    ax.set_aspect('equal')
    fig.colorbar(c, label="Log contact probability")

    ax.set_title(r'1 bead $\sim$ 200-400 bp')
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Bead index")
    plt.tight_layout()

    return fig, ax


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
    axarr = axarr.flatten()

    zero_pos = .5 * (com_arr[rel_ids[0]] + com_arr[rel_ids[1]])

    for i in range(n_rows):
        if i == 0:
            if n_rows == 1:
                seg_com_arr = com_arr
                seg_zpos_arr = zero_pos
            else:
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
        axarr[i * 3].pcolorfast(X, Y, e0_e1_hist.T)
        axarr[i * 3].set_xlabel(r'$x-x_0$ ($\mu$m)')

        rho_z_hist, rho_edges, z_edges = cylin_distr_hists(
            seg_com_arr, zero_pos, z_uvec, hist_max=hist_max)
        Rho, Z = np.meshgrid(rho_edges, z_edges)
        axarr[i * 3 + 1].pcolorfast(Rho, Z, rho_z_hist.T)
        axarr[i * 3 + 1].set_xlabel(r'$\rho - \rho_0$ ($\mu$m)')

        rad_hist, rad_edges = rad_distr_hists(
            seg_com_arr, zero_pos, hist_max=hist_max)
        rad_mean = np.average(
            rad_edges[:-1] + np.diff(rad_edges), weights=rad_hist)
        axarr[i * 3 + 2].bar(rad_edges[:-1], rad_hist,
                             width=np.diff(rad_edges), align='edge')
        axarr[i * 3 + 2].set_ylabel(r'Probability density ($\mu$m$^{-1}$)')
        axarr[i * 3 + 2].set_xlabel(r'$|{\bf r} - {\bf r}_0|$ ($\mu$m)')
        axarr[i * 3 + 2].axvline(rad_mean, color='r',
                                 label=f'mean = {rad_mean:.4} $\mu$m')

        axarr[i * 3 + 2].legend()

        for ax in axarr[i * 3:i * 3 + 2]:
            ax.set_aspect('equal')
            ax.set_ylabel(r'$z-z_0$ ($\mu$m)')

    plt.tight_layout()

    return fig, axarr


def make_summed_contact_kymo_graph(
        contact_mat, time_arr, contact_type="", vmin=-30, vmax=10):
    fig, axarr = plt.subplots(1, 2, figsize=(16, 6))

    nbeads = contact_mat.shape[0]
    x = np.arange(nbeads + 1)[::int((nbeads) / contact_mat.shape[0])]
    X, Y = np.meshgrid(x, x)
    c0 = axarr[0].pcolorfast(X, Y, np.log(contact_mat[:, :, -1]), vmin=vmin)
    axarr[0].set_aspect('equal')
    _ = fig.colorbar(c0, ax=axarr[0], label="Log contact probability")
    _ = axarr[0].set_xlabel(r'Bead index')
    _ = axarr[0].set_ylabel(r'Bead index')

    contact_kymo = np.sum(contact_mat, axis=0) - 1
    y = np.arange(contact_kymo.shape[0] + 1)
    # Add extra time point
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)
    if contact_type == "log":
        c1 = axarr[1].pcolorfast(X, Y, np.log(contact_kymo))
        _ = fig.colorbar(
            c1,
            ax=axarr[1],
            label="Log sum contact \n probability")
    else:
        c1 = axarr[1].pcolorfast(X, Y, contact_kymo, vmax=vmax)
        _ = fig.colorbar(c1, ax=axarr[1], label="Sum contact probability")
    axarr[1].set_xlabel("Time $t$ (sec)")
    axarr[1].set_ylabel("Bead index")
    fig.tight_layout()

    return fig, axarr, contact_kymo


def make_rog_vs_time_graph(time_arr, com_arr, label=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    rog_arr = calc_rad_of_gyration(com_arr)
    _ = ax.plot(time_arr, rog_arr, label=label)
    _ = ax.set_xlabel('Time (sec)')
    _ = ax.set_ylabel(r'Radius of gyration $R_g$ ($\mu$m)')
    return fig, ax


def plot_rog_vs_time_graph(ax, time_arr, com_arr, label=None):
    rog_arr = calc_rad_of_gyration(com_arr)
    _ = ax.plot(time_arr, rog_arr, label=label)
    return


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
