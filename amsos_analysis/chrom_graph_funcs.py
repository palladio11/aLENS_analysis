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
from scipy.signal import savgol_filter

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
                             calc_rad_of_gyration, get_contact_kymo_data,
                             get_pos_kymo_data, get_pos_cond_data,
                             get_sep_dist_mat, get_contact_mat_analysis,
                             get_link_tension, gauss_weighted_contact,
                             get_contact_cond_data,
                             )

from .chrom_condensate_analysis import (Condensate,
                                        get_max_and_total_cond_size,
                                        gen_condensate_track_info,
                                        extract_condensates)


def make_all_condensate_graphs(h5_data, opts, overwrite=False):
    """TODO: Docstring for make_all_condensate_graphs.

    @param h5_path TODO
    @param **kwargs TODO
    @return: TODO

    """
    cond_sty = {
        "axes.titlesize": 20,
        "axes.labelsize": 24,
        "lines.linewidth": 2,
        "lines.markersize": 2,
        "xtick.labelsize": 24,
        "ytick.labelsize": 24,
        "font.size": 20,
        "font.sans-serif": 'Helvetica',
        "text.usetex": False,
        'mathtext.fontset': 'cm',
    }
    plt.style.use(cond_sty)

    if overwrite and 'analysis' in h5_data.keys():
        print('Deleting analysis')
        del h5_data['analysis']
    analysis_grp = h5_data.require_group('analysis')

    # Start and end of data arrays
    ss_ind = 1
    end_ind = -1
    start_bead = 0
    end_bead = None
    analysis_grp.attrs['timestep_range'] = [ss_ind, end_ind]

    # Basic data
    sy_dat = h5_data['raw_data']['sylinders'][start_bead:end_bead,
                                              :, ss_ind:end_ind]
    com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])
    nbeads = com_arr.shape[0]

    # Make combined position kymo graph and condensate graph
    fig1, axarr1 = plt.subplots(1, 2, figsize=(20, 8))

    if 'pos_kymo' not in analysis_grp.keys():
        time_arr, cond_hist_arr, bin_edges = get_pos_kymo_data(
            h5_data, ts_range=(ss_ind, end_ind),
            bead_range=(start_bead, end_bead),
            bins=200,
            analysis=analysis_grp)
    else:
        time_arr = h5_data['time'][ss_ind:end_ind]
        cond_hist_arr = analysis_grp['pos_kymo'][...]
        bin_edges = analysis_grp['pos_kymo_bin_edges'][...]
    bin_centers = .5 * (bin_edges[:-1] + bin_edges[1:])
    plot_pos_kymo(fig1, axarr1[0], time_arr, cond_hist_arr, bin_edges)

    # Make position condensate graph
    if 'pos_cond_edges' not in analysis_grp:
        pos_cond_edge_coords, pos_cond_num_arr = get_pos_cond_data(
            time_arr, cond_hist_arr, bin_centers, 10, bin_win=0,
            time_win=1001, analysis=analysis_grp)
    else:
        pos_cond_edge_coords = analysis_grp['pos_cond_edges'][...]
        pos_cond_num_arr = analysis_grp['pos_cond_num'][...]
    plot_condensate_kymo(axarr1[1], pos_cond_edge_coords)
    axarr1[1].set_ylim(bin_centers[0], bin_centers[-1])

    fig1.savefig(opts.analysis_dir / f'pos_kymo.png')

    # register_cmaps()
    plt.rcParams['image.cmap'] = 'YlOrRd'

    # Make average hic plot
    contact_mat = None
    if 'avg_contact_mat' not in analysis_grp:
        log_avg_contact_mat, contact_mat, contact_kymo = get_contact_mat_analysis(
            com_arr, avg_block_step=1, analysis=analysis_grp)
    else:
        log_avg_contact_mat = analysis_grp['avg_contact_mat'][...]

    fig2, ax2 = make_hic_plot(com_arr, log_avg_contact_mat, vmin=-7)
    fig2.savefig(opts.analysis_dir / f'average_log_contatct.png')

    # Make contact kymograph and last image
    fig4, axarr4 = plt.subplots(1, 3, figsize=(24, 6))
    if 'contact_kymo' not in analysis_grp:
        if contact_mat is None:
            log_avg_contact_mat, contact_mat, contact_kymo = get_contact_mat_analysis(
                com_arr, avg_block_step=1)
    else:
        contact_kymo = analysis_grp['contact_kymo'][...]
    # TODO better handling of this
    if contact_mat is not None:
        fig3, axarr3 = make_summed_contact_kymo_graph(
            contact_mat, time_arr, contact_type="", vmin=-25, vmax=7,
            avg_contact_mat=log_avg_contact_mat, avg_vmin=-7)

    plot_contact_kymo(fig4, axarr4[0], time_arr, contact_kymo, vmax=7)

    # Make contact condensate analysis
    if 'contact_cond_edges' not in analysis_grp:
        contact_cond_edges, contact_cond_num = get_contact_cond_data(
            time_arr, contact_kymo, 3.5, bead_win=101, time_win=1001,
            analysis=analysis_grp)
    else:
        contact_cond_edges = analysis_grp['contact_cond_edges'][...]
        contact_cond_num = analysis_grp['contact_cond_num'][...]

    plot_condensate_kymo(axarr4[1], contact_cond_edges,
                         ylims=[start_bead, nbeads],
                         ylabel='Bead index')

    if 'max_contact_cond_size' not in analysis_grp:
        max_width_arr, total_width_arr = get_max_and_total_cond_size(
            time_arr, contact_cond_edges, contact_cond_num,
            analysis=analysis_grp)
    else:
        max_contact_cond_size = analysis_grp['max_contact_cond_size'][...]
        total_contact_cond_beads = analysis_grp['total_contact_cond_beads'][...]

    plot_condensate_characterize(axarr4[2], time_arr,
                                 max_contact_cond_size,
                                 total_contact_cond_beads,
                                 contact_cond_num)
    fig4.tight_layout()
    fig4.savefig(opts.analysis_dir / f'contact_cond_charact.png')

    fig5, axarr5 = plt.subplots(1, 2, figsize=(21, 6))
    if 'condensates' not in analysis_grp:
        cond_grp = analysis_grp.create_group('condensates')
        cond_lst = gen_condensate_track_info(h5_data, cond_grp)
    else:
        cond_lst = extract_condensates(analysis_grp['condensates'])

    plot_condensate_kymo(axarr5[0], contact_cond_edges,
                         ylims=[start_bead, nbeads],
                         ylabel='Bead index')

    plot_condensate_tracks(axarr5[1], time_arr, cond_lst,
                           ylims=[start_bead, nbeads])
    fig5.tight_layout()
    fig5.savefig(opts.analysis_dir / f'contact_cond_track.png')

    plt.rcParams['image.cmap'] = 'coolwarm'
    # Make tension kymograph
    fig6, ax6 = make_tension_kymo(h5_data, ss_ind, end_ind, time_win=1001)
    fig6.savefig(opts.analysis_dir / f'tension_kymo.png')


def make_contact_condensate_characterize_graphs(
        contact_kymo, time_arr, threshold, bead_win, time_win):
    """TODO: Docstring for make_contact_condensate_characterize_graphs.

    @param contact_kymo TODO
    @param time_arr TODO
    @return: TODO

    """
    fig, axarr = plt.subplots(1, 3, figsize=(20, 8))
    cond_edge_coords, cond_num_arr = get_contact_cond_data(
        time_arr, contact_kymo, threshold, bead_win, time_win)
    plot_condensate_kymo(axarr[0], cond_edge_coords, ylabel='Beads')

    return fig, axarr


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


def make_hic_plot(com_arr, log_contact_avg, vmin=-7, vmax=None):
    fig, ax = plt.subplots(figsize=(10, 8))

    nbeads = com_arr.shape[0]
    x = np.arange(nbeads + 1)[::int((nbeads) / log_contact_avg.shape[0])]
    X, Y = np.meshgrid(x, x)
    c = ax.pcolorfast(X, Y, log_contact_avg, vmax=vmax, vmin=vmin)
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


def make_summed_contact_kymo_graph(contact_mat,
                                   time_arr, contact_type="",
                                   vmin=-25, vmax=10,
                                   avg_contact_mat=None, avg_vmin=-7):

    plt.rcParams['image.cmap'] = 'YlOrRd'
    nbeads = contact_mat.shape[0]
    if isinstance(avg_contact_mat, np.ndarray):
        fig, axarr = plt.subplots(1, 3, figsize=(22, 6))
        x = np.arange(nbeads + 1)[::int((nbeads) / avg_contact_mat.shape[0])]
        X, Y = np.meshgrid(x, x)
        c = axarr[2].pcolorfast(X, Y, avg_contact_mat, vmin=avg_vmin)
        axarr[2].set_aspect('equal')
        _ = fig.colorbar(c, ax=axarr[2], label="Log contact probability")

        _ = axarr[2].set_title('Time average contact matrix')
        _ = axarr[2].set_xlabel("Bead index")
        _ = axarr[2].set_ylabel("Bead index")
    else:
        fig, axarr = plt.subplots(1, 2, figsize=(16, 6))

    # Last frame
    x = np.arange(nbeads + 1)[::int((nbeads) / contact_mat.shape[0])]
    X, Y = np.meshgrid(x, x)
    c1 = axarr[0].pcolorfast(X, Y,
                             np.log(contact_mat[:, :, -1]), vmin=vmin)
    axarr[0].set_aspect('equal')
    _ = fig.colorbar(c1, ax=axarr[0], label="Log contact probability")
    _ = axarr[0].set_xlabel(r'Bead index')
    _ = axarr[0].set_ylabel(r'Bead index')
    _ = axarr[0].set_title(
        'Last frame contact mat \n 1 bead $\sim$ 200-400 bp')

    # Contact kymograph
    contact_kymo = get_contact_kymo_data(contact_mat)
    plot_contact_kymo(fig, axarr[1], time_arr, contact_kymo,
                      contact_type=contact_type, vmax=vmax)
    # contact_kymo = np.sum(contact_mat, axis=0) - 1
    # y = np.arange(contact_kymo.shape[0] + 1)
    # # Add extra time point
    # x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    # X, Y = np.meshgrid(x, y)
    # if contact_type == "log":
    #     c0 = axarr[1].pcolorfast(X, Y, np.log(contact_kymo))
    #     _ = fig.colorbar(
    #         c0,
    #         ax=axarr[1],
    #         label="Log sum contact \n probability")
    # else:
    #     c0 = axarr[1].pcolorfast(X, Y, contact_kymo, vmax=vmax)
    #     _ = fig.colorbar(c0, ax=axarr[1],
    #                      label=r"Contact probability (sec$^{-1}$)")
    # _ = axarr[1].set_title(
    #     "Contact probabilty 'kymograph'")
    # axarr[1].set_xlabel("Time $t$ (sec)")
    # axarr[1].set_ylabel("Bead index")

    fig.tight_layout()

    return fig, axarr


def make_tension_kymo(h5_data, ss_ind, end_ind, time_win=1001):
    """TODO: Docstring for make_tension_kymo.

    @param h5_data TODO
    @return: TODO

    """
    time_arr = h5_data['time'][ss_ind:end_ind]
    tension_arr = get_link_tension(h5_data)[:, ss_ind:end_ind]
    tension_arr = savgol_filter(tension_arr, time_win, 3, axis=-1)
    fig, ax = plt.subplots(figsize=(10, 8), )
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    y = np.arange(tension_arr.shape[0] + 1)
    X, Y = np.meshgrid(x, y)
    c0 = ax.pcolorfast(X, Y, tension_arr)
    _ = fig.colorbar(c0, ax=ax, label=r"Tension (pN)")
    return fig, ax


# def make_pos_kymo_graph(h5_data, pos_hist_kymo, vmax=60):
#     """TODO: Docstring for make_pos_kymo_graph.
#     @return: TODO

#     """
#     time_arr, cond_hist_arr, bin_edges = get_pos_kymo_data(
#         h5_data, ts_range=(ss_ind, end_ind), bins=200)
#     bin_centers = .5 * (bin_edges[:-1] + bin_edges[1:])
#     pass


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


def plot_link_energy_vs_time(ax, h5_data, ss_frac=.75):
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


def plot_condensate_kymo(ax, edge_coords, xlims=None,
                         ylims=None, axis_label_flag=True, ylabel=''):
    """TODO: Docstring for plot_contact_conden_kymo.

    @param ax TODO
    @param edge_coords TODO
    @return: None

    """
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    if len(edge_coords) > 0:
        ax.scatter(edge_coords[:, 0], edge_coords[:, 1], label='start')
        ax.scatter(edge_coords[:, 0], edge_coords[:, 2], label='end')
        ax.scatter(edge_coords[:, 0], .5 *
                   (edge_coords[:, 2] + edge_coords[:, 1]), label='center')
        ax.vlines(edge_coords[:, 0], edge_coords[:, 1],
                  edge_coords[:, 2], color='k', alpha=.01)
    # ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.05))
    if axis_label_flag:
        ax.set_xlabel("Time $t$ (sec)")
        ax.set_ylabel(ylabel)


def plot_condensate_characterize(
        ax0, time_arr, max_cond_size, tot_cond_beads, cond_num_arr):
    """TODO: Docstring for plot_condensate_characterize.

    @param ax TODO
    @param edge_coords TODO
    @param cond_num_arr TODO
    @return: TODO

    """
    ax1 = ax0.twinx()

    ax0.set_xlabel('Time $t$ (sec)')

    ax0.plot(time_arr, tot_cond_beads, color='b', label='Total')
    ax0.plot(time_arr,
             max_cond_size,
             color='b',
             linestyle='--',
             label='Largest')
    ax0.set_ylabel('Bead number', color='b')
    ax0.tick_params(axis='y', labelcolor='b')
    ax0.legend()
    ax0.set_ylim(0)

    ax1.plot(time_arr, cond_num_arr, color='r')
    ax1.set_ylabel('Condensate number', color='r')
    ax1.set_ylim(0)
    ax1.tick_params(axis='y', labelcolor='r')
    plt.tight_layout()


def plot_contact_kymo(fig, ax, time_arr, contact_kymo,
                      contact_type="", vmax=10):
    """TODO: Docstring for plot_contact_kymo.

    @param ax TODO
    @param contact_mat TODO
    @return: TODO

    """
    y = np.arange(contact_kymo.shape[0] + 1)
    # Add extra time point
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)
    if contact_type == "log":
        c = ax.pcolorfast(X, Y, np.log(contact_kymo))
        _ = fig.colorbar(c, ax=ax, label="Log sum contact \n probability")
    else:
        c = ax.pcolorfast(X, Y, contact_kymo, vmax=vmax)
        _ = fig.colorbar(c, ax=ax, label=r"Contact probability")
    _ = ax.set_title("Contact probabilty 'kymograph'")
    ax.set_xlabel("Time $t$ (sec)")
    ax.set_ylabel("Bead index")
    return


def plot_pos_kymo(fig, ax, time_arr, pos_hist_kymo, bin_edges, vmax=60):
    """TODO: Docstring for plot_contact_kymo.

    @param ax TODO
    @param contact_mat TODO
    @return: TODO

    """
    y = bin_edges
    # Add extra time point
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)

    ax.set_title("Position kymograph")
    c0 = ax.pcolorfast(X, Y, pos_hist_kymo, vmax=vmax)
    fig.colorbar(c0, ax=ax, label=r"Number of beads")
    ax.set_ylabel(r"Position $x$ ($\mu$m)")
    ax.set_xlabel("Time $t$ (sec)")


def plot_condensate_tracks(ax, time_arr, cond_lst, ylims=None):
    """TODO: Docstring for plot_condesate_tracks.

    @param ax TODO
    @param time_arr TODO
    @param cond_lst TODO
    @param ylims TODO
    @return: TODO

    """
    for cond in cond_lst:
        ax.plot(cond.time_arr,
                .5 * np.asarray(cond.edge_coord_arr).sum(axis=1),
                label=f'id = {cond.id}')
    ax.legend(loc='center left', bbox_to_anchor=(1.05, .5))
    ax.set_xlim(time_arr[0], time_arr[-1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])

    ax.set_ylabel('Bead index')
    ax.set_xlabel('Time $t$ (sec)')


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
