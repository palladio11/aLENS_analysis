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

from .chrom_analysis import (get_energy_arrays, total_distr_hists,
                             get_all_rog_stats)


def graph_link_energy_vs_time(ax, h5_data, ss_frac=.75):
    """TODO: Docstring for plot_energy_vs_time.

    @param h5_data TODO
    @return: TODO

    """
    mean_energy, sem_energy, expt_energy = get_energy_arrays(h5_data)
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


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
