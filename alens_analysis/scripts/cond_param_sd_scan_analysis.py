#!/usr/bin/env python

"""@package docstring
File: cond_param_sd_scan_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
import re
import time
import yaml
from pathlib import Path
import h5py

# Data manipulation
import numpy as np
from scipy.special import erf
from scipy.integrate import quad

# Visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import (
    Circle, RegularPolygon, FancyArrowPatch, ArrowStyle)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, NullFormatter)
import matplotlib.colors as colors

# From alens_analysis.py
import alens_analysis as aa
from alens_analysis.colormaps import register_cmaps


def seed_scan_graphs(paramset_dir):
    # try:
    #     sd_h5_data_lst = [h5py.File(h5p, 'r+')
    #                       for h5p in paramset_dir.glob('**/analysis/*.h5')]
    #     ss_ind = sd_h5_data_lst[0]['analysis/pos_kymo'].attrs['timestep_range'][0]
    #     end_ind = sd_h5_data_lst[0]['analysis/pos_kymo'].attrs['timestep_range'][1]
    #     time_arr = sd_h5_data_lst[0]['time'][ss_ind:end_ind]
    #     start_bead = 0
    #     end_bead = None
    #     nbeads = sd_h5_data_lst[0]['raw_data']['sylinders'][start_bead:end_bead, 0, 0].shape[0]

    #     fig1, axarr1 = plt.subplots(1, 3, figsize=(24, 6))
    #     cond_num_arr, max_width_arr, total_bead_arr = aa.get_scan_cond_data(
    #         sd_h5_data_lst)

    #     aa.plot_condensate_num_sd_scan(axarr1[0], time_arr, cond_num_arr)
    #     aa.plot_condensate_size_sd_scan(
    #         axarr1[1:], time_arr, max_width_arr, total_bead_arr)
    #     fig1.tight_layout()

    #     log_avg_contact_mat = aa.get_scan_avg_contact_mat(sd_h5_data_lst)
    #     fig2, ax2 = aa.make_hic_plot(nbeads, log_avg_contact_mat, vmin=-7.)

    #     fig3, ax3 = plt.subplots(figsize=(8, 6))
    #     avg_contact_kymo = aa.get_scan_avg_kymo(sd_h5_data_lst)
    #     aa.plot_contact_kymo(fig3, ax3, time_arr, avg_contact_kymo, vmax=7.)
    # except:
    #     raise
    # finally:
    #     for h5d in sd_h5_data_lst:
    #         h5d.close()
    pass


def main():
    """TODO: Docstring for main.
    @return: TODO

    """
    # TODO make analysis dir
    # TODO
    pass


##########################################
if __name__ == "__main__":
    main()
