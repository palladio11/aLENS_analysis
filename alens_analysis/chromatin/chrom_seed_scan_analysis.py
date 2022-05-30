#!/usr/bin/env python

"""@package docstring
File: chrom_seed_scan_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import yaml
from copy import deepcopy

# Data manipulation
import numpy as np
import scipy.stats as stats
from scipy.signal import savgol_filter


def get_scan_cond_data(sd_h5_data_lst, analysis=None):
    ss_ind = sd_h5_data_lst[0]['analysis/pos_kymo'].attrs['timestep_range'][0]
    end_ind = sd_h5_data_lst[0]['analysis/pos_kymo'].attrs['timestep_range'][1]

    sd_cond_num_arr = None  # Row=time, Column=Seed
    sd_max_width_arr = None
    sd_total_bead_arr = None
    for h5d in sd_h5_data_lst:
        time_arr = h5d['time'][ss_ind:end_ind]

        cond_num_arr = h5d['analysis']['contact_cond_num'][...]
        edge_coords = h5d['analysis']['contact_cond_edges'][...]

        # Find the largest condensate by bead size at each time step
        if len(edge_coords) > 0:
            cond_widths_arr = edge_coords[:, 2] - edge_coords[:, 1]
        else:
            cond_widths_arr = np.asarray([])

        i_ec = 0  # index of edge_coord
        max_width_arr = []
        total_width_arr = []
        for i, t in np.ndenumerate(time_arr):
            max_width_arr += [0]
            total_width_arr += [0]
            # If the number of condensates at time step is zero, leave max width 0
            if cond_num_arr[i] == 0:
                continue
            # Iteratively check which is the largest condensate at a time step
            while edge_coords[i_ec, 0] < t and i_ec < cond_widths_arr.size:
                max_width_arr[-1] = max(max_width_arr[-1],
                                        cond_widths_arr[i_ec])
                total_width_arr[-1] += cond_widths_arr[i_ec]
                i_ec += 1

        if sd_cond_num_arr is None:
            sd_cond_num_arr = cond_num_arr[..., np.newaxis]
            sd_max_width_arr = np.asarray(max_width_arr)[..., np.newaxis]
            sd_total_bead_arr = np.asarray(total_width_arr)[..., np.newaxis]
        else:
            sd_cond_num_arr = np.hstack(
                (sd_cond_num_arr,
                 cond_num_arr[..., np.newaxis]))
            sd_max_width_arr = np.hstack(
                (sd_max_width_arr,
                 np.asarray(max_width_arr)[..., np.newaxis]))
            sd_total_bead_arr = np.hstack(
                (sd_total_bead_arr,
                 np.asarray(total_width_arr)[..., np.newaxis]))

    if analysis is not None:
        cond_num_dset = analysis.create_dataset(
            'cond_num', data=sd_cond_num_arr)
        cond_num_dset.attrs['timestep_range'] = [ss_ind, end_ind]

        cond_max_width_dset = analysis.create_dataset(
            'cond_max_width', data=sd_max_width_arr)
        cond_max_width_dset.attrs['timestep_range'] = [ss_ind, end_ind]

        cond_total_bead_dset = analysis.create_dataset(
            'cond_total_bead', data=sd_total_bead_arr)
        cond_total_bead_dset.attrs['timestep_range'] = [ss_ind, end_ind]
    return sd_cond_num_arr, sd_max_width_arr, sd_total_bead_arr


def get_scan_avg_contact_mat(sd_h5_data_lst, analysis=None):
    num_seeds = len(sd_h5_data_lst)
    avg_contact_mat = None
    for h5d in sd_h5_data_lst:
        is_log = h5d['analysis']['avg_contact_mat'].attrs['log']
        if avg_contact_mat is None:
            avg_contact_mat = h5d['analysis']['avg_contact_mat'][...] if not is_log else np.exp(
                h5d['analysis']['avg_contact_mat'][...])
        else:
            avg_contact_mat += h5d['analysis']['avg_contact_mat'][...] if not is_log else np.exp(
                h5d['analysis']['avg_contact_mat'][...])
    log_avg_contact_mat = np.log(avg_contact_mat) - np.log(num_seeds)
    return log_avg_contact_mat


def get_scan_avg_kymo(sd_h5_data_lst, analysis=None):
    num_seeds = len(sd_h5_data_lst)
    avg_contact_kymo = None

    for h5d in sd_h5_data_lst:
        if avg_contact_kymo is None:
            avg_contact_kymo = h5d['analysis']['contact_kymo'][...]
        else:
            avg_contact_kymo += h5d['analysis']['contact_kymo'][...]

    return avg_contact_kymo/num_seeds


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
