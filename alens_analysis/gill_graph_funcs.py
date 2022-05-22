#!/usr/bin/env python

"""@package docstring
File: gill_graph_funcs.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import re
import time
import yaml
from pprint import pprint
from pathlib import Path
import h5py
import numpy as np


def find_avg_val_arr(time_arr_lst, val_arr_lst, n_timesteps=100):
    # Calculate an average trajectory plot
    max_t_lst = [t_arr[-1] for t_arr in time_arr_lst]
    t_max = np.max(max_t_lst)
    # Set up average arrays
    avg_time_arr = np.linspace(0, t_max, n_timesteps+1)
    avg_val_arr = np.zeros(n_timesteps+1)

    for time_arr, val_arr in zip(time_arr_lst, val_arr_lst):
        avg_val_arr += np.interp(avg_time_arr, time_arr, val_arr.flatten())

    avg_val_arr /= len(val_arr_lst)
    return avg_time_arr, avg_val_arr
