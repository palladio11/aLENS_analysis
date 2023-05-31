#!/usr/bin/env python

import sys

import re
import time
import yaml
from pprint import pprint
from pathlib import Path
import h5py

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


def make_motion_graph(h5_data):
    time_arr = h5_data['time'][:]
    #print(time_arr.size)
    sy_dat = h5_data['raw_data']['sylinders'][...]
    com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])
    fig, ax = plt.subplots()
    ax.plot(time_arr, np.linalg.norm(sy_dat[1,:,:],axis=0))
    return fig, ax
    

if __name__ == "__main__":
    sim_path = sys.argv[1]
    with h5py.File(next(Path(sys.argv[1]).glob('analysis/*.h5')), 'r+') as h5_data:
        fig, ax = make_motion_graph(h5_data)
        fig.savefig("something.png")