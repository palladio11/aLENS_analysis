#!/usr/bin/env python

"""@package docstring
File: create_cond_dict.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
from alens_analysis.colormaps import register_cmaps
import alens_analysis.chromatin.chrom_graph_funcs as cgf
import alens_analysis.chromatin.chrom_condensate_analysis as cca
import alens_analysis.chromatin.chrom_analysis as ca
import alens_analysis.chromatin as aac
import alens_analysis as aa
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
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import (
    Circle,
    RegularPolygon,
    FancyArrowPatch,
    ArrowStyle)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, NullFormatter)
import matplotlib.colors as colors

# Clustering stuff
from sklearn.cluster import MeanShift, estimate_bandwidth, DBSCAN, OPTICS
from itertools import cycle

colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")

# From alens_analysis.py


ss_ind = 0
end_ind = -1
start_bead = 0
end_bead = None


##########################################
if __name__ == "__main__":
    with h5py.File(next(Path.cwd().glob('*.h5')), 'r+') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        # analysis_grp = h5_data['analysis']
        print(time_arr.shape)

        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])
        clust_cent_list = []
        clust_label_list = []
        for i in range(time_arr.size):
            clust, cluster_centers, cluster_label_inds = cca.identify_spatial_clusters(
                com_arr[:, :, i], thresh=40)
            clust_cent_list += [cluster_centers]
            clust_label_list += [cluster_label_inds]
    data_dict = {
        "time_arr": time_arr.tolist(),
        "cluster_center_list": [[c.tolist() for c in t] for t in clust_cent_list],
        "cluster_label_list": [[c.tolist() for c in t] for t in clust_label_list],
    }
    with (Path.cwd() / 'clust_data.yaml').open('w') as yf:
        yaml.dump(data_dict, yf)
