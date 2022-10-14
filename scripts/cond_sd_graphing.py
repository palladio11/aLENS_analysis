#!/usr/bin/env python

"""@package docstring
File: cond_sd_graphing.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

# Basic useful imports
import re
import time
import yaml
import sys
from pprint import pprint
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
import alens_analysis.chromatin as aac
import alens_analysis.chromatin.chrom_analysis as ca
import alens_analysis.chromatin.chrom_condensate_analysis as cca
import alens_analysis.chromatin.chrom_graph_funcs as cgf
from alens_analysis import cluster_analysis as cla

# # Locations
# ws_path = Path('/home/alamson/DATA/Chromatin/')
# mnt_path = Path.home() / 'projects/DATA/Chromatin/'
# ceph_path = Path.home() / 'ceph/DATA/Chromatin/'

graph_sty = {
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
plt.style.use(graph_sty)


# Important analysis parameters
THRESH = .1  # Threshold for similarity between clusters to be considered
# for descendant/progenitor
NSKIP = 20  # Number of time steps that can be skipped before a cluster is
# considered to be different from a progenitor
TREEMINSIZE = 20  # Number of members a tree needs to have
MINPROGS = 3  # Number of progenitors a branch needs to be added to a tree.


def collect_data(run_path):
    ss_ind = 1
    end_ind = None
    start_bead = 0
    end_bead = None

    # Get bead position information
    with h5py.File(next(run_path.glob('analysis/raw*.h5')), 'r') as h5_data:
        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    # Get cluster information
    h5_clust_file = next(run_path.glob('analysis/cluster*.h5'))
    with h5py.File(h5_clust_file, 'r') as h5_data:
        cluster_grp = h5_data['clusters']
        time_arr = h5_data['time'][...]
        time_grp_list = sorted(cluster_grp.values(),
                               key=lambda x: x.attrs['time'])
        clusters = []
        for tg in time_grp_list:
            clusters += [[cla.Cluster(h5_data=c) for c in tg.values()]]

    return time_arr, com_arr, clusters


def make_trees(clusters):
    root_clusters = cla.find_descendants(clusters, thresh=THRESH, nskip=NSKIP)

    trees = []
    tree_id_gen = aa.helpers.gen_id()
    for root in root_clusters:
        tree = cla.ClusterTree(next(tree_id_gen))
        tree.add_recursive(root)
        if len(tree.clusters) > TREEMINSIZE:
            trees += [tree]

    # Prune smaller branches
    for tree in trees:
        tree.prune_branches(min_n_progs=MINPROGS)
    return trees


def make_seed_graphs(path):
    analysis_dir = path / 'analysis'
    time_arr, com_arr, clusters = collect_data(path)
    trees = make_trees(clusters)
    fig0, axarr0 = plt.subplots(2, 2, figsize=(20, 15))
    cgf.graph_cluster_and_tree_info_vs_time(axarr0, time_arr, trees)
    fig0.tight_layout()
    fig0.savefig(analysis_dir / 'cluster_and_tree_info_vs_time.png')

def make_all_seed_graphs(path, nthreads=5):
    # TODO Multi thread this 
    pass

def collect_all_seed_simulations(arg):
    pass

def make_seed_scan_analysis(arg):
    pass



##########################################
if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
        assert path.exists()
        make_seed_graphs(path)
    else:
        make_seed_graphs(Path.cwd())
