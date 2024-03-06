
# Basic useful imports
import re
import time
import yaml
from pprint import pprint
from pathlib import Path
import h5py
import warnings

# Data manipulation
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
import scipy.stats as stats
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull

# Visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import (Circle, RegularPolygon, FancyArrowPatch, ArrowStyle)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, NullFormatter)
import matplotlib.colors as mcolors

# Clustering stuff
from sklearn.cluster import MeanShift, estimate_bandwidth, DBSCAN, OPTICS
from itertools import cycle
# plt.cm.tab20.colors

# From alens_analysis.py
import alens_analysis as aa
import alens_analysis.chromatin as aac
import alens_analysis.chromatin.chrom_analysis as ca
import alens_analysis.chromatin.chrom_condensate_analysis as cca
import alens_analysis.chromatin.chrom_graph_funcs as cgf
from alens_analysis import cluster_analysis as cla

from alens_analysis.colormaps import register_cmaps

# Locations
ws_path = Path('/home/alamson/DATA/Chromatin/')
mnt_path = Path.home() / 'projects/DATA/Chromatin/'
ceph_path = Path.home() / 'ceph/DATA/Chromatin/'

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

cluster_similarity_threshold = .4 
nskip = 10 # Time snapshot skips for cluster finding. = 5 secs
vmax = 40 # Max colorbar value in kymographs
tree_length = 30 # min length of a cluster tree in time snapshots. = 15 sec

colors = cycle(mcolors.XKCD_COLORS.keys())

register_cmaps()
#plt.rcParams['image.cmap'] = 'emct8'
#plt.rcParams['image.cmap'] = 'warm'
plt.rcParams['image.cmap'] = 'YlOrRd'
#plt.rcParams['image.cmap'] = 'twilight'
#plt.rcParams['image.cmap'] = 'coolwarm'
#plt.rcParams['image.cmap'] = 'RdYlBu_r'


def cluster_analysis_graph(sim_path, part_min=40):

    flat_time_arr = []
    flat_clust_cent_arr = []
    flat_clust_ind_arr = []
    num_clusters_list = []
    num_cluster_beads_list = []

    # Length of chain and time to look at
    ss_ind = 1
    end_ind = None
    start_bead = 0
    end_bead = None
    with h5py.File(next(sim_path.glob('analysis/raw*.h5')), 'r') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        print(time_arr.shape)
        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                            :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    h5_clust_file = sim_path / 'analysis/cluster_analysis.h5'
    with h5py.File(h5_clust_file, 'r') as h5_data:
        cluster_grp = h5_data['clusters']
        time_arr = h5_data['time'][...]
        time_grp_list = sorted(cluster_grp.values(), key=lambda x: x.attrs['time'])
        clusters = []
        for tg in time_grp_list:
            clusters += [[cla.Cluster(h5_data = c) for c in tg.values()]]

    bead_ind_arr = np.zeros((com_arr.shape[0], com_arr.shape[2]))
    one_mask = np.ones(com_arr.shape[0])
    for c, clust_grp in enumerate(clusters):
        # Secondary thresholding
        clust_grp = [clust for clust in clust_grp if len(clust.part_ids) > part_min]
        num_clusters_list += [len(clust_grp)]
        num_beads = 0
        for i, clust in enumerate(clust_grp):
            flat_time_arr += [clust.time]
            flat_clust_cent_arr += [clust.center]
            flat_clust_ind_arr += [clust.part_ids]
            num_beads += len(clust.part_ids)
            bead_ind_arr[clust.part_ids,c] += one_mask[clust.part_ids]

        num_cluster_beads_list += [num_beads]
        
    # for c, clust_grp in enumerate(clusters):
    #     for clust in clust_grp:
    #         bead_ind_arr[clust.part_ids,c] += one_mask[clust.part_ids]
    # bead_ind_arr[:,:-1] += bead_ind_arr[:,1:]
    # bead_ind_arr *= .5


    fig, axarr = plt.subplots(2,2, figsize=(20, 15))
    X, Y = np.meshgrid(time_arr, np.arange(com_arr.shape[0]))
    c = axarr[1,0].pcolor(X, Y, bead_ind_arr, shading='nearest')

    flat_clust_cent_arr = np.asarray(flat_clust_cent_arr)
    _ = axarr[0,0].plot(flat_time_arr, flat_clust_cent_arr[:,0], '.' )
    _ = axarr[0,1].scatter(time_arr, num_clusters_list)
    _ = axarr[1,1].scatter(time_arr, num_cluster_beads_list)
    _ = axarr[0,0].set_ylabel('$x$-position ($\mu$m)')
    _ = axarr[0,1].set_ylabel('Number of clusters')
    _ = axarr[1,0].set_ylabel('Bead index')
    _ = axarr[1,1].set_ylabel('Number of beads in clusters')
    for ax in axarr.flatten():
        _ = ax.set_xlabel('Time (sec)')
    fig.tight_layout()


def kymo_and_cluster_graph(sim_path):
    h5_contact_file = sim_path / 'analysis/contact_analysis.h5'

    fig, axarr = plt.subplots(1,1, figsize=(10, 6))
    _ = axarr.set_ylim(0, 1600)

    with h5py.File(h5_contact_file, 'r') as h5_data:
        time_arr = h5_data['time'][...]
        contact_kymo = h5_data['contact_kymo'][...]
    
        cgf.plot_contact_kymo(fig, axarr, time_arr, contact_kymo, vmax=vmax)


    # Cluster analysis
    h5_clust_file = sim_path / 'analysis/cluster_analysis.h5'
    with h5py.File(h5_clust_file, 'r') as h5_data:
        cluster_grp = h5_data['clusters']
        time_arr = h5_data['time'][...]
        time_grp_list = sorted(cluster_grp.values(), key=lambda x: x.attrs['time'])
        clusters = []
        for tg in time_grp_list:
            clusters += [[cla.Cluster(h5_data = c) for c in tg.values()]]
    t4 = time.time()
    root_clusters = cla.find_descendants(clusters, thresh=cluster_similarity_threshold, nskip=nskip)

    trees = []
    tree_id_gen = aa.helpers.gen_id()
    for root in root_clusters:
        tree = cla.ClusterTree(next(tree_id_gen))
        tree.add_recursive(root)
        trees += [tree]

    fig, axarr = plt.subplots(1,2, figsize=(16, 6))

    # Graph all clusters
    for tree, color in zip(trees, colors):
        if len(tree.clusters) < tree_length:
            continue
        for clust in tree.clusters:
            _ = axarr[0].plot([clust.time]*len(clust.part_ids), clust.part_ids, color = color, markersize= .1, marker='.', linestyle='None')

# Graph largest branch
    biggest_tree = max(trees, key=lambda x: x.get_main_clust_branch()[0].mass_hist)
    biggest_tree.update_branch_roots()
    branch_roots = biggest_tree.get_branch_roots()

    for root, color in zip(branch_roots, colors):
        branch_clusters = root.get_largest_branch()
        for clust in branch_clusters:
            _ = axarr[1].plot([clust.time]*len(clust.part_ids), clust.part_ids, color = color, markersize= .1, marker='.', linestyle='None')


    _ = axarr[0].set_ylabel("Bead index")
    for ax in axarr:
        _ = ax.set_ylim(0, 1600)
        _ = ax.set_xlim(0, time_arr[-1])
        _ = ax.set_xlabel("Time")
