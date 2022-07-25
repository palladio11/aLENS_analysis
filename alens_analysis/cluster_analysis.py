#!/usr/bin/env python

"""@package docstring
File: cluster_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
import yaml
from copy import deepcopy
from pathlib import Path
import h5py

# Data manipulation
import numpy as np
import scipy.stats as stats
from scipy.signal import savgol_filter
# import alens_analysis.chromatin.chrom_condensate_analysis as cca

# Clustering stuff
from itertools import cycle

# from .helpers import contiguous_regions


class Cluster:
    """Class for a cluster of particles at a single point in time.
    Used in a cluster history tree to track a cluster through time.
    """

    def __init__(self, id, time, part_ids):
        self.id = id
        self.time = time
        self.part_ids = part_ids

        # Values given later
        self.clust_id = 0
        self.desc_id = 0
        # Variables for searching history
        self.first_progenitor = None
        self.next_progenitor = None
        self.descendant = None
        self.last_progenitor = None
        self.main_leaf_progenitor = None
        self.root_descendant = None
        self.history_tree = None
        self.visited = False

        # Possible descendants
        self.poss_descendants = []  # List of tuples of time and cluster to compare against
        self.poss_progenitors = []  # List of tuples of time and cluster to compare against

    def get_root(self):
        return (self.root_descendant if not self.root_descendant is None
                else self.descendant.get_root())

    def get_main_leaf_prog(self):
        return (self.main_leaf_progenitor if not self.main_leaf_progenitor is None
                else self.first_progenitor.get_main_leaf_prog())

    def compare(self, cluster_b):
        return len(set(self.part_ids).intersection(set(cluster_b.part_ids)))


class ClusterHistoryTree(object):
    def __init__(self, id=0):
        self.tree_id = id
        self.clusters = []

    def add_recursive(self, clust):
        clust.history_tree = self
        self.clusters += [clust]

        if clust.get_root() is None:
            clust.root_descendant = clust
        else:
            clust.root_descendant = clust.get_root()

        if clust.first_progenitor is None:
            clust.main_leaf_progenitor = clust
        else:
            self.add_recursive(clust.first_progenitor)

        clust.main_leaf_progenitor = clust.get_main_leaf_prog()
        next_prog = clust.next_progenitor
        while(not next_prog is None):
            self.add_recursive(next_prog)
            next_prog = clust.next_progenitor

        clust.last_progenitor = self.clusters[-1]


class AllClusterTrees(object):
    def __init__(self, all_clusters):
        self.trees = []
        self.all_clusters = all_clusters
        print("I made a tree")

    def merge_tree(self):
        pass

    def build_trees(self):
        pass


def testing():
    ceph_path = Path.home() / 'ceph/DATA/Chromatin/'
    # Get cluster data
    data_path = ceph_path / "22-07-15_aLc1_line16000_100umconf"


if __name__ == "__main__":
    testing()

#     ss_ind = 1
#     end_ind = -1
#     start_bead = 0
#     end_bead = None
#     data_path = (
#         ceph_path / "22_aLchr700_sticky_runs/22-01-02_aLchr1_scan.12_line700_2xsticky_3umconf_eq/simulations/s6")

#     with h5py.File(next(data_path.glob('analysis/*.h5')), 'r+') as h5_data:
#         time_arr = h5_data['time'][ss_ind:end_ind]
#         analysis_grp = h5_data['analysis']

#         sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
#                                                :, ss_ind:end_ind]
#         com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

#         clust_cent_list = []
#         clust_label_list = []
#         for i in range(time_arr.size):
#             clust, cluster_centers, cluster_label_inds = cca.identify_spatial_clusters(
#                 com_arr[:, :, i], thresh=40)
#             clust_cent_list += [cluster_centers]
#             clust_label_list += [cluster_label_inds]

#     # Construct cluster objects with possible progenitors
#     #   Store in a flatten list of tuples that contain the time
#     # Make trees
#     pass


# testing()
