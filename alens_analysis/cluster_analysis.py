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
from alens_analysis.chromatin.chrom_condensate_analysis import identify_spatial_clusters

# Clustering stuff
from itertools import cycle

# from .helpers import contiguous_regions


class Cluster:
    """Class for a cluster of particles at a single point in time.
    Used in a cluster history tree to track a cluster through time.
    """

    def __init__(self, id=-1, time=-1, part_ids=[], center=[], h5_data=None):
        self.id = id
        self.time = time
        self.part_ids = part_ids
        self.center = center

        if h5_data is not None:
            self.read_clust_from_hdf5_dset(h5_data)

        self.descendant = None
        self.progenitors = []  # List of cluster objects to compare against
        self.mass_hist = 0

        # Possible descendants
        # self.poss_descendants = []  # List of tuples of time and cluster to compare against

        # Variables for searching history
        # self.first_progenitor = None
        # self.next_progenitor = None
        # self.last_progenitor = None
        # self.main_leaf_progenitor = None
        # self.root_descendant = None
        # self.clust_tree = None

    def get_root(self):
        if self.root_descendant is not None:
            return self.root_descendant
        elif self.descendant is None:
            self.root_descendant = self
            return self
        else:
            self.root_descendant = self.descendant.get_root()
            return self.root_descendant

    def compare(self, cluster_b):
        return len(set(self.part_ids).intersection(set(cluster_b.part_ids)))

    def get_largest_branch(self):
        cur = self
        main_progs = [cur]
        while cur.progenitors:
            cur = cur.progenitors[0]
            main_progs += [cur]
        return main_progs

    def read_clust_from_hdf5_dset(self, h5_dset):
        self.id = h5_dset.attrs['id']
        self.time = h5_dset.attrs['time']
        self.center = h5_dset.attrs['center'][...]
        self.part_ids = h5_dset[...]

    def write_clust_to_hdf5_dset(self, h5_grp):
        dset = h5_grp.create_dataset(f'clust_{self.id}', data=self.part_ids)
        dset.attrs['id'] = self.id
        dset.attrs['time'] = self.time
        dset.attrs['center'] = self.center


class ClusterTree(object):
    def __init__(self, id=0):
        self.tree_id = id
        self.clusters = []
        self.main_clust_branch = []
        self.branch_roots = []

    def add_recursive(self, clust):
        clust.clust_tree = self
        self.clusters += [clust]

        for prog in clust.progenitors:
            assert(prog != clust)
            self.add_recursive(prog)
            clust.mass_hist += prog.mass_hist

        clust.progenitors.sort(key=lambda x: x.mass_hist, reverse=True)
        clust.mass_hist += len(clust.part_ids)

    def get_main_clust_branch(self):
        if self.main_clust_branch:
            return self.main_clust_branch
        self.update_main_clust_branch()
        return self.main_clust_branch

    def update_main_clust_branch(self):
        # First element in clusters is always the root
        self.main_clust_branch = self.clusters[0].get_largest_branch()

    def get_branch_roots(self):
        if self.branch_roots:
            return self.branch_roots
        self.update_branch_roots()
        return self.branch_roots

    def update_branch_roots(self):
        self.branch_roots = [self.clusters[0]]
        for clust in self.clusters:
            # Skip the first branch for that is the main branch
            for progs in clust.progenitors[1:]:
                self.branch_roots += [progs]


class AllClusterTrees(object):
    def __init__(self, all_clusters):
        self.trees = []
        self.all_clusters = all_clusters

    def merge_tree(self):
        pass

    def build_trees(self):
        pass


def find_descendants(clusters, thresh=.6, nskip=1):
    root_clusters = []
    for i in range(len(clusters)-1):
        for cur in clusters[i]:  # cur = current cluster
            best_score = 0
            best_cand = None
            cand_size = 0
            # Look at clusters 1 snapshot ahead and find descendant
            for s in range(1, nskip+1):
                for cand in clusters[i+s]:  # cand = descendant candidate cluster
                    score = cur.compare(cand)
                    if score > best_score:
                        best_score = score
                        best_cand = cand
                        cand_size = len(cand.part_ids)
                # If you have found a cluster that is over the threshold, make it the descendant of current node. Break out of skipping loop
                if best_score > thresh * min(len(cur.part_ids), cand_size):
                    cur.descendant = best_cand
                    best_cand.progenitors += [cur]
                    assert(best_cand != cur)
                    break

                # Don't look for a snapshot that is beyond the time index i+s
                if i+s == len(clusters)-1:
                    break
            if cur.descendant is None:
                # assert(cur.descendant is None)
                root_clusters += [cur]
    # For clusters existing at last time point, add them to root clusters
    for clust in clusters[-1]:
        assert(clust.descendant is None)
        root_clusters += [clust]

    return root_clusters


#############
## Helpers ##
#############


def create_cluster_dict(analysis_path, start_bead=0, end_bead=None, ss_ind=1, end_ind=None):
    with h5py.File(next(analysis_path.glob('*.h5')), 'r+') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        print(time_arr.shape)

        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])
        clust_cent_list = []
        clust_label_list = []
        for i in range(time_arr.size):
            clust, cluster_centers, cluster_label_inds = identify_spatial_clusters(
                com_arr[:, :, i], thresh=40)
            clust_cent_list += [cluster_centers]
            clust_label_list += [cluster_label_inds]
    data_dict = {
        "time_arr": time_arr.tolist(),
        "cluster_center_list": [[c.tolist() for c in t] for t in clust_cent_list],
        "cluster_label_list": [[c.tolist() for c in t] for t in clust_label_list],
    }
    with (analysis_path / 'clust_data.yaml').open('w') as yf:
        yaml.dump(data_dict, yf)


def testing():
    ss_ind = 1
    end_ind = -1
    start_bead = 0
    end_bead = None
    # data_path = Path(
    #     '/home/alamson/DATA/Chromatin/22-04-28_aLchr1_scan.12_line800_sticky55nm_eps1_Ka30_5umconf/simulations/s10')
    data_path = (Path.home() /
                 'ceph/DATA/Chromatin/22_aLchr700_sticky_runs/22-01-02_aLchr1_scan.12_line700_2xsticky_3umconf_eq/simulations/s6')
    # Get cluster data
    with h5py.File(next(data_path.glob('analysis/*.h5')), 'r+') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        analysis_grp = h5_data['analysis']
        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])


if __name__ == "__main__":
    testing()

    # testing()

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
