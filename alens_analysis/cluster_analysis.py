#!/usr/bin/env python

"""@package docstring
File: cluster_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
import yaml
import sys
from copy import deepcopy
from pathlib import Path
import h5py

# Data manipulation
import numpy as np
import scipy.stats as stats
from scipy.signal import savgol_filter
from alens_analysis.helpers import gen_id
from sklearn.cluster import MeanShift, estimate_bandwidth, DBSCAN, OPTICS

# Clustering stuff
from itertools import cycle


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

        self.reset_history()

        # Possible descendants
        # self.poss_descendants = []  # List of tuples of time and cluster to compare against

        # Variables for searching history
        # self.first_progenitor = None
        # self.next_progenitor = None
        # self.last_progenitor = None
        # self.main_leaf_progenitor = None
        # self.root_descendant = None
        # self.clust_tree = None

    def reset_history(self):
        self.descendant = None
        # self.descendants = [] # List of cluster objects to compare against
        self.progenitors = []  # List of cluster objects to compare against
        self.mass_hist = 0

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

    def get_all_progenitors(self):
        """Recursively collect and return a list of all progenitors

        Returns
        -------
        list
            Cluster objects that where progentors of cluster
        """
        progs = []
        for prog in self.progenitors:
            progs += [prog]
            progs += prog.get_all_progenitors()
        return progs
        # # Alternative generator implementation
        # for prog in self.progenitors:
        #     yield prog
        #     yield from prog.get_all_progenitors()

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
        # Roots of branches that merged into the main branch
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
        if not self.main_clust_branch:
            self.update_main_clust_branch()
        return self.main_clust_branch

    def update_main_clust_branch(self):
        # First element in clusters is always the root
        self.main_clust_branch = self.clusters[0].get_largest_branch()

    def get_branch_roots(self):
        if not self.branch_roots:
            self.update_branch_roots()
        return self.branch_roots

    def update_branch_roots(self):
        self.branch_roots = [self.clusters[0]]
        for clust in self.clusters:
            # Skip the first branch for that is the main branch
            for progs in clust.progenitors[1:]:
                self.branch_roots += [progs]

    def prune_branches(self, min_n_progs=3):
        roots = self.get_branch_roots()
        del_roots = []
        for root in roots:
            progs = root.get_all_progenitors()
            if len(progs) < min_n_progs:
                del_roots += [root]
        for root in del_roots:
            for prog in root.get_all_progenitors():
                self.clusters.remove(prog)
            root.descendant.progenitors.remove(root)
            root.descendant.mass_hist -= root.mass_hist
            try:
                self.clusters.remove(root)
            except:
                "Warning could not remove cluster {root.id} from tree."
        self.update_branch_roots()


class AllClusterTrees(object):
    def __init__(self, all_clusters):
        self.trees = []
        self.all_clusters = all_clusters

    def merge_tree(self):
        pass

    def build_trees(self):
        pass


def find_descendants(clusters, thresh=.6, nskip=1):
    for t in clusters:
        for clust in t:
            clust.reset_history()
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
                # If you have found a cluster that is over the threshold,
                # make it the descendant of current node.
                # Break out of skipping loop
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


def identify_spatial_clusters(com_arr,
                              eps=0.05, min_samples=12, thresh=20, verbose=True, 
                              **kwargs):
    clust = OPTICS(min_samples=min_samples, eps=eps, cluster_method='dbscan')
    clust.fit(com_arr)
    labels = clust.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    if verbose:
        print("number of estimated clusters : %d" % n_clusters_)

    # Collect a list of cluster centers and the list of their labels
    cluster_centers = []
    cluster_label_inds = []
    for k in range(n_clusters_):
        cli = np.where(labels == k)[0]
        if cli.size < thresh:
            continue
        cluster_label_inds += [cli]
        cluster_centers += [com_arr[cli, :].mean(axis=0)]
    if verbose:
        print("number of thresholded clusters : %d" % len(cluster_centers))

    return clust, cluster_centers, cluster_label_inds


#############
## Helpers ##
#############

def collect_cluster_data(run_path,
                         ss_ind=1, end_ind=None, start_bead=0, end_bead=None,
                         **kwargs):

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
            clusters += [[Cluster(h5_data=c) for c in tg.values()]]

    return time_arr, com_arr, clusters


def get_sd_scan_cluster_num_and_bead_lst(param_dir_path):
    print(param_dir_path)
    seed_paths = [path for path in param_dir_path.glob(
        '**/s*') if path.is_dir()]
    assert seed_paths
    sd_cluster_num_lst = []
    sd_total_bead_lst = []

    for sp in seed_paths:
        time_arr, com_arr, t_grp_clusters = collect_cluster_data(sp)
        trees = make_cluster_trees(t_grp_clusters)
        num_cluster_arr = np.zeros(time_arr.shape)
        num_cluster_beads_arr = np.zeros(time_arr.shape)
        for tree in trees:
            for clust in tree.clusters:
                t_idx = np.where(time_arr == clust.time)
                num_cluster_arr[t_idx] += 1
                num_cluster_beads_arr[t_idx] += len(clust.part_ids)

        sd_cluster_num_lst += [num_cluster_arr]
        sd_total_bead_lst += [num_cluster_beads_arr]
    # TODO Pad this array with empty values
    return time_arr, sd_cluster_num_lst, sd_total_bead_lst


def make_cluster_trees(clusters,
                       thresh=.1, nskip=20, tree_min_size=20, min_progs=3,
                       **kwargs):
    root_clusters = find_descendants(clusters, thresh=thresh, nskip=nskip)

    trees = []
    tree_id_gen = gen_id()
    for root in root_clusters:
        tree = ClusterTree(next(tree_id_gen))
        tree.add_recursive(root)
        if len(tree.clusters) > tree_min_size:
            trees += [tree]

    # Prune smaller branches
    for tree in trees:
        tree.prune_branches(min_n_progs=min_progs)
    return trees


def create_cluster_hdf5(anal_file_path,
                        ss_ind=1, end_ind=None, start_bead=0, end_bead=None,
                        force=True, **kwargs):

    # Create path for cluster data file
    clust_path = (anal_file_path.parent /
                  f'cluster_{anal_file_path.parent.stem}.h5')
    if clust_path.exists():
        if not force:
            print(
                f"Warning: cluster data file {clust_path.name} exists and was not overwritten.")
            return
        clust_path.unlink()

    # Load analysis data to get particle positions for cluster algorithms
    id_gen = gen_id()
    with h5py.File(anal_file_path, 'r+') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        print(time_arr.shape)
        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    # Write cluster and write out data
    with h5py.File(clust_path, 'w') as h5_clust:
        h5_clust.create_dataset('time', data=time_arr)
        clust_grp = h5_clust.create_group('clusters')
        for i, t in enumerate(time_arr):
            time_grp = clust_grp.create_group(f'time_{t}')
            time_grp.attrs['time'] = t
            clust, cluster_centers, cluster_label_inds = identify_spatial_clusters(
                com_arr[:, :, i], **kwargs)
            for cli, cc in zip(cluster_label_inds, cluster_centers):
                cluster = Cluster(next(id_gen), t, cli, cc)
                cluster.write_clust_to_hdf5_dset(time_grp)


def create_cluster_yaml(anal_file_path, ss_ind=1, end_ind=-1, start_bead=0,
                        end_bead=None):
    with h5py.File(anal_file_path, 'r+') as h5_data:
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
    with (anal_file_path.parent / 'clust_data.yaml').open('w') as yf:
        yaml.dump(data_dict, yf)


def testing():
    ss_ind = 1
    end_ind = -1
    start_bead = 0
    end_bead = None
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
    create_cluster_hdf5(Path(sys.argv[1]))
    # testing()
