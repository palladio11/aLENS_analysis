#!/usr/bin/env python

"""@package docstring
File: create_cond_dict.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# Basic useful imports
#import alens_analysis.chromatin.chrom_condensate_analysis as cca
#import alens_analysis.chromatin.chrom_analysis as ca
from alens_analysis.cluster_analysis import Cluster
from alens_analysis.chromatin.chrom_condensate_analysis import identify_spatial_clusters
from alens_analysis.helpers import gen_id
import sys
import yaml
import numpy as np
# from pprint import pprint
from pathlib import Path
import h5py


def create_cluster_hdf5(anal_file_path, ss_ind=1, end_ind=-1, start_bead=0,
                        end_bead=None):
    id_gen = gen_id()
    with h5py.File(anal_file_path, 'r+') as h5_data:
        time_arr = h5_data['time'][ss_ind:end_ind]
        print(time_arr.shape)
        sy_dat = h5_data['raw_data/sylinders'][start_bead:end_bead,
                                               :, ss_ind:end_ind]
        com_arr = .5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    clust_path = anal_file_path.parent / 'clust_data.h5'
    if clust_path.exists():
        clust_path.unlink()
    with h5py.File(clust_path, 'w') as h5_clust:
        print(anal_file_path.parent)
        h5_clust.create_dataset('time', data=time_arr)
        clust_grp = h5_clust.create_group('clusters')
        for i, t in enumerate(time_arr):
            time_grp = clust_grp.create_group(f'time_{t}')
            time_grp.attrs['time'] = t
            clust, cluster_centers, cluster_label_inds = identify_spatial_clusters(
                com_arr[:, :, i], thresh=40)
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


##########################################
if __name__ == "__main__":
    create_cluster_hdf5(Path(sys.argv[1]))
