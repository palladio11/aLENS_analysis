#!/usr/bin/env python

"""@package docstring
File: chrom_condensate_analysis.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import numpy as np
import warnings

from ..helpers import gen_id


class Condensate(object):
    def __init__(self, id, edge_info, merged_from=None, split_from=None):
        self.id = id
        self.time_arr = [edge_info[0]]
        self.edge_coord_arr = [edge_info[1:]]

        # First index is larger condensate
        self.merged_from = merged_from if merged_from is not None else []
        # Should only have one item
        self.split_from = split_from if split_from is not None else []

        self.merged_to = []  # Should only have one item
        self.split_to = []  # First index is larger condensate

        # self.com_arr = [self.get_edge_com(0)]
        self.edge_added = False

    def get_edge_com(self, i):
        return .5 * (self.edge_coord_arr[i].sum())

    def add_edge(self, edge_coords):
        self.time_arr += [edge_coords[0]]
        self.edge_coord_arr += [edge_coords[1:]]
        # self.com_arr += [self.get_edge_com(-1)]
        self.edge_added = True

    def __repr__(self):
        return (f'(id = {self.id}, merged_from={self.merged_from}, split_from={self.split_from}, '
                f'merged_to={self.merged_to}, split_to={self.split_to})')

    def write_analysis(self, h5_grp):
        """TODO: Docstring for write_analysis.

        @param h5_grp TODO
        @param name TODO
        @return: TODO

        """
        combined = np.hstack((np.asarray(self.time_arr)[:, np.newaxis],
                              np.asarray(self.edge_coord_arr)))
        cond_dset = h5_grp.create_dataset(f'condensate_{self.id}',
                                          data=combined)
        cond_dset.attrs['id'] = self.id
        cond_dset.attrs['merged_from'] = self.merged_from
        cond_dset.attrs['split_from'] = self.split_from
        cond_dset.attrs['merged_to'] = self.merged_to
        cond_dset.attrs['split_to'] = self.split_to

    def set_cond_from_hdf5(self, h5_data):
        """TODO: Docstring for set_cond_from_hdf5.

        @param h5_data TODO
        @return: TODO

        """
        self.time_arr = h5_data[:, 0]
        self.edge_coord_arr = h5_data[:, 1:]
        self.id = h5_data.attrs['id']
        self.merged_from = h5_data.attrs['merged_from']
        self.split_from = h5_data.attrs['split_from']
        self.merged_to = h5_data.attrs['merged_to']
        self.split_to = h5_data.attrs['split_to']


def gen_condensate_track_info(time_arr, cond_edge_arr, cond_num_arr,
                              analysis=None):
    # ts_range = h5_data['analysis']['pos_kymo'].attrs['timestep_range']
    # time_arr = h5_data['time'][ts_range[0]:ts_range[1]]
    id_gen = gen_id()

    # Contact condensate dataset
    # ind 0= time, ind 1= lower edge bead index, ind 2= higher edge bed index
    # cond_edge_dset = h5_data['analysis']['contact_cond_edges']
    # cond_num_arr = h5_data['analysis']['contact_cond_num'][...]
    n_edge_coords = cond_edge_arr[...].size

    i_ec = 0  # index of edge_coord
    stored_condensates = {}
    current_condensates = {}
    for i, t in np.ndenumerate(time_arr):
        # List of lists storing all current condensates ids
        # that have COMs within the range of prospective condensate.
        # No id is given propsective condensate so index is the identifier.
        ids_of_cur_coms_in_new = []

        # Dictionary where the keys are the ids of the current condensates
        # and values a list of all prospective condensates (by index) with COMs
        # in current condensate bead range.
        inds_of_new_coms_in_cur_dict = {}
        for cond_id in current_condensates:
            inds_of_new_coms_in_cur_dict[cond_id] = []

        # If there are no prospective condensates at this time step,
        # store all current condensates if there are any since all
        # condensates have ended
        if cond_num_arr[i] == 0:
            if current_condensates:
                stored_condensates.update(current_condensates)
                current_condensates.clear()
            continue

        # Collect connectivity between all prospective edges and current
        # condensates
        test_edges = []
        te_ind = 0
        while cond_edge_arr[i_ec, 0] < t and i_ec < n_edge_coords:
            test_edges += [cond_edge_arr[i_ec, :]]
            edge_com = .5 * (test_edges[-1][1] + test_edges[-1][2])
            ids_of_cur_coms_in_new += [[]]
            for cond_id, cond in current_condensates.items():
                com = cond.get_edge_com(-1)
                if com >= test_edges[-1][1] and com <= test_edges[-1][2]:
                    # com of condensate is in range of prospective edge
                    ids_of_cur_coms_in_new[-1] += [cond_id]
                if (edge_com >= cond.edge_coord_arr[-1][0]
                        and edge_com <= cond.edge_coord_arr[-1][1]):
                    # com of prospective edge is in range of condensate
                    inds_of_new_coms_in_cur_dict[cond_id] += [te_ind]
            te_ind += 1
            i_ec += 1

        # Continuation logic for condensates
        conds_to_add_cur = {}
        for te_i, cur_cond_ids in enumerate(ids_of_cur_coms_in_new):
            assert(len(cur_cond_ids) >= 0), "No negative events"
            # assert(len(cur_cond_ids) < 3), "No super merging events allowed"

            # No current condensates of COM in prospective edge
            #   either split or spontaneous generation event
            if len(cur_cond_ids) == 0:
                # Check to see this edge's COM is in any current condensate
                split_found = False
                for cur_id, te_inds in inds_of_new_coms_in_cur_dict.items():
                    if te_i in te_inds:  # We have a splitting event!
                        # Create new cond and add to conds_to_add_to_cur
                        new_cond = Condensate(
                            next(id_gen), test_edges[te_i], split_from=[cur_id])
                        conds_to_add_cur[new_cond.id] = new_cond
                        # Set split_to array in cur_id (always add to back
                        # split_to arr)
                        current_condensates[cur_id].split_to += [new_cond.id]
                        split_found = True
                        break

                if not split_found:  # Spontaneous condensate generation
                    # Create new condensate and add to conds_to_add_to_cur
                    new_cond = Condensate(next(id_gen), test_edges[te_i])
                    conds_to_add_cur[new_cond.id] = new_cond

            # One condensate has its COM in prospective edge
            #   either continuing a current condensate or splitting event
            elif len(cur_cond_ids) == 1:
                cur_cond = current_condensates[cur_cond_ids[0]]
                n_tes = len(inds_of_new_coms_in_cur_dict[cur_cond_ids[0]])
                if n_tes == 0:
                    new_cond = Condensate(next(id_gen), test_edges[te_i])
                    conds_to_add_cur[new_cond.id] = new_cond
                    print(f"One way edge found for cond {new_cond.id}")
                    # warnings.warn(f"One way edge found for cond {new_cond.id}")
                # assert(n_tes < 3), " Super splitter? "
                if n_tes == 1:  # The condensate continues
                    cur_cond.add_edge(test_edges[te_i])
                    continue
                elif n_tes > 1:  # Splitting event
                    if n_tes > 2:
                        print(
                            f"Super splitter event occured for cond {cur_cond.id}")
                    # Create new condensate and add to conds_to_add_to_cur
                    new_cond = Condensate(
                        next(id_gen), test_edges[te_i], split_from=[cur_cond.id])
                    conds_to_add_cur[new_cond.id] = new_cond
                    # Add new condensate id to from of split from id because
                    # com of current condensate is in new condensate range
                    # (always first by convention)
                    cur_cond.split_to = [new_cond.id] + cur_cond.split_to
            elif len(cur_cond_ids) > 1:  # Merging event
                # Create new cond and add to conds_to_add_to_cur
                new_cond = Condensate(next(id_gen), test_edges[te_i])
                conds_to_add_cur[new_cond.id] = new_cond
                if len(cur_cond_ids) > 2:
                    print(f'Super merging event occured at cond {new_cond.id}')
                for cc_id in cur_cond_ids:
                    current_condensates[cc_id].merged_to = [new_cond.id]
                    if te_i in inds_of_new_coms_in_cur_dict[cc_id]:
                        new_cond.merged_from = [cc_id] + new_cond.merged_from
                    else:
                        new_cond.merged_from += [cc_id]

        # Move all condensates that have split, merged or ended to stored
        # condensates
        conds_to_remove = []
        for i, cond in current_condensates.items():
            if cond.merged_to or cond.split_to or not cond.edge_added:
                stored_condensates[i] = cond
                conds_to_remove += [i]
            cond.edge_added = False
        for i in conds_to_remove:
            current_condensates.pop(i)

        current_condensates.update(conds_to_add_cur)
    # Add all condenstates that were still existing at end of simulation
    stored_condensates.update(current_condensates)
    if analysis is not None:
        for cond in stored_condensates.values():
            cond.write_analysis(analysis)
    return sorted(stored_condensates.values(), key=lambda cond: cond.id)


def extract_condensates(h5_grp):
    """TODO: Docstring for extract_condensates.

    @param h5_grp TODO
    @return: TODO

    """
    cond_lst = []
    for cond_dset in h5_grp.values():
        cond = Condensate(0, (0, 0, 0))  # Make dummy condensate object
        cond.set_cond_from_hdf5(cond_dset)
        cond_lst += [cond]
    return sorted(cond_lst, key=lambda cond: cond.id)


def get_max_and_total_cond_size(
        time_arr, edge_coords, cond_num_arr, analysis=None):
    """TODO: Docstring for get_max_and_total_cond_size.

    @param time_arr TODO
    @param edge_coords TODO
    @param cond_num_arr TODO
    @return: TODO

    """
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
            max_width_arr[-1] = max(max_width_arr[-1], cond_widths_arr[i_ec])
            total_width_arr[-1] += cond_widths_arr[i_ec]
            i_ec += 1
    if analysis is not None:
        max_width_dset = analysis.create_dataset(
            'max_contact_cond_size', data=max_width_arr)
        total_width_dset = analysis.create_dataset(
            'total_contact_cond_beads', data=total_width_arr)

    return np.asarray(max_width_arr), np.asarray(total_width_arr)


def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i


def get_auto_corr_fast(pos_arr):
    """Get the autocorrelation function for positions.

    @param pos_arr  TODO
    @return: TODO

    """
    nsteps = pos_arr.size
    n = next_pow_two(nsteps)

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fftn(pos_arr, s=[2 * n], axes=[-1])
    pos_corr = np.fft.ifftn(f * np.conjugate(f))[:nsteps].real

    pos_corr /= 4 * n
    return pos_corr


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
