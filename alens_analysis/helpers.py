#!/usr/bin/env python

"""@package docstring
File: helpers.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import numpy as np


def gen_id():
    i = 0
    while True:
        yield i
        i += 1


def contiguous_regions(condition):
    """By Joe Kington, Finds contiguous True regions of the boolean array
    "condition". Returns a 2D array where the first column is the start
    index of the region and the second column is the end index.
    http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array"""

    # Find the indices of changes in "condition"
    diff = np.diff(condition.astype(int))
    idx, = diff.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll  shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # if the end of the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # if the end of condtion is True, append the last index
        idx = np.r_[idx, condition.size - 1]

    idx.shape = (-1, 2)
    return idx


def collect_contiguous_intervals(arr, delta):
    """ Collect and return different contiguous regions of an array. """
    arr_deriv = np.gradient(arr, delta)

    # Get a matrix of indices. First column is starting point of positive
    # interval, second column is ending point of positive interval
    contig_idx_arr = contiguous_regions(arr_deriv > 0)
    start_pos = contig_idx_arr[:, 0]
    end_pos = contig_idx_arr[:, 1]
    pos_lengths = (end_pos - start_pos) * 1.
    neg_lengths = (start_pos[1:] - end_pos[:-1]) * 1.

    return arr_deriv, contig_idx_arr, pos_lengths, neg_lengths


def find_steady_state_ind(arr, avg_inv=(0, None)):
    """Find the first time an array reaches an average value given an interval
    to average over (usually the end of the array).

    Parameters
    ----------
    arr : _type_
        _description_
    avg_inv : tuple, optional
        _description_, by default (0,None)
    """
    arr = np.array(arr)
    avg = arr[avg_inv[0]:avg_inv[1]].mean()
    std = arr[avg_inv[0]:avg_inv[1]].std()
    # Reaching steady-state from a smaller value
    if avg > arr[0]:
        return (arr >= (avg-std)).nonzero()[0]
    # Reaching steady-steady from a larger value
    return (arr <= (avg+std)).nonzero()[0][0]
