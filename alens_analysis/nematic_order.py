#!/usr/bin/env python

"""@package docstring
File: nematic_order.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import numpy as np


def calc_nematic_order(syls):
    """Calculate the nematic order parameter for a set of syls

    :syls: Raw sylinder data from HDF5 file. PBCs should be applied.
    :returns: Array of nematic order parameters for each frame

    """
    
    # Get necessary derived functions
    directions = syls[:,5:8, :] - syls[:, 2:5, :] 
    lengths = np.linalg.norm(directions, axis=1)
    unit_dirs = directions / lengths[:, None, :]
    n_syls = syls.shape[0]

    # nematic_tensor averaged over all sylinders
    nematic_tensor = (np.einsum('ijk,ilk->jlk', unit_dirs, unit_dirs) 
                      - np.eye(3)[:,:,None]/3.)/n_syls
        
    nematic_order = []
    for i in range(nematic_tensor.shape[2]):
        nematic_order.append(np.max(np.linalg.eigvals(nematic_tensor[:,:,i])))

    ## Possible vectorized version
    # ufunc_eigvals = np.vectorize(np.linalg.eigvals, signature='(n)->()')
    # nematic_order = np.max(ufunc_eigvals(nematic_tensor, axis=(0,1)), axis=0)

    return np.array(nematic_order)

def nematic_analysis(direct_arr):
    """Takes in a list of orientations and calculates the nematic order parameter and director

    Parameters
    ----------
    directors : _type_
        _description_
    """
    n = direct_arr.shape[0]
    lengths = np.linalg.norm(direct_arr, axis=1)
    unit_dirs = direct_arr / lengths[:, None, :]
    nematic_tensor = (np.einsum('ij,il->jl', unit_dirs, unit_dirs) 
                      - np.eye(3)[:,:,None]/3.)/n

    eigvals, eigvecs = np.linalg.eig(nematic_tensor)
    sort_inds = np.argsort(eigvals)
    nem_order = eigvals[sort_inds[-1]]
    nem_director = eigvecs[:, sort_inds[-1]]
    return nem_order, nem_director


