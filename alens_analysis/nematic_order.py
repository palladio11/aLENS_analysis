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
    

