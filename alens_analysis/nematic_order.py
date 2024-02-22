#!/usr/bin/env python

"""@package docstring
File: nematic_order.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


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

def make_nematic_tensor_arr(direct_arr):
    """Make the nematic tensor from a list of directors
    N x 3 array of directors

    :directors: List of directors
    :returns: Nematic tensor

    """
    lengths = np.linalg.norm(direct_arr, axis=1)
    unit_dirs = direct_arr / lengths[:, None]
    nematic_tensor_arr = (np.einsum('ij,il->ijl', unit_dirs, unit_dirs) 
                      - np.eye(3)[:,:]/3.)
    return nematic_tensor_arr

def nematic_analysis(direct_arr):
    """Takes in a list of orientations and calculates the average nematic order parameter and director

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

def fourier_transform_nematic_tensor_arr(nematic_tensor_arr, r_arr, k=[0.,0.,1.]):
    """Take the fourier transform of the nematic tensor array

    :nematic_tensor_arr: N x 3 x N array of nematic tensors
    :returns: Fourier transform of the nematic tensor array

    """
    kr = np.einsum('ni,i->n', r_arr, k)
    fQ = nematic_tensor_arr * (np.cos(kr) + 1j * np.sin(kr))[:, None, None]

    return fQ

def make_structure_factor(Q_arr, r_arr, k):
    """Make the structure factor from the fourier transform of the nematic tensor array

    :fQ: Fourier transform of the nematic tensor array
    :returns: Structure factor

    """

    fQ_arr = fourier_transform_nematic_tensor_arr(Q_arr, r_arr, k=k)

    # Find structure factor
    Snm = np.einsum('nij,mij->nm', fQ_arr, fQ_arr.conj()) 
    Snm -= np.diag((2./3.)*np.ones(Snm.shape[0]))
    return np.mean(Snm).real