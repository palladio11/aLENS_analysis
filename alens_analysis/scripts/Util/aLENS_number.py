import re
import os
import glob
import argparse as agp

import numpy as np
import scipy.sparse as ss
import scipy.io as sio

import vtk
from vtk.util.numpy_support import vtk_to_numpy

import yaml
import numba as nb

import aLENS


@nb.njit
def findMove(x0: float, x1: float, L: float):
    '''x0,x1,L must be scalar FP numbers'''
    dx = np.abs(x1-x0)
    if dx > L*0.5:  # jumped across pbc boundary
        if x1 > x0:
            return x1-L-x0
        else:
            return x1+L-x0
    else:
        return x1-x0


@nb.njit(parallel=False)
def calcNematicS_numba(PList):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    # calc orientation
    # mean
    # print('Entries in list: ', len(PList))
    assert PList.shape[1] == 3
    N = PList.shape[0]
    QList = np.zeros(shape=(N, 3, 3))
    nematicOrder = np.zeros((3, 3))
    for i in range(N):
        p = PList[i]
        QList[i] = np.outer(p, p)
        nematicOrder += QList[i]

    nematicOrder *= (1.0/float(N))
    nematicOrder -= np.identity(3)/3
    # This is the correct S
    prod = 0
    for i in range(3):
        for j in range(3):
            prod += nematicOrder[i, j]*nematicOrder[i, j]
    # S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    S = np.sqrt(prod*1.5)
    return S

@nb.njit(parallel=False)
def calcPolarP_numba(PList):
    assert PList.shape[1] == 3
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    # calc orientation
    # mean
    # print('Entries in list: ', len(PList))
    N = PList.shape[0]
    polarOrder = np.zeros(3)
    for i in range(N):
        p = PList[i]
        polarOrder = polarOrder + p

    polarOrder *= (1.0/float(N))
    return polarOrder

@nb.njit
def calcCenterOrient(TList):
    '''TList must be a numpy array with shape (N,8), gid, radius, end0, end1'''
    minus_ends = TList[:, 2:5]
    plus_ends = TList[:, 5:8]
    centers = 0.5*(minus_ends+plus_ends)
    orients = plus_ends-minus_ends
    N = orients.shape[0]

    for i in nb.prange(N):
        p = orients[i]
        orients[i] = normalize(p)
    return centers, orients