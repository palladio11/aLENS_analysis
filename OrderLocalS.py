import sys
import vtk
import glob
import re
import os

import numba as nb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py

import PointCloud as pc

pbcL = 600  # um
pbcY = 10
pbcZ = 10
rcut = 1.0  # um
box_size = np.array([pbcL, pbcY, pbcZ])

# overwrite existing file
h5filename = 'OrderLocalS.hdf5'
h5data = h5py.File(h5filename, 'w')
h5data.close()


@nb.njit
def normalize(vec):
    return vec/np.sqrt(vec.dot(vec))


@nb.njit(parallel=True)
def orientOrder(PList):
    '''PList must be a numpy array with shape (N,3)'''
    # calc orientation
    # mean
    # print('Entries in list: ', len(PList))
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


class Frame:
    def __init__(self, filename):
        data = np.loadtxt(filename,
                          skiprows=2, usecols=(1, 2, 3, 4, 5, 6, 7, 8))
        self.TList = data[data[:, 0].argsort()]
        print(self.TList[:10])
        # gid, radius, end0, end1
        self.filename = filename


def get_basename(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


@nb.njit
def getOrient(TList):
    minus_ends = TList[:, 2:5]
    plus_ends = TList[:, 5:8]
    centers = 0.5*(minus_ends+plus_ends)
    orients = plus_ends-minus_ends
    N = orients.shape[0]

    for i in nb.prange(N):
        p = orients[i]
        orients[i] = normalize(p)
    return centers, orients


@nb.njit
def calcOrder(pairs, orients):
    N = orients.shape[0]
    order = np.zeros((N, 2))
    for id in nb.prange(N):
        neighbors = pairs[pairs[:, 0] == id][:, 1]
        neighbors = np.append(neighbors, id)
        PList_local = orients[neighbors, :]
        # print(PList_local.shape)
        # if PList_local.shape == 0:
        #     print(neighbors)
        #     exit()
        S = orientOrder(PList_local)
        order[id, 0] = S  # nematic order S
        order[id, 1] = len(neighbors)  # number of rods averaged

    return order


def calcLocalOrderS(TList):
    centers, orients = getOrient(TList)
    # find neighbors for each rod
    pairs = pc.get_pair(centers, box_size, rcut)

    order = calcOrder(pairs, orients)

    return order


def main():

    SylinderFileList = glob.glob('./result*/SylinderAscii_*.dat')
    SylinderFileList.sort(key=getFrameNumber_lambda)
    print(SylinderFileList)
    for file in SylinderFileList[:1]:
        frame = Frame(file)
        order = calcLocalOrderS(frame.TList[::100])
        print(order)
        h5data = h5py.File(h5filename, 'a')
        grp = h5data.create_group(get_basename(frame.filename))
        dset = grp.create_dataset(
            "OrderLocalS", order.shape, dtype='float64', chunks=True)
        dset[...] = order[...]
        h5data.close()


if __name__ == '__main__':
    main()
