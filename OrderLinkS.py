import sys
import vtk
import glob
import re
import os

import numba as nb
import numpy as np
import scipy as sp
import scipy.sparse as ss
import scipy.io as sio
import h5py

pbcL = 600  # um
pbcY = 10
pbcZ = 10
rcut = 1.0  # um
box_size = np.array([pbcL, pbcY, pbcZ])

# overwrite existing file
h5filename = 'OrderLinkS.hdf5'
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
        self.filename = filename
        # MT
        data = np.loadtxt(filename,
                          skiprows=2, usecols=(1, 2, 3, 4, 5, 6, 7, 8))
        self.TList = data[data[:, 0].argsort()]  # sort by gid
        print(self.TList[:10])
        # Protein
        filename = filename.replace('Sylinder', 'Protein')
        print(filename)
        data = np.loadtxt(filename, skiprows=2, usecols=(9, 10), dtype=np.int)
        self.PList = data[data[:, 0].argsort()]  # sort by gid
        print(self.PList[:10])


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


def calcOrder(pairs, orients):
    print(len(pairs))
    pairs = pairs[np.logical_and(pairs[:, 0] >= 0, pairs[:, 1] >= 0)]
    N = orients.shape[0]  # number of rods
    Npair = pairs.shape[0]  # number of pairs
    print(N, Npair)

    # adjacency matrix
    nbMat = ss.coo_matrix(
        (np.ones(Npair), (pairs[:, 0], pairs[:, 1])), shape=(N, N), dtype=np.int)
    nbMat = (nbMat+nbMat.transpose())
    # sio.mmwrite('nbMat.mtx', nbMat)

    order = np.zeros((N, 2))
    for id in nb.prange(N):
        # get the column indices of the nnz in row id of nbMat
        nnz = nbMat.getrow(id).nonzero()
        neighbors = np.append(nnz[1], id)
        PList_local = orients[neighbors]
        # print(PList_local)

        S = orientOrder(PList_local)
        order[id, 0] = S  # nematic order S
        order[id, 1] = len(neighbors)  # number of rods averaged

    return order


def calcLinkOrderS(TList, PList):
    centers, orients = getOrient(TList)
    pairs = PList[:, -2:]
    # find neighbors for each rod
    order = calcOrder(pairs, orients)

    return order


def main():

    SylinderFileList = glob.glob('./result*/SylinderAscii_*.dat')
    SylinderFileList.sort(key=getFrameNumber_lambda)
    print(SylinderFileList)
    for file in SylinderFileList[:1]:
        frame = Frame(file)
        order = calcLinkOrderS(frame.TList, frame.PList)
        print(order)
        h5data = h5py.File(h5filename, 'a')
        grp = h5data.create_group(get_basename(frame.filename))
        dset = grp.create_dataset(
            "OrderLocalS", order.shape, dtype='float64', chunks=True)
        dset[...] = order[...]
        h5data.close()


if __name__ == '__main__':
    main()
