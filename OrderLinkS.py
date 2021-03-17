import glob

import numba as nb
import numpy as np
import scipy as sp
import h5py

import AMSOS as am

# overwrite existing file
h5filename = 'OrderLinkS.hdf5'
h5data = h5py.File(h5filename, 'w')
h5data.close()


def calcOrder(pairs, orients):

    N = orients.shape[0]  # number of rods
    Npair = pairs.shape[0]  # number of pairs
    print(N, Npair)

    nbMat = am.getAdjacencyMatrixFromPairs(pairs, N)

    order = np.zeros((N, 2))
    for id in nb.prange(N):
        # get the column indices of the nnz in row id of nbMat
        nnz = nbMat.getrow(id).nonzero()
        neighbors = np.append(nnz[1], id)
        PList_local = orients[neighbors]

        S = am.calcNematicS(PList_local)
        order[id, 0] = S  # nematic order S
        order[id, 1] = len(neighbors)  # number of rods averaged

    return order


def calcLinkOrderS(TList, PList):
    centers, orients = am.calcCenterOrient(TList)
    pairs = PList[:, -2:]
    # find neighbors for each rod
    order = calcOrder(pairs, orients)

    return order


def main():
    SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

    for file in SylinderFileList:
        frame = am.FrameAscii(file, readProtein=True, sort=True, info=True)
        order = calcLinkOrderS(frame.TList, frame.PList)
        print(order)
        h5data = h5py.File(h5filename, 'a')
        grp = h5data.create_group(am.get_basename(frame.filename))
        dset = grp.create_dataset(
            "OrderLinkS", order.shape, dtype='float64', chunks=True)
        dset[...] = order[...]
        h5data.close()


if __name__ == '__main__':
    main()
