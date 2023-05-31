import numba as nb
import numpy as np

import Util.aLENS as am
import Util.HDF5_Wrapper as h5

# overwrite existing file
h5filename = 'OrderLinkS.hdf5'
h5.newFile(h5filename)


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
        path = am.get_basename(frame.filename)
        h5.saveData(h5filename, order, path, 'OrderLinkS', float)


if __name__ == '__main__':
    main()
