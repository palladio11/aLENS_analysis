import numba as nb
import numpy as np
import Util.aLENS as am
import Util.HDF5_Wrapper as h5
from numpy.lib.recfunctions import structured_to_unstructured

# overwrite existing file
h5filename = 'OrderLinkSP'
h5.newFile(h5filename)


def calcOrder(pairs, orients):

    N = orients.shape[0]  # number of rods
    Npair = pairs.shape[0]  # number of pairs
    print(N, Npair)

    nbMat = am.getAdjacencyMatrixFromPairs(pairs, N)

    order = np.zeros((N, 7))
    for id in nb.prange(N):
        # get the column indices of the nnz in row id of nbMat
        nnz = nbMat.getrow(id).nonzero()
        neighbors = np.append(nnz[1], id)
        PList_local = orients[neighbors]

        P = am.calcPolarP(PList_local)
        S = am.calcNematicS(PList_local)
        # print(PList_local, S, P)
        order[id, 0] = len(neighbors)  # number of rods averaged
        order[id, 1:4] = S  # nematic order S
        order[id, 4:7] = P

    return order


def calcLinkOrderSP(TList, PList):
    pairs = structured_to_unstructured(PList[['idbind0', 'idbind1']])
    centers, orients = am.calcCenterOrient(TList)
    # find neighbors for each rod
    order = calcOrder(pairs, orients)

    return order


def main():
    SylinderFileList = am.getFileListSorted(
        './result*-*/SylinderAscii_*.dat')

    for file in SylinderFileList:
        frame = am.FrameAscii(file, readProtein=True, sort=True, info=True)
        order = calcLinkOrderSP(frame.TList, frame.PList)
        print(order)
        basename = am.get_basename(frame.filename)
        path = basename.split('_')
        path[1] = path[1].zfill(8)
        print(path)
        h5.saveData(h5filename, order, "_".join(path), 'OrderLinkSP', float)


if __name__ == '__main__':
    main()
