import re
import os
import glob
import argparse as agp

import numpy as np
import scipy as sp
import scipy.sparse as ss
import scipy.optimize as so
import scipy.io as sio

import vtk
import yaml
import numba as nb


def get_basename(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def getFileListSorted(files, info=True):
    files = glob.glob(files)
    files.sort(key=getFrameNumber_lambda)
    if info:
        print(files)
    return files


def getDefaultArgParser(info):
    '''default argparser'''
    parser = agp.ArgumentParser(description=info)
    parser.add_argument('-c', '--config', type=str, default='../RunConfig.yaml',
                        help='path to config yaml file')
    parser.add_argument('-p', '--protein', type=str, default='../ProteinConfig.yaml',
                        help='path to protein yaml file')
    # examples
    # parser.add_argument('ngrid', type=int,
    #                     help='number of samples along X axis')
    # parser.add_argument('--rcut', type=float,
    #                     help='cut-off radius of g(r), default 0.1um', default=0.1)
    return parser


def parseConfig(yamlFile):
    config = yaml.load(open(yamlFile, 'r'), Loader=yaml.FullLoader)
    print('Config: ', config)
    return config


def getAdjacencyMatrixFromPairs(pairs, N, info=False, save=False, symmetrize=True):
    '''pairs is a list of [i,j] pairs. 0<=i,j<N'''
    if info:
        print(len(pairs))
    pairs = pairs[np.logical_and(pairs[:, 0] >= 0, pairs[:, 1] >= 0)]
    Npair = pairs.shape[0]  # number of pairs
    nbMat = ss.coo_matrix(
        (np.ones(Npair), (pairs[:, 0], pairs[:, 1])), shape=(N, N), dtype=np.int)
    if symmetrize:
        nbMat = (nbMat+nbMat.transpose())
    if save:
        sio.mmwrite('nbMat.mtx', nbMat)

    return nbMat


@nb.njit(parallel=True)
def normalize(vec):
    '''vec must be a numpy array'''
    return vec/np.sqrt(vec.dot(vec))


@nb.njit(parallel=True)
def findMove(x0, x1, L):
    '''x0,x1,L must be scalar FP numbers'''
    dx = np.abs(x1-x0)
    if dx > L*0.5:  # jumped across pbc boundary
        if x1 > x0:
            return x1-L-x0
        else:
            return x1+L-x0
    else:
        return x1-x0


def orientOrder(orientList, count=False):
    '''orientList is a list of 3D vecs'''
    # calc orientation
    # mean
    if count:
        print('Entries in list: ', len(orientList))
    PList = np.array(orientList)
    QList = np.array([np.outer(p, p) for p in PList])
    polarOrder = np.mean(PList, axis=0)
    nematicOrder = np.mean(QList, axis=0) - np.identity(3)/3
    # This is the correct S
    S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    return np.array([polarOrder[0], S])


@nb.njit(parallel=True)
def calcNematicS(PList):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
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


@nb.njit(parallel=True)
def calcPolarP(PList):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    # calc orientation
    # mean
    # print('Entries in list: ', len(PList))
    N = PList.shape[0]
    polarOrder = np.array([0, 0, 0])
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


class FrameAscii:
    '''Load Ascii.dat data'''

    def __init__(self, filename, readProtein=False, sort=True, info=False):
        self.filename = filename
        # MT
        data = np.loadtxt(filename,
                          skiprows=2, usecols=(1, 2, 3, 4, 5, 6, 7, 8))
        if sort:
            self.TList = data[data[:, 0].argsort()]  # sort by gid
        else:
            self.TList = data
        if info:
            print(self.TList[:10])
        if readProtein:
            filename = filename.replace('Sylinder', 'Protein')
            data = np.loadtxt(filename, skiprows=2,
                              usecols=(9, 10), dtype=np.int)
            if sort:
                self.PList = data[data[:, 0].argsort()]  # sort by gid
            else:
                self.PList = data
            if info:
                print(self.PList[:10])


class FrameVTK:
    '''Load VTK pvtp data. datafields are dynamically loaded.'''

    def __init__(self, sylinderFile):
        self.sylinders = []
        self.parseSylinderFile(sylinderFile)

    def parseFile(self, dataFile, objType, objList):
        print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        print("parsing data for ", nObj, " sylinders")
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

        print("-------------------------------------")
        self.sylinders.sort(key=lambda x: x.gid, reverse=False)

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)

    def printData(self):
        # output all data for debug
        for s in self.sylinders[:10]:
            # print(s.end0, s.end1)
            attrs = vars(s)
            print('*************************************')
            print('\n'.join("%s: %s" % item for item in attrs.items()))
            print('*************************************')
