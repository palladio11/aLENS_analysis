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


def cart2sph(xyz):
    '''xyz.shape==(N,3), xyz => r, theta, phi'''
    assert xyz.shape[1] == 3
    xy = xyz[:, 0]**2 + xyz[:, 1]**2

    ptsnew = np.zeros(xyz.shape)
    ptsnew[:, 0] = np.sqrt(xy + xyz[:, 2]**2)  # r
    ptsnew[:, 1] = np.arctan2(np.sqrt(xy), xyz[:, 2])  # theta
    ptsnew[:, 2] = np.arctan2(xyz[:, 1], xyz[:, 0])  # phi
    return ptsnew


def e_sph(xyz):
    '''compute spherical basis vectors at vec on a spherical surface'''
    assert xyz.shape[1] == 3
    sph_coord = am.cart2sph(xyz)
    theta = sph_coord[:, 1]
    phi = sph_coord[:, 2]
    er = np.vstack([np.sin(theta)*np.cos(phi), np.sin(theta)
                   * np.sin(phi), np.cos(theta)])
    et = np.vstack([np.cos(theta)*np.cos(phi), np.cos(theta)
                   * np.sin(phi), -np.sin(theta)])
    ep = np.vstack([-np.sin(phi), np.cos(phi), np.zeros(phi.shape[0])])
    return er.T, et.T, ep.T


def volCyl(rad, h):
    '''cylinder volume'''
    return np.pi*(rad**2)*h


def volMT(rad, h):
    '''spherocylinder volume'''
    return volCyl(rad, h) + (4.0/3.0)*np.pi*(rad**3)


def mkdir(foldername):
    '''mkdir, skip if existing'''
    try:
        print('mkdir '+foldername)
        os.mkdir(foldername)
    except FileExistsError:
        print('folder already exists')
    return


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
    parser.add_argument('-c', '--config', type=str,
                        default='../RunConfig.yaml',
                        help='path to config yaml file')
    parser.add_argument('-p', '--protein', type=str,
                        default='../ProteinConfig.yaml',
                        help='path to protein yaml file')
    # examples
    # parser.add_argument('ngrid', type=int,
    #                     help='number of samples along X axis')
    # parser.add_argument('--rcut', type=float,
    #                     help='cut-off radius of g(r), default 0.1um',
    #  default=0.1)
    return parser


def parseConfig(yamlFile):
    file = open(yamlFile, 'r')
    config = yaml.load(file, Loader=yaml.FullLoader)
    print('Config: ', config)
    file.close()
    return config


def getAdjacencyMatrixFromPairs(pairs, N, info=False,
                                save=False, symmetrize=True):
    '''pairs is a list of [i,j] pairs. 0<=i,j<N'''
    if info:
        print(len(pairs))
    pairs = pairs[np.logical_and(pairs[:, 0] >= 0, pairs[:, 1] >= 0)]
    Npair = pairs.shape[0]  # number of pairs
    nbMat = ss.coo_matrix(
        (np.ones(Npair), (pairs[:, 0], pairs[:, 1])),
        shape=(N, N), dtype=np.int)
    if symmetrize:
        nbMat = (nbMat+nbMat.transpose())
    if save:
        sio.mmwrite('nbMat.mtx', nbMat)

    return nbMat


@nb.njit(parallel=True)
def normalize(vec):
    '''vec must be a numpy array'''
    return vec/np.sqrt(vec.dot(vec))


def normalize_all(vec):
    '''vec.shape == [N, dim]'''
    return vec/np.linalg.norm(vec, axis=1)[:, np.newaxis]


@nb.njit
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


def calcNematicS(PList, weight=None):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    assert PList.shape[1] == 3
    N = PList.shape[0]
    nematicOrder = np.zeros(shape=(3, 3))
    for i in range(3):
        for j in range(3):
            nematicOrder[i, j] = np.average(
                PList[:, i]*PList[:, j], axis=0, weights=weight)
    nematicOrder -= np.identity(3)/3.0
    S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    return S


@nb.njit(parallel=True)
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


def calcPolarP(PList, weight=None):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    assert PList.shape[1] == 3
    polarOrder = np.average(PList, axis=0, weights=weight)
    return polarOrder


@nb.njit(parallel=True)
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


def parseSylinderAscii(filename,  sort=True, info=False):
    fields = [('gid', np.int32), ('radius', np.float64),
              ('mx', np.float64), ('my', np.float64), ('mz', np.float64),
              ('px', np.float64), ('py', np.float64), ('pz', np.float64),
              ('group', np.int32)
              ]
    data = np.loadtxt(filename, skiprows=2,
                      usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9), dtype=fields)
    if info:
        print(data[:10])

    if sort:
        data = np.sort(data, order='gid')  # sort by gid

    return data


def parseProteinAscii(filename, sort=True, info=False):
    fields = [('gid', np.int32), ('tag', np.int32),
              ('mx', np.float64), ('my', np.float64), ('mz', np.float64),
              ('px', np.float64), ('py', np.float64), ('pz', np.float64),
              ('idbind0', np.int32), ('idbind1', np.int32)
              ]
    data = np.loadtxt(filename, skiprows=2,
                      usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), dtype=fields)
    if info:
        print(data[:10])

    if sort:
        data = np.sort(data, order='gid')  # sort by gid

    return data


class FrameAscii:
    '''Load Ascii.dat data'''

    def __init__(self, filename, readProtein=False, sort=True, info=False):
        self.filename = filename
        self.TList = parseSylinderAscii(filename, sort, info)

        if readProtein:
            filename = filename.replace('Sylinder', 'Protein')
            self.PList = parseProteinAscii(filename, sort, info)


class FrameVTK:
    '''Load VTK pvtp data. datafields are dynamically loaded.'''

    def __init__(self, dataFile):
        self.data = {}  # dict, dataname -> np.array
        self.filename = dataFile
        self.parseFile(dataFile)

    def parseFile(self, dataFile):
        print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        points = data.GetPoints()
        self.data["points"] = vtk_to_numpy(
            points.GetData())

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            print("Parsing Cell Data", dataName)
            self.data[dataName] = vtk_to_numpy(cdata)

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            print("Parsing Point Data", dataName)
            self.data[dataName] = vtk_to_numpy(pdata)

    def printData(self):
        # output all data for debug
        for attr in self.data.keys():
            print(attr, self.data[attr])
