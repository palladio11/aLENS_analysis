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
    sph_coord = cart2sph(xyz)
    theta = sph_coord[:, 1]
    phi = sph_coord[:, 2]
    er = np.vstack([np.sin(theta)*np.cos(phi), np.sin(theta)
                   * np.sin(phi), np.cos(theta)])
    et = np.vstack([np.cos(theta)*np.cos(phi), np.cos(theta)
                   * np.sin(phi), -np.sin(theta)])
    ep = np.vstack([-np.sin(phi), np.cos(phi), np.zeros(phi.shape[0])])
    return np.ascontiguousarray(er.T), np.ascontiguousarray(et.T), np.ascontiguousarray(ep.T)


def point_line_proj(point, p0, p1):
    '''find projection of point on p0-p1'''
    u = point-p0
    v = p1-p0
    v_norm = np.sqrt(v.dot(v))
    proj_of_u_on_v = (np.dot(u, v)/v_norm**2)*v
    proj = p0+proj_of_u_on_v  # projection of point to p line
    return proj


def find_closest_mt(mt, point, pbc, box):
    ''''''
    assert len(pbc) == 3
    assert len(box) == 3
    proj = point_line_proj(point, mt[0], mt[1])
    shift = np.zeros(3)
    for k in range(3):
        if not pbc[k]:  # ignore non-periodic direction
            continue
        candidates = [(proj[k]-box[k], -1), (proj[k], 0), (proj[k]+box[k], 1)]
        candidates.sort(key=lambda x: np.linalg.norm(x[0]-point[k]))
        shift[k] = candidates[0][1]

    return (mt[0]+shift, mt[1]+shift)


def check_inline(p0, p1, p2, eps=1e-5):
    proj = point_line_proj(p2, p0, p1)
    if np.linalg.norm(p2-proj) < eps:
        return True
    else:
        return False


class ParamBase:
    def __init__(self, text):
        parser = agp.ArgumentParser(description=text)
        parser.add_argument('--config', type=str,
                            default='../RunConfig.yaml',
                            help='path to config yaml file')
        parser.add_argument('--pconfig', type=str,
                            default='../ProteinConfig.yaml',
                            help='path to protein yaml file')
        parser.add_argument('--data_root', type=str,
                            default='.',
                            help='path to result*-* folders')
        parser.add_argument('--stride', type=int,
                            default=100,
                            help='snapshot stride')
        parser.add_argument('--start', type=int,
                            default=0,
                            help='snapshot start')
        parser.add_argument('--end', type=int,
                            default=-1,
                            help='snapshot end')
        parser.add_argument('--nworkers', type=int,
                            default=4,
                            help='number of parallel workers')

        self.add_argument(parser)

        args = parser.parse_args()
        for k, v in vars(args).items():
            setattr(self, k, v)

        self.config = parseConfig(args.config, False)
        self.protein = parseConfig(args.pconfig, False)
        self.add_param()

        print(', \n'.join("%s: %s" % item for item in vars(self).items()))

        self.syfiles = getFileListSorted(
            self.data_root+"/result*-*/SylinderAscii_*.dat", False)[self.start:self.end:self.stride]
        self.ptfiles = getFileListSorted(
            self.data_root+"/result*-*/ProteinAscii_*.dat", False)[self.start:self.end:self.stride]

        print("SylinderFiles", self.syfiles[:10])
        print("ProteinFiles", self.ptfiles[:10])

        return

    def add_argument(self, parser):
        return

    def add_param(self):
        return


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


def parseConfig(yamlFile, info=True):
    file = open(yamlFile, 'r')
    config = yaml.load(file, Loader=yaml.FullLoader)
    if info:
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


def normalize(vec):
    '''vec must be a numpy array'''
    return vec/np.linalg.norm(vec)


def normalize_all(vec):
    '''vec.shape == [N, dim]'''
    return vec/np.linalg.norm(vec, axis=1)[:, np.newaxis]


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
    w, v = np.linalg.eig(nematicOrder)
    director = normalize(v[:, np.argmax(np.abs(w))])
    return S*director


def calcPolarP(PList, weight=None):
    '''PList must be a numpy array with shape (N,3), each row normalized'''
    assert PList.shape[1] == 3
    polarOrder = np.average(PList, axis=0, weights=weight)
    return polarOrder


def calcCenterOrient(TList):
    '''TList must be a numpy array with shape (N,8), gid, radius, end0, end1'''
    assert TList.shape[1] == 8
    minus_ends = TList[:, 2:5]
    plus_ends = TList[:, 5:8]
    centers = 0.5*(minus_ends+plus_ends)
    orients = normalize_all(plus_ends-minus_ends)
    N = orients.shape[0]
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
