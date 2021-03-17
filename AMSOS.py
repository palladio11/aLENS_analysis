import re
import os

import numpy as np
import scipy as sp

import vtk
import numba as nb


def get_basename(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


@nb.njit(parallel=True)
def normalize(vec):
    return vec/np.sqrt(vec.dot(vec))


@nb.njit(parallel=True)
def calcNematicS(PList):
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


@nb.njit
def calcCenterOrient(TList):
    '''TList must be a numpy array with shape (N,6)'''
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

    def __init__(self, filename, readProtein=False, sort=False, info=False):
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
