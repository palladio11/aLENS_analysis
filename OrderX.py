import sys
import vtk
import glob
import re
import os

import numba as nb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

pbcL = 4  # um
gridX = np.linspace(0, pbcL, num=40)

@nb.njit
def normalize(vec):
    return vec/np.sqrt(vec.dot(vec))


def findMove(x0, x1, L):
    dx = np.abs(x1-x0)
    if dx > L*0.5:  # jumped across pbc boundary
        if x1 > x0:
            return x1-L-x0
        else:
            return x1+L-x0
    else:
        return x1-x0

# member variables are dynamically added by parsing data files
# for Sylinder, Protein, and ConBlock classes


class Sylinder(object):
    end0 = None
    end1 = None
    pass


class Frame:

    def __init__(self, sylinderFile=None):
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


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


# @nb.njit
def orientOrder(orientList):
    # calc orientation
    # mean
    print('Entries in list: ', len(orientList))
    PList = np.array(orientList)
    QList = np.array([np.outer(p, p) for p in PList])
    polarOrder = np.mean(PList, axis=0)
    nematicOrder = np.mean(QList, axis=0) - np.identity(3)/3
    # This is the correct S
    S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    return np.array([polarOrder[0], S])


def calcOrderX(sylinders):
    sypos = []
    for sy in sylinders:
        sypos.append(np.array([sy.end0, sy.end1]).flatten())

    sypos = np.array(sypos)  # N x 6 array. N: sylinder number
    order = []
    for x in gridX:
        # find those intersect with x plane
        indices0 = (sypos[:, 0] - x) * (sypos[:, 3] - x) < 0  # x
        indicesm1 = (sypos[:, 0] - (x-pbcL)) * \
            (sypos[:, 3] - (x-pbcL)) < 0  # x - pbcL
        indicesp1 = (sypos[:, 0] - (x+pbcL)) * \
            (sypos[:, 3] - (x+pbcL)) < 0  # x + pbcL
        indices = np.logical_or(
            indices0, np.logical_or(indicesm1, indicesp1))
        # print(np.count_nonzero(indices))
        # convert end0/end1 to orientation
        enddata = sypos[indices]
        orientdata = np.array([normalize(sy[3:6]-sy[0:3]) for sy in enddata])
        order.append(orientOrder(orientdata))

    return np.array(order)


def plotOrderX(dataOrderX):
    folderName = 'order_x'
    try:
        os.mkdir(folderName)
    except FileExistsError:
        pass
    nframes = dataOrderX.shape[0]
    for i in range(nframes):
        plt.clf()
        plt.plot(gridX, dataOrderX[i, :, 0], label='p_x')
        plt.plot(gridX, dataOrderX[i, :, 1], label='S')
        plt.xlim(0, pbcL)
        plt.ylim(-1, 1)
        plt.xlabel('position um')
        plt.legend(loc='lower right')
        plt.savefig(folderName+os.sep +
                    'order_x_{:06d}'.format(i)+'.png', dpi=300)
    # generate movies
    path = os.path.abspath(os.path.curdir)
    names = path.split(os.sep)
    casename = names[-2]
    cmd = "ffmpeg -y -framerate 30 -pattern_type glob -i '"+folderName+os.sep+'order_x_' + \
        "*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 " + \
        casename+"_orderx.mp4"
    print(cmd)
    os.system(cmd)


def main():

    if not os.path.isfile('order_x.npy'):
        # if True:
        SylinderFileList = glob.glob('./result*/Sylinder_*.pvtp')
        SylinderFileList.sort(key=getFrameNumber_lambda)
        print(SylinderFileList)
        data = []
        for file in SylinderFileList[:]:
            data.append(calcOrderX(Frame(file).sylinders))
        data = np.array(data)
        print(data)
        print(data.shape)
        np.save('order_x', data)

    order_x = np.load('order_x.npy')
    plotOrderX(order_x)


if __name__ == '__main__':
    main()
