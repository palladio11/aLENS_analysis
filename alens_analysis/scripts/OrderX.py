import sys
import vtk
import glob
import re
import os

import numba as nb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Util.aLENS as am

pbcL = 600  # um
Ngrid = int(pbcL/0.1)+1
print("N Grid: ", Ngrid)
gridX = np.linspace(0, pbcL, num=Ngrid)


def calcOrderX(TList):
    center, orient = am.calcCenterOrient(TList)
    order = np.zeros(shape=(Ngrid, 4))  # px,py,pz,S
    ends = TList[:, 2:8]
    for i in nb.prange(Ngrid):
        x = gridX[i]
        # find those intersect with x plane
        indices0 = (ends[:, 0] - x) * (ends[:, 3] - x) < 0  # x
        indicesm1 = (ends[:, 0] - (x-pbcL)) * \
            (ends[:, 3] - (x-pbcL)) < 0  # x - pbcL
        indicesp1 = (ends[:, 0] - (x+pbcL)) * \
            (ends[:, 3] - (x+pbcL)) < 0  # x + pbcL
        indices = np.logical_or(indices0, np.logical_or(indicesm1, indicesp1))
        P, S = am.orientOrder(orient[indices, :])
        order[i, 0:3] = P
        order[i, 3] = S

    return order


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
        plt.plot(gridX, dataOrderX[i, :, 3], label='S')
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
        SylinderFileList = am.getFileListSorted(
            'result*-*/SylinderAscii_*.dat')

        data = []
        for file in SylinderFileList:
            data.append(calcOrderX(am.FrameAscii(file, info=True).TList))
        data = np.array(data)
        print(data)
        print(data.shape)
        np.save('order_x', data)

    order_x = np.load('order_x.npy')
    plotOrderX(order_x)


if __name__ == '__main__':
    main()
