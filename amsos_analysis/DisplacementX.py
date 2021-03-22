import sys
import vtk
import glob
import re
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numba as nb
import scipy.optimize as so

import yaml

# TODO: This script requires updates

# read box from RunConfig.yaml
config_file = open('../RunConfig.yaml', 'r')
config = yaml.load(config_file, Loader=yaml.FullLoader)
boxLow = config['simBoxLow']
boxHigh = config['simBoxHigh']
pbcX = boxHigh[0]-boxLow[0]  # um
pbcY = boxHigh[1]-boxLow[1]  # um
pbcZ = boxHigh[2]-boxLow[2]  # um

avgWindow = [10, 20, 50, 100]  # moving average for 10 snapshots
deltat = config['timeSnap']  # time between two snapshots
pbcL = pbcX  # um
vel_max = 0.1  # um/s for histogram

config_file.close()

print('average window', avgWindow)
print('delta t', deltat)
print('pbc L', pbcL)
print('vel_max', vel_max)


@nb.njit
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


def genPosHistory():
    SylinderFileList = glob.glob('./result*/Sylinder_*.pvtp')
    # sort as numerical order
    SylinderFileList.sort(key=getFrameNumber_lambda)
    print(SylinderFileList)

    frame = Frame(SylinderFileList[0])
    posHistory = []
    for i in range(1, len(SylinderFileList)):
        # for i in range(1, 10):
        # gid, cx, cy, cz, left/right
        pos = np.zeros((len(frame.sylinders), 5))
        # get frame
        frame = Frame(SylinderFileList[i])
        for j in range(len(frame.sylinders)):
            sy = frame.sylinders[j]
            center = (np.array(sy.end0)+np.array(sy.end1))*0.5
            pos[j] = np.array(
                [sy.gid[0], center[0], center[1], center[2], np.sign(sy.end1[0]-sy.end0[0])])
        posHistory.append(pos)

    # generate file
    posHistory = np.array(posHistory)
    posHistory = np.swapaxes(posHistory, 0, 1)
    np.save('posHistory', posHistory)


def genDispHistory():
    # data: 3D array, N * steps * 5
    # shift to zero, filter PBC jump
    data = np.load('posHistory.npy')
    ntraj = data.shape[0]
    steps = data.shape[1]
    disp = np.zeros([ntraj, steps])
    for i in range(ntraj):
        for j in range(1, steps):
            x0 = data[i, j-1, 1]
            x1 = data[i, j, 1]
            dx = findMove(x0, x1, pbcL)
            disp[i, j] = disp[i, j-1]+dx
    np.save('dispHistory', disp)


def plot_fit(disp, dt, title):
    ntraj = disp.shape[0]
    nsteps = disp.shape[1]
    frame_transient = int(nsteps/2)  # find fit using the last half of data
    if ntraj < 1:
        plt.clf()
        plt.xlabel('time s')
        plt.ylabel('displacement um')
        plt.legend(title='nsamples {:6d}\n'.format(ntraj))
        plt.title(title)
        filename = title.replace(' ', '_')
        plt.savefig(filename+'.png', dpi=300)
        return
    # mean
    mean = np.mean(disp, axis=0)
    time = np.linspace(0, dt*nsteps, nsteps)
    # std
    stddev = np.std(disp, axis=0)
    # fit
    f = lambda t, *v: v[0] * t + v[1]
    popt, pcov = so.curve_fit(
        f, time[frame_transient:], mean[frame_transient:], [1, 0])
    vel = popt[0]
    offset = popt[1]
    print(popt, pcov)
    f = lambda t, *d: 2*d[0]*(t-dt*frame_transient)+d[1]
    popt, pcov = so.curve_fit(
        f, time[frame_transient:], stddev[frame_transient:]**2, [1, 0])
    print(popt, pcov)
    diff = popt[0]
    # plot
    plt.clf()
    plt.plot(time, mean, label='mean displacement')
    plt.fill_between(time, (mean-stddev), (mean+stddev),
                     alpha=0.2, label='stddev')
    plt.plot(time, time*vel+offset, '--', label='linear fit')
    plt.xlabel('time s')
    plt.ylabel('displacement um')
    plt.legend(title='nsamples {:6d}\n'.format(
        ntraj)+'V = {:6e}'.format(vel)+' um/s\n'+'D = {:6e}'.format(diff)+' um^2/s')
    plt.title(title)
    filename = title.replace(' ', '_')
    plt.savefig(filename+'.png', dpi=300)


def plotVelDiff():
    pos = np.load('posHistory.npy')
    disp = np.load('dispHistory.npy')
    # toward right, moving right
    rr = np.logical_and(disp[:, -1] > disp[:, 0], pos[:, 0, 4] > 0)
    rl = np.logical_and(disp[:, -1] > disp[:, 0], pos[:, 0, 4] < 0)
    lr = np.logical_and(disp[:, -1] < disp[:, 0], pos[:, 0, 4] > 0)
    ll = np.logical_and(disp[:, -1] < disp[:, 0], pos[:, 0, 4] < 0)
    dispRR = disp[rr]
    dispRL = disp[rl]
    dispLR = disp[lr]
    dispLL = disp[ll]
    plot_fit(dispRR, deltat, 'Orient right moving right')
    plot_fit(dispRL, deltat, 'Orient right moving left')
    plot_fit(dispLR, deltat, 'Orient left moving right')
    plot_fit(dispLL, deltat, 'Orient left moving left')


def genMovingAverage(avgsteps):
    # calc moving average
    pos = np.load('posHistory.npy')
    disp = np.load('dispHistory.npy')
    print(pos.shape)
    ntraj = pos.shape[0]
    nsteps = pos.shape[1]
    meanVelHistory = np.zeros([ntraj, nsteps])
    for i in range(ntraj):
        for j in range(nsteps-avgsteps):
            meanVelHistory[i, j] = (
                disp[i, j+avgsteps]-disp[i, j])/(deltat*avgsteps)

    print(meanVelHistory)
    np.save('meanVelHistory_'+str(avgsteps), meanVelHistory)


def plotMovingAverage(avgsteps):
    folderName = 'velx_'+str(avgsteps)
    try:
        os.mkdir(folderName)
    except FileExistsError:
        pass

    pos = np.load('posHistory.npy')
    disp = np.load('dispHistory.npy')
    meanVelHistory = np.load('meanVelHistory_'+str(avgsteps)+'.npy')
    ntraj = pos.shape[0]
    nsteps = pos.shape[1]
    # plot vx history of 5 rods
    plt.clf()
    for i in range(5):
        plt.plot(meanVelHistory[i, :], label="Gid"+str(i))
    plt.xlabel('frame index')
    plt.ylabel('Vel X, um/s')
    plt.savefig(folderName+os.sep+'VelXHistory.png', dpi=300)

    # plot vx histogram for each frame
    bins = np.linspace(-vel_max, vel_max, num=41)
    ortr = pos[:, 0, 4] > 0  # orient right
    ortl = pos[:, 0, 4] < 0  # orient left

    for i in range(nsteps):
        plt.clf()
        plt.hist(meanVelHistory[ortr][:, i], bins=bins,
                 label='orient right', density=True, alpha=0.4)
        plt.hist(meanVelHistory[ortl][:, i], bins=bins,
                 label='orient left', density=True, alpha=0.4)
        plt.legend(loc='upper left')
        plt.xlabel('Vel X, um/s')
        plt.ylabel('PDF')
        plt.ylim(0, 4/vel_max)
        plt.xlim(-vel_max, vel_max)
        plt.tick_params(axis="x", direction="in")
        plt.tick_params(axis="y", direction="in")
        plt.annotate('Snapshot {:06d}'.format(
            i), (0.7, 0.9), xycoords='axes fraction')
        plt.savefig(folderName+os.sep +
                    'velxhist_{:06d}'.format(i)+'.png', dpi=300)
    # generate movies
    path = os.path.abspath(os.path.curdir)
    names = path.split('/')
    casename = names[-2]
    cmd = "ffmpeg -y -framerate 30 -pattern_type glob -i '"+folderName+os.sep+'velxhist_' + \
        "*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 " + \
        casename+'_velxhist_'+str(avgsteps)+".mp4"
    print(cmd)
    os.system(cmd)


def main():
    if not os.path.isfile('posHistory.npy'):
        genPosHistory()
    if not os.path.isfile('dispHistory.npy'):
        genDispHistory()
    plotVelDiff()
    for steps in avgWindow:
        if not os.path.isfile('meanVelHistory_'+str(steps)+'.npy'):
            genMovingAverage(steps)
        plotMovingAverage(steps)


if __name__ == '__main__':
    main()
