import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special as ss

import os
import glob
import argparse
import re

outputName = 'AsciiOrientOrder'

parser = argparse.ArgumentParser(
    description='Calculate sylinder orientational order from saved Ascii files.')
parser.add_argument('-w', type=int, dest='window', default=50,
                    help='moving average mean window')

args = parser.parse_args()
meanWindow = args.window
print('mean window: ', meanWindow)


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


class Tubule:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.minus = np.array([float(data[2]), float(data[3]), float(data[4])])
        self.plus = np.array([float(data[5]), float(data[6]), float(data[7])])
        xi = (self.plus-self.minus)
        self.orientation = xi/np.sqrt(xi.dot(xi))


class Cell:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.radius = float(data[2])
        self.minus = np.array([float(data[3]), float(data[4]), float(data[5])])
        self.plus = np.array([float(data[6]), float(data[7]), float(data[8])])
        xi = (self.plus-self.minus)
        self.orientation = xi/np.sqrt(xi.dot(xi))


class Protein:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.posBind1 = np.array(
            [float(data[2]), float(data[3]), float(data[4])])
        self.posBind2 = np.array(
            [float(data[5]), float(data[6]), float(data[7])])
        self.idBind1 = int(data[8])
        self.idBind2 = int(data[9])


class Frame:
    def __init__(self, filename):
        self.TList = []
        self.PList = []
        file = open(filename, 'r')
        for line in file:
            if line.startswith('T '):  # parse tubule data
                self.TList.append(Tubule(linestring=line))
            elif line.startswith('C '):  # parse 'cell' data, also put into TList
                self.TList.append(Cell(linestring=line))
            elif line.startswith('P '):
                self.PList.append(Protein(linestring=line))


def orientOrder(frame):
    # calc orientation
    orientList = []
    for T in frame.TList:
        orientList.append(T.orientation)
    # mean
    print('Entries in frame: ', len(frame.TList))
    PList = np.array(orientList)
    QList = np.array([np.outer(p, p) for p in PList])
    polarOrder = np.mean(PList, axis=0)
    nematicOrder = np.mean(QList, axis=0) - np.identity(3)/3
    # This is the correct S
    S = np.sqrt(np.tensordot(nematicOrder, nematicOrder)*1.5)
    return np.hstack([S, polarOrder, nematicOrder.flatten()])


files = glob.glob('./result/result*/SylinderAscii_*.dat')
# sort as numerical order
files.sort(key=getFrameNumber_lambda)
print(files)


data = []

for f in files:
    frame = Frame(f)
    # S, polarOrder, nematicOrder = orientOrder(frame)
    frame_data = orientOrder(frame)
    print(frame_data)
    S = frame_data[0]
    polarOrder = frame_data[1:4]
    nematicOrder = frame_data[4:]
    print('S: \n', S)
    print('polar order: \n', polarOrder)
    print('nematic order: \n', nematicOrder)
    data.append(frame_data)

data = np.array(data)
np.savetxt('record_'+outputName+'.csv', data, fmt='%6g', delimiter=',',
           header='S,px,py,pz,Qxx,Qxy,Qxz,Qyx,Qyy,Qyz,Qzx,Qzy,Qzz')

print("mean S: \n", np.mean(data[-meanWindow:, 0]))
print("mean polar order: \n", np.mean(data[-meanWindow:, 1:4], axis=0))
print("mean nematic order: \n", np.mean(data[-meanWindow:, 4:], axis=0))

plt.plot(np.linalg.norm(data[:, 1:4], axis=1), 'r-', label=r'$|p|$')
plt.plot(data[:, 0], 'b-',
         label=r'$S=<P_2(\cos\theta)>=\sqrt{\frac{3}{2}Q:Q}$')
plt.legend()
plt.savefig('record_'+outputName+'.png')
plt.show()
