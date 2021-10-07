import os
import numba as nb
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import structured_to_unstructured

import Util.AMSOS as am
import Util.HDF5_Wrapper as h5

param = am.ParamBase('Calculate net displacement along certain direction')

config = param.config
boxsize = np.array(config['simBoxHigh'])-np.array(config['simBoxLow'])
pbc = np.array(config['simBoxPBC'])
deltat = config['timeSnap']*param.stride  # time between two snapshots


h5name = 'Trajectory'
h5.newFile(h5name)


def calcTc(TList):
    Tm = structured_to_unstructured(TList[['mx', 'my', 'mz']])
    Tp = structured_to_unstructured(TList[['px', 'py', 'pz']])
    Tc = (Tm+Tp)*0.5
    return Tc


def genTrajectory(files, start=0, end=-1):
    '''calculate center of mass trajectory and save as hdf5'''
    frame = am.FrameAscii(
        files[start], readProtein=False, sort=True, info=True)

    nMT = frame.TList.shape[0]
    Tc = calcTc(frame.TList)  # center of mass

    traj = Tc.copy()  # traj from t=0
    h5.saveData(h5name, traj, '/t_{:08d}'.format(start), 'traj', float)

    prev_Tc = None
    end = len(files)+end+1 if end < 0 else end
    for j in range(start+1, end):
        prev_Tc = Tc
        Tc = calcTc(am.FrameAscii(
            files[j], readProtein=False, sort=True, info=True).TList)
        disp = np.zeros((nMT, 3))  # disp per snapshot
        assert traj.shape == disp.shape
        for i in nb.prange(nMT):
            for k in range(3):
                x0 = prev_Tc[i, k]
                x1 = Tc[i, k]
                # filter PBC jump
                dx = am.findMove(x0, x1, boxsize[k]) if pbc[k] else x1-x0
                disp[i, k] = dx
        traj += disp
        h5.saveData(h5name, traj, '/t_{:08d}'.format(j), 'traj', float)

    return


SylinderFileList = param.syfiles
genTrajectory(SylinderFileList)
