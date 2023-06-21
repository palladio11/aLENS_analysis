import numpy as np
import scipy as sp
import numba as nb

import os

import Util.aLENS as am
import Util.HDF5_Wrapper as h5
import point_cloud.PointCloud as pc

h5FileName = 'ptsvel'
h5.newFile(h5FileName)


def process_frame(frame):
    path = am.get_basename(frame.filename)
    minus_pts = frame.data['points'][::2, :]
    plus_pts = frame.data['points'][1::2, :]
    mpvec = (plus_pts-minus_pts)
    vel = frame.data['vel']
    omega = frame.data['omega']
    minus_vel = np.zeros(minus_pts.shape)
    plus_vel = np.zeros(plus_pts.shape)
    N = minus_vel.shape[0]
    # minus_vel[i] = vel[i]-np.cross(omega[i], mpvec[i])*0.5
    # plus_vel[i] = vel[i]+np.cross(omega[i], mpvec[i])*0.5
    rot = np.cross(omega, mpvec)*0.5
    minus_vel = vel - rot
    plus_vel = vel+rot

    h5.saveData(h5FileName, minus_pts, path, 'minus_pts', float)
    h5.saveData(h5FileName, minus_vel, path, 'minus_vel', float)
    h5.saveData(h5FileName, plus_pts, path, 'plus_pts', float)
    h5.saveData(h5FileName, plus_vel, path, 'plus_vel', float)

    return


files = am.getFileListSorted('result*-*/Sylinder_*.pvtp')

for f in files[:1000]:
    frame = am.FrameVTK(f)
    process_frame(frame)
