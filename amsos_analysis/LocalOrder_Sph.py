import numpy as np
import os
import scipy.spatial as ss
import meshzoo
import meshio
import numba as nb

import Util.AMSOS as am

center = np.array([100.0, 100.0, 100.0])
Ri = 5.0
Ro = 5.102
Rc = (Ri+Ro)*0.5
LMT = 0.25
radAve = LMT

mesh_order = 100
nseg = 20  # split each MT into nseg segments

foldername = 'LocalOrder'


# a cylinder with height Ro-Ri, approximate
volAve = np.pi*radAve*radAve*(Ro-Ri)
volMT = (np.pi*(0.0125**2)*LMT)+4*np.pi*(0.0125**3)/3
volSeg = volMT/nseg

try:
    os.mkdir(foldername)
except FileExistsError:
    pass

# create a spherical mesh on center and Rc
# points, cells = meshzoo.uv_sphere(
    # num_points_per_circle=mesh_order, num_circles=mesh_order)
points, cells = meshzoo.icosa_sphere(mesh_order)
for i in range(points.shape[0]):
    p = points[i, :]
    p = p*Rc
    points[i, :] = p+center


def calcLocalOrder(frame, pts, rad):
    '''pts: sample points, rad: average radius'''
    # step1: build cKDTree with TList center
    # step2: sample the vicinity of every pts
    # step3: compute average vol, P, S for every point
    TList = frame.TList
    Tm = TList[:, 2:5]
    Tp = TList[:, 5:8]
    Tvec = Tp-Tm  # vector
    Tlen = np.linalg.norm(Tvec, axis=1)  # length
    Tdct = Tvec/Tlen[:, np.newaxis]  # unit vector
    NMT = TList.shape[0]
    centers = np.zeros((nseg*NMT, 3))
    vecs = np.zeros((nseg*NMT, 3))

    for i in range(nseg):
        centers[i*NMT:(i+1)*NMT, :] = Tm+(i*1.0/nseg+0.5) * Tvec
        vecs[i*NMT:(i+1)*NMT, :] = Tdct

    tree = ss.cKDTree(centers)
    search = tree.query_ball_point(pts, rad, workers=-1, return_sorted=False)
    N = pts.shape[0]
    volfrac = np.zeros(N)
    nematic = np.zeros(N)
    polarity = np.zeros((N, 3))
    for i in range(N):
        idx = search[i]
        if len(idx) == 0:
            volfrac[i] = 0
            polarity[i, :] = np.array([0, 0, 0])
            nematic[i] = 0
        else:
            PList = vecs[idx]
            volfrac[i] = len(idx)*volSeg/volAve
            polarity[i, :] = am.calcPolarP(PList)
            nematic[i] = am.calcNematicS(PList)

    name = am.get_basename(frame.filename)
    meshio.write_points_cells(foldername+"/sphere_{}.vtu".format(name), points,
                              cells=[("triangle", cells)],
                              point_data={'volfrac': volfrac,
                                          'nematic': nematic,
                                          'polarity': polarity})


SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

for file in SylinderFileList:
    frame = am.FrameAscii(file, readProtein=False, sort=True, info=True)
    calcLocalOrder(frame, points, radAve)
