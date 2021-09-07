import numpy as np
from numpy.lib.recfunctions import repack_fields, structured_to_unstructured
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

mesh_order = 30
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
    Tm = structured_to_unstructured(TList[['mx', 'my', 'mz']])
    Tp = structured_to_unstructured(TList[['px', 'py', 'pz']])
    Tvec = Tp-Tm  # vector
    Tlen = np.linalg.norm(Tvec, axis=1)  # length
    Tdct = Tvec/Tlen[:, np.newaxis]  # unit vector
    NMT = TList.shape[0]
    centers = np.zeros((nseg*NMT, 3))
    vecs = np.zeros((nseg*NMT, 3))

    for i in range(nseg):
        centers[i*NMT:(i+1)*NMT, :] = Tm+((i+0.5)*1.0/nseg) * Tvec
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
            vecList = vecs[idx]
            volfrac[i] = len(idx)*volSeg/volAve
            polarity[i, :] = am.calcPolarP(vecList)
            nematic[i] = am.calcNematicS(vecList)

    PList = frame.PList
    Pm = structured_to_unstructured(PList[['mx', 'my', 'mz']])
    Pp = structured_to_unstructured(PList[['px', 'py', 'pz']])
    Pbind = structured_to_unstructured(PList[['idbind0', 'idbind1']])
    xlinker_n_all = np.zeros(N)
    xlinker_n_db = np.zeros(N)
    centers = 0.5*(Pm+Pp)
    tree = ss.cKDTree(centers)
    search = tree.query_ball_point(pts, rad, workers=-1, return_sorted=False)
    for i in range(N):
        idx = search[i]
        if len(idx) == 0:
            xlinker_n_all[i] = 0
            xlinker_n_db[i] = 0
        else:
            xlinker_n_all[i] = len(idx)/volAve
            xList = Pbind[idx]
            xlinker_n_db[i] = np.count_nonzero(np.logical_and(
                xList[:, 0] != -1, xList[:, 1] != -1))/volAve

    name = am.get_basename(frame.filename)
    meshio.write_points_cells(foldername+"/sphere_{}.vtu".format(name), points,
                              cells=[("triangle", cells)],
                              point_data={'volfrac': volfrac,
                                          'nematic': nematic,
                                          'polarity': polarity,
                                          'xlinker_n_all': xlinker_n_all,
                                          'xlinker_n_db': xlinker_n_db
                                          })


SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

for file in SylinderFileList[:5]:
    frame = am.FrameAscii(file, readProtein=True, sort=False, info=True)
    calcLocalOrder(frame, points, radAve)
