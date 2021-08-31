import numpy as np
import scipy as sp
import scipy.spatial as ss
import meshzoo
import meshio


import point_cloud.PointCloud as pc

center = np.array([100.0, 100.0, 100.0])
Ri = 5.0
Ro = 5.102
Rc = (Ri+Ro)*0.5
icosa_order = 10
LMT = 0.25
radAve = 0.25

# a cylinder with height Ro-Ri, approximate
volAve = 4*np.pi*radAve*radAve*(Ro-Ri)
volMT = (np.pi*(0.0125**2)*LMT)+4*np.pi*(0.0125**3)/3

# create a spherical mesh on center and Rc
points, cells = meshzoo.icosa_sphere(icosa_order)
for i in range(points.shape[0]):
    p = points[i, :]
    p = p*Rc
    points[i, :] = p+center


def calcLocalOrder(TList, pts, rad):
    '''pts: sample points, rad: average radius'''
    # step1: build cKDTree with TList center
    # step2: sample the vicinity of every pts
    # step3: compute average vol, P, S for every point
    minus = TList[:, :3]
    plus = TList[:, 3:6]
    centers = 0.5*(minus+plus)
    tree = ss.cKDTree(centers)
    search = tree.query_ball_point(pts, rad, workers=-1)
    N = pts.shape[0]
    volfrac = np.zeros(N)
    nematic = np.zeros(N)
    polarity = np.zeros(N, 3)
    for i in range(N):
        idx = search[i]
        print(idx)
        volfrac[i] = len(idx)*volMT/volAve
        PList = TList[idx]
        PList = PList / np.linalg.norm(PList, axis=1)
        polarity[i, :] = am.calcPolarP(PList)
        nematic[i] = am.calcNematicS(PList)

    meshio.write_points_cells("icosa_sphere.vtu", points,
                              cells=[("triangle", cells)],
                              point_data=[('volfrac', volfrac),
                                          ('nematic', nematic),
                                          ('polarity'), polarity])


SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

for file in SylinderFileList:
    frame = am.FrameAscii(file, readProtein=False, sort=True, info=True)
    calcLocalOrder(frame.TList, points, radAve)
