from mpl_toolkits.mplot3d import proj3d
from scipy.spatial import SphericalVoronoi, geometric_slerp
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import meshzoo
import meshio

import point_cloud.PointCloud as pc

center = np.array([100.0, 100.0, 100.0])
Ri = 5.0
Ro = 5.102
Rc = (Ri+Ro)*0.5
icosa_order = 10

# create a spherical mesh on center and Rc
points, cells = meshzoo.icosa_sphere(icosa_order)
for i in range(points.shape[0]):
    p = points[i, :]
    p = p*Rc
    points[i, :] = p+center

# output mesh
meshio.write_points_cells("icosa_sphere.vtu", points,
                          cells=[("triangle", cells)])

# generate voronoi cell
sv = SphericalVoronoi(points, Rc, center)
areas = sv.calculate_areas()  # areas um to 4 pi Rc^2
vpts = sv.vertices  # voronoi points
vrgs = sv.regions  # voronoi regions

# output voronoi mesh
for rg in vrgs:
    print(rg)
    print(vpts[rg, :])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# plot generator points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')
# plot Voronoi vertices
ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')
plt.savefig('voronoi_sphere.jpg', dpi=300)

# TODO: finish this
# set input data


# def calcLocalOrder(mesh, TList):


# SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

# for file in SylinderFileList:
#     frame = am.FrameAscii(file, readProtein=False, sort=True, info=True)
#     order = calcLocalOrder(mesh, frame.TList)
