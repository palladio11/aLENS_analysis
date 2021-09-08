import numpy as np
from numpy.lib.recfunctions import repack_fields, structured_to_unstructured
import os
import scipy.spatial as ss
import meshzoo
import meshio
import numba as nb

import Util.AMSOS as am


parser = am.getDefaultArgParser('calc local stat on a spherical shell')

parser.add_argument('-r', '--rad', type=float,
                    default=0.25,
                    help='average radius')
parser.add_argument('-n', '--nseg', type=int,
                    default=20,
                    help='number of segments per each MT')
parser.add_argument('-m', '--mesh', type=int,
                    default=50,
                    help='order of icosa mesh')

args = parser.parse_args()

config = am.parseConfig(args.config)

R0 = config['boundaries'][0]['radius']
R1 = config['boundaries'][1]['radius']

center = np.array(config['boundaries'][0]['center'])
Rc = (R0+R1)*0.5
radAve = args.rad
nseg = args.nseg  # split each MT into nseg segments
mesh_order = args.mesh
# a cylinder with height R1-R0, approximate
volAve = np.pi*radAve*radAve*np.abs(R1-R0)
foldername = 'LocalOrder'

print(center, Rc, radAve, nseg, mesh_order, foldername, volAve)


try:
    os.mkdir(foldername)
except FileExistsError:
    pass


def calcOrient(vec):
    '''vec must be on a sphere centered at [0,0,0]'''
    rxy = np.linalg.norm(vec[:2])
    z = vec[2]
    return np.array([z*vec[0]/rxy, z*vec[1]/rxy, -rxy])


points, cells = meshzoo.icosa_sphere(mesh_order)

etheta = np.zeros((points.shape[0], 3))  # e_theta norm vectors
for i in range(points.shape[0]):
    p = points[i, :]
    etheta[i, :] = calcOrient(p)
    p = p*Rc
    points[i, :] = p+center

etheta = am.normalize_all(etheta)  # unit vectors


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
    seg_center = np.zeros((nseg*NMT, 3))
    seg_vec = np.zeros((nseg*NMT, 3))
    seg_len = np.zeros(nseg*NMT)

    for i in range(nseg):
        seg_center[i*NMT:(i+1)*NMT, :] = Tm+((i+0.5)*1.0/nseg) * Tvec
        seg_vec[i*NMT:(i+1)*NMT, :] = Tdct
        seg_len[i*NMT:(i+1)*NMT] = Tlen/nseg

    tree = ss.cKDTree(seg_center)
    search = tree.query_ball_point(pts, rad, workers=-1, return_sorted=False)
    N = pts.shape[0]
    volfrac = np.zeros(N)
    nematic = np.zeros(N)
    polarity = np.zeros((N, 3))
    polarity_theta = np.zeros(N)
    for i in range(N):
        idx = search[i]
        if len(idx) == 0:
            volfrac[i] = 0
            polarity[i, :] = np.array([0, 0, 0])
            polarity_theta[i] = 0
            nematic[i] = 0
        else:
            vecList = seg_vec[idx]
            volfrac[i] = am.volMT(0.0125, np.sum(seg_len[idx]))/volAve
            polarity[i, :] = am.calcPolarP(vecList)
            polarity_theta[i] = np.dot(polarity[i], etheta[i])
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
                                          'polarity_theta': polarity_theta,
                                          'xlinker_n_all': xlinker_n_all,
                                          'xlinker_n_db': xlinker_n_db
                                          })


SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

for file in SylinderFileList:
    frame = am.FrameAscii(file, readProtein=True, sort=False, info=True)
    calcLocalOrder(frame, points, radAve)
