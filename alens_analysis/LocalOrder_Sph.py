import numpy as np
from numpy.lib.recfunctions import structured_to_unstructured
import scipy.spatial as ss
import meshzoo
import meshio
import dask
import dask.distributed as dd

import codetiming as ct

import Util.AMSOS as am


class Param:
    def __init__(self):
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
        parser.add_argument('--stride', type=int,
                            default=100,
                            help='snapshot stride')
        parser.add_argument('--start', type=int,
                            default=0,
                            help='snapshot start')
        parser.add_argument('--end', type=int,
                            default=-1,
                            help='snapshot end')
        parser.add_argument('--nworkers', type=int,
                            default=4,
                            help='number of parallel dask workers')

        args = parser.parse_args()
        self.start = args.start
        self.end = args.end
        self.stride = args.stride
        self.nworkers = args.nworkers
        config = am.parseConfig(args.config)

        R0 = config['boundaries'][0]['radius']
        R1 = config['boundaries'][1]['radius']
        center = np.array(config['boundaries'][0]['center'])
        Rc = (R0+R1)*0.5
        self.radAve = args.rad
        self.volAve = np.pi*(self.radAve**2)*np.abs(R1-R0)
        self.nseg = args.nseg  # split each MT into nseg segments
        mesh_order = args.mesh
        # a cylinder with height R1-R0, approximate
        self.foldername = 'LocalOrder'
        am.mkdir(self.foldername)

        points, self.cells = meshzoo.icosa_sphere(mesh_order)
        self.er, self.etheta, self.ephi = am.e_sph(
            points)   # e_theta norm vectors
        # scale and shift
        points = points*Rc
        self.points = points + center[np.newaxis, :]

        print(', \n'.join("%s: %s" % item for item in vars(self).items()))

        return


def calcLocalOrder(file, param):
    '''pts: sample points, rad: average radius'''
    # step1: build cKDTree with TList center
    # step2: sample the vicinity of every pts
    # step3: compute average vol, P, S for every point

    timer = ct.Timer(name='CalcLocalOrder')
    timer.start()

    rad = param.radAve
    volAve = param.volAve
    nseg = param.nseg
    foldername = param.foldername
    pts = param.points
    cells = param.cells
    etheta = param.etheta

    print(file)
    frame = am.FrameAscii(file, readProtein=True, sort=False, info=False)

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
    nematic = np.zeros((N, 3))
    polarity = np.zeros((N, 3))
    for i in range(N):
        idx = search[i]
        if len(idx) != 0:
            vecList = seg_vec[idx]
            volfrac[i] = am.volMT(0.0125, np.sum(seg_len[idx]))/volAve
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
        if len(idx) != 0:
            xlinker_n_all[i] = len(idx)/volAve
            xList = Pbind[idx]
            xlinker_n_db[i] = np.count_nonzero(np.logical_and(
                xList[:, 0] != -1, xList[:, 1] != -1))/volAve

    name = am.get_basename(frame.filename)
    meshio.write_points_cells(foldername+"/sphere_{}.vtu".format(name), pts,
                              cells=[("triangle", cells)],
                              point_data={'er': param.er,
                                          'etheta': param.etheta,
                                          'ephi': param.ephi,
                                          'volfrac': volfrac,
                                          'nematic': nematic,
                                          'polarity': polarity,
                                          'xlinker_n_all': xlinker_n_all,
                                          'xlinker_n_db': xlinker_n_db
                                          })
    timer.stop()

    return


if __name__ == '__main__':
    param = Param()
    SylinderFileList = am.getFileListSorted(
        './result*-*/SylinderAscii_*.dat', info=False)
    # for file in SylinderFileList[::param.stride]:
    #     calcLocalOrder(file, param)
    client = dd.Client(threads_per_worker=1,
                       n_workers=param.nworkers, processes=True)
    fp = client.scatter(param, broadcast=True)

    files = SylinderFileList[param.start:param.end:param.stride]
    future = client.map(calcLocalOrder,
                        [file for file in files],
                        [fp for file in files]
                        )
    dd.wait(future)

    # lazy_results = []
    # for file in SylinderFileList[::param.stride]:
    #     lazy_result = dask.delayed(calcLocalOrder)(file, fp)
    #     lazy_results.append(lazy_result)
    # dask.compute(*lazy_results)
