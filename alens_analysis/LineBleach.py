import numba as nb
import numpy as np
import Util.aLENS as am
import Util.HDF5_Wrapper as h5
from Util.aLENS import ParamBase, check_inline, normalize, point_line_proj
from numpy.lib.recfunctions import structured_to_unstructured


class Param(ParamBase):
    def add_argument(self, parser):
        parser.add_argument('--num_xcuts', type=int, dest='num_xcuts',
                            default=10, help='number of x position cuts')
        return

    # def add_param(self):
    #     self.eps = 1e-5
    #     self.outputFolder = self.data_root+'/BindDensity'
    #     am.mkdir(self.outputFolder)
    #     if self.plot_each:
    #         self.outputPlotFolder = self.data_root+'/BindDensityPlot'
    #         am.mkdir(self.outputPlotFolder)
    #     return


yzlist = np.array([[0.5, 0.5], [0.5, 0.6], [0.5, 0.4], [0.4, 0.5], [0.6, 0.5]])
yzr = 0.1


def process_frame(fileidx, param, trajorient):
    # for each cut, determine polarity, left moving group velocity and right moving group velocity
    xlow = param.config['simBoxLow'][0]
    xhigh = param.config['simBoxHigh'][0]
    pbcL = xhigh-xlow
    print(fileidx)
    xgrid = np.linspace(xlow, xhigh, param.num_xcuts)
    file = param.syfiles[fileidx]
    frame = am.FrameAscii(file, readProtein=False, sort=True, info=False)
    TList = frame.TList
    N = len(TList)
    ends = np.zeros((N, 6))
    center, orient = am.calcCenterOrient(TList)
    ends[:, :3] = structured_to_unstructured(TList[['mx', 'my', 'mz']])
    ends[:, 3:] = structured_to_unstructured(TList[['px', 'py', 'pz']])

    data = np.zeros((param.num_xcuts*len(yzlist), 9))  # polarity, vL, vR
    for i, x in enumerate(xgrid):
        for jyz in range(len(yzlist)):
            yzpoint = yzlist[jyz]
            y = yzpoint[0]
            z = yzpoint[1]

            # find those intersect with x plane and close to
            indices0 = (ends[:, 0] - x) * (ends[:, 3] - x) < 0  # x
            indicesm1 = (ends[:, 0] - (x-pbcL)) * \
                (ends[:, 3] - (x-pbcL)) < 0  # x - pbcL
            indicesp1 = (ends[:, 0] - (x+pbcL)) * \
                (ends[:, 3] - (x+pbcL)) < 0  # x + pbcL
            yc = (ends[:, 1]+ends[:, 4])/2
            zc = (ends[:, 2]+ends[:, 5])/2
            indicesy = (yc-(y-yzr)) * (yc-(y+yzr)) < 0
            indicesz = (zc-(z-yzr)) * (zc-(z+yzr)) < 0

            indices = np.logical_and(
                np.logical_or(
                    indices0,
                    np.logical_or(indicesm1, indicesp1)),
                np.logical_and(indicesy, indicesz))

            if np.count_nonzero(indices) == 0:
                continue
            # print(len(orient[indices, :]))
            P = am.calcPolarP(orient[indices, :])
            vel = (trajorient[fileidx+500, indices, :3] -
                   trajorient[fileidx-500, indices, :3]) / (1000*param.config['timeSnap'])
            velR = vel[vel[:, 0] > 0]
            velL = vel[vel[:, 0] < 0]

            if velR.shape[0] == 0 or velL.shape[0] == 0:
                continue
            data[i*len(yzlist)+jyz, :3] = P
            data[i*len(yzlist)+jyz, 3:6] = np.mean(velL, axis=0)
            data[i*len(yzlist)+jyz, 6:9] = np.mean(velR, axis=0)
            # print(P, velR, velL)

    return data


if __name__ == '__main__':
    param = Param('calculate motor binding density')
    trajorient = np.load('traj_orient.npy')
    data = []
    for i in range(600, len(param.syfiles)-600):
        data.append(process_frame(i, param, trajorient))

    np.save('LineBleachX', np.array(data))

    # # distribute jobs by dask
    # files = param.syfiles
    # client = dd.Client(threads_per_worker=1,
    #                    n_workers=param.nworkers, processes=True)
    # fp = client.scatter(param, broadcast=True)

    # future = client.map(process_frame,
    #                     [file for file in files],
    #                     [fp for _ in files]
    #                     )
    # dd.wait(future)
