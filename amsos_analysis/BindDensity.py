from numpy.lib.recfunctions import structured_to_unstructured
import numpy as np
import matplotlib.pyplot as plt
import json as js
from Util.AMSOS import ParamBase, check_inline, normalize, point_line_proj
import Util.AMSOS as am
import dask.distributed as dd
import matplotlib
matplotlib.use('Agg')


class Param(ParamBase):
    def add_argument(self, parser):
        parser.add_argument('--plot_each', type=bool, dest='plot_each',
                            default=False, help='plot hist for each frame')
        return

    def add_param(self):
        self.eps = 1e-5
        self.outputFolder = self.data_root+'/BindDensity'
        am.mkdir(self.outputFolder)
        if self.plot_each:
            self.outputPlotFolder = self.data_root+'/BindDensityPlot'
            am.mkdir(self.outputPlotFolder)
        return


def find_dbl_x(mt0, mt1, p0, p1):
    '''p0 bound on mt0, p1 bound on mt1'''
    assert check_inline(mt0[0], mt0[1], p0)
    assert check_inline(mt1[0], mt1[1], p1)

    foot0on1 = point_line_proj(p0, mt1[0], mt1[1])
    foot1on0 = point_line_proj(p1, mt0[0], mt0[1])
    norm_mt0 = normalize(mt0[1]-mt0[0])
    norm_mt1 = normalize(mt1[1]-mt1[0])

    return (p1-foot0on1).dot(norm_mt1), (p0-foot1on0).dot(norm_mt0)


def process_frame(file, param):
    frame = am.FrameAscii(file, readProtein=True, sort=False, info=False)
    name = am.get_basename(frame.filename)
    number = am.getFrameNumber_lambda(frame.filename)
    pbc = param.config['simBoxPBC']
    box = np.array(param.config['simBoxHigh']) - \
        np.array(param.config['simBoxLow'])

    TList = frame.TList
    Tid = structured_to_unstructured(TList[['gid']])
    Tm = structured_to_unstructured(TList[['mx', 'my', 'mz']])
    Tp = structured_to_unstructured(TList[['px', 'py', 'pz']])
    NT = len(TList)
    assert NT == Tm.shape[0] and NT == Tp.shape[0]
    gid_to_index = dict()
    for i, v in enumerate(Tid[:, 0]):
        if v in gid_to_index:
            print('duplicate gid detected in T list')
            exit()
        gid_to_index[v] = i

    PList = frame.PList
    Pa = structured_to_unstructured(PList[['mx', 'my', 'mz']])
    Pb = structured_to_unstructured(PList[['px', 'py', 'pz']])
    Pbind = structured_to_unstructured(PList[['idbind0', 'idbind1']])

    NP = len(PList)
    assert NP == Pa.shape[0] and NP == Pb.shape[0] and NP == Pbind.shape[0]

    uCount = 0
    saCount = 0
    sbCount = 0
    dCount = 0
    x01 = []  # distance of 1 from projection of head 0 on mt
    x10 = []  # distance of 0 from projection of head 1 on mt
    for i in range(NP):
        if Pbind[i][0] == -1 and Pbind[i][1] == -1:
            uCount += 1
        elif Pbind[i][0] != -1 and Pbind[i][1] == -1:
            saCount += 1
        elif Pbind[i][1] != -1 and Pbind[i][0] == -1:
            sbCount += 1
        else:
            dCount += 1
            # compute conformation
            mt0idx = gid_to_index[Pbind[i][0]]
            mt1idx = gid_to_index[Pbind[i][1]]
            mt0 = (Tm[mt0idx], Tp[mt0idx])
            mt1 = (Tm[mt1idx], Tp[mt1idx])
            Pcenter = 0.5*(Pa[i]+Pb[i])
            mt0 = am.find_closest_mt(mt0, Pcenter, pbc, box)
            mt1 = am.find_closest_mt(mt1, Pcenter, pbc, box)
            xab = find_dbl_x(mt0, mt1, Pa[i], Pb[i])
            x01.append(xab[0])
            x10.append(xab[1])

    tosave = {'uCount': uCount,
              'saCount': saCount,
              'sbCount': sbCount,
              'dCount': dCount,
              'head1_to_foot_of0': x01
              #   'head0_to_foot_of1': x10,
              }
    with open(param.outputFolder+'/'+name+'.json', 'w') as f:
        js.dump(tosave, f, ensure_ascii=True, indent=2)

    if param.plot_each:
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_subplot(111)
        xlim = 0.2
        if len(x01) > 0:
            ax.hist(x01, bins=100, range=[-xlim, xlim],
                    density=True, alpha=0.5)
        # if len(x10) > 0:
        #     ax.hist(x10, bins=100, range=[-xlim, xlim],
        #             density=True, alpha=0.5)
        plt.savefig(param.outputPlotFolder +
                    '/density_{:08d}'.format(number)+'.png', dpi=150)
        plt.close(fig)
        plt.cla()
        plt.clf()

    return


if __name__ == '__main__':
    param = Param('calculate motor binding density')
    files = param.syfiles

    # distribute jobs by dask
    client = dd.Client(threads_per_worker=1,
                       n_workers=param.nworkers, processes=True)
    fp = client.scatter(param, broadcast=True)

    future = client.map(process_frame,
                        [file for file in files],
                        [fp for _ in files]
                        )
    dd.wait(future)
