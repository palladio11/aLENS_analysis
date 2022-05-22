import os
import numba as nb
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import h5py

import Util.aLENS as am
import Util.HDF5_Wrapper as amh5
from Util.aLENS import ParamBase


class Param(ParamBase):
    def add_argument(self, parser):
        parser.add_argument('--calcmean', type=bool,  default=False,
                            help='if calculate moving average')
        parser.add_argument('--velmax', type=float, dest='velmax', default=1.0,
                            help='vel max (um/s) for histogram')
        parser.add_argument('--navgs', type=int, nargs='+', default=[10, 20, 50, 100],
                            help='moving average for every ... snapshots')
        parser.add_argument('--axis', type=int, default=0,
                            help='0,1,2 -> x,y,z axis')
        parser.add_argument('--fittime', type=float, default=0,
                            help='fit between time = [fittime/2, fittime]')
        return

    def add_param(self):
        self.eps = 1e-5

        for navg in self.navgs:
            outputFolder = self.data_root + os.sep + \
                self.case_foldername(navg)
            am.mkdir(outputFolder)
        return

    def case_casename(self, navg):
        return 'Displacement_ax{}_{}'.format(self.axis, navg)

    def case_foldername(self, navg):
        return self.case_casename(navg)


def convert():
    '''convert hdf5 file to npy file'''
    if os.path.exists('traj_orient.npy'):
        return np.load('traj_orient.npy')

    f = h5py.File(h5file, 'r')
    keys = list(f.keys())
    nsteps = len(keys)
    npars = len(f[keys[0]]['traj'])
    traj = np.zeros((nsteps, npars, 6))
    for i in range(nsteps):
        traj[i, :, :3] = f[keys[i]]['traj']
        traj[i, :, 3:] = f[keys[i]]['vec']
    f.close()

    np.save('traj_orient.npy', traj, allow_pickle=False)

    return traj


def process(data, index, navg, params):
    assert index >= navg
    vel = (data[index, :, :3] - data[index-navg, :, :3]) / \
        (navg*params.config['timeSnap'])
    foldername = params.case_foldername(navg)
    amh5.saveData(foldername+os.sep+'Displacement',
                  vel, "{:08d}".format(index), 'vel_avg', float)

    assert vel.shape[1] == 3
    bins = np.linspace(-params.velmax, params.velmax, num=41)
    plt.clf()
    axplus = data[index, :, 3+params.axis] > 0
    axminus = data[index, :, 3+params.axis] < 0
    plt.hist(vel[axplus, params.axis], bins=bins,
             label='orient toward axis +, count = {}'.format(
                 np.count_nonzero(axplus)),
             density=True, alpha=0.4)
    plt.hist(vel[axminus, params.axis], bins=bins,
             label='orient toward axis -, count = {}'.format(
                 np.count_nonzero(axminus)),
             density=True, alpha=0.4)
    plt.legend(loc='upper left')
    plt.xlabel('Vel axis {:d}, um/s'.format(params.axis))
    plt.ylabel('PDF')
    plt.ylim(0, 4/params.velmax)
    plt.xlim(-params.velmax, params.velmax)
    plt.tick_params(axis="x", direction="in")
    plt.tick_params(axis="y", direction="in")
    plt.annotate('Snapshot {:06d}'.format(
        index), (0.7, 0.9), xycoords='axes fraction')
    plt.savefig(foldername+os.sep +
                'vel_ax{:d}_hist_{:06d}'.format(params.axis, index)+'.png', dpi=300)

    return


def plot_fit(data, params, mask, title):
    nsteps = data.shape[0]
    ntraj = data.shape[1]
    fittime = params.fittime
    # find fit using the last half of data

    disp = data[:, mask, params.axis]
    disp = disp - disp[0, ...]
    print(disp.shape)

    if disp.shape[1] < 1:
        plt.clf()
        plt.xlabel('time s')
        plt.ylabel('displacement um')
        plt.legend(title='nsamples {:6d}\n'.format(ntraj))
        plt.title(title)
        filename = title.replace(' ', '_')
        plt.savefig(filename+'.png', dpi=300)
        return

    dt = params.config['timeSnap']

    if fittime > 0:
        fitframe = int(fittime/dt)
    else:
        fitframe = ntraj

    frame_transient = fitframe//2

    # mean
    mean = np.mean(disp, axis=1)
    time = np.linspace(0, dt*nsteps, nsteps, endpoint=False)
    # std
    stddev = np.std(disp, axis=1)
    # fit
    f = lambda t, *v: v[0] * t + v[1]
    # print(time[frame_transient:fitframe], mean[frame_transient:fitframe])
    popt, pcov = so.curve_fit(
        f, time[frame_transient:fitframe], mean[frame_transient:fitframe], [1, 0])
    vel = popt[0]
    offset = popt[1]
    print(popt, pcov)
    f = lambda t, *d: 2*d[0]*(t-dt*frame_transient)+d[1]
    # print(time[frame_transient:fitframe], stddev[frame_transient:fitframe])
    popt, pcov = so.curve_fit(
        f, time[frame_transient:fitframe], stddev[frame_transient:fitframe]**2, [1, 0])
    print(popt, pcov)
    diff = popt[0]
    # plot
    plt.clf()
    plt.plot(time, mean, label='mean displacement')
    plt.fill_between(time, (mean-stddev), (mean+stddev),
                     alpha=0.2, label='stddev')
    plt.plot(time, time*vel+offset, '--', label='linear fit')
    plt.xlabel('time s')
    plt.ylabel('displacement um')
    plt.legend(title='nsamples {:6d}\n V = {:g} um/s\n D = {:g} um^2/s'.format(
        disp.shape[1], vel, diff))
    plt.title(title)
    filename = title.replace(' ', '_')
    plt.savefig(filename+'.png', dpi=300)

    return


def plotVelDiff(data, params):
    '''axis=0,1,2 means x,y,z axis'''
    # initial toward right, moving right
    disp_t0 = data[0, :, :3]
    disp_end = data[-1, :, :3]
    ort_t0 = data[0, :, 3:]

    axis = params.axis

    rr = np.logical_and(ort_t0[:, axis] > 0,
                        disp_end[:, axis] > disp_t0[:,  axis])
    rl = np.logical_and(ort_t0[:, axis] > 0,
                        disp_end[:, axis] < disp_t0[:,  axis])
    lr = np.logical_and(ort_t0[:, axis] < 0,
                        disp_end[:, axis] > disp_t0[:,  axis])
    ll = np.logical_and(ort_t0[:, axis] < 0,
                        disp_end[:, axis] < disp_t0[:,  axis])

    plot_fit(data, params, rr,
             'Orient right moving right axis {:d}'.format(axis))
    plot_fit(data, params, rl,
             'Orient right moving left axis {:d}'.format(axis))
    plot_fit(data, params, lr,
             'Orient left moving right axis {:d}'.format(axis))
    plot_fit(data, params, ll,
             'Orient left moving left axis {:d}'.format(axis))
    return


def plotTraj(traj, ntrajs, axis):

    def x(i): return [traj[_, i, axis] - traj[0, i, axis]
                      for _ in range(traj.shape[0])]
    for i in range(ntrajs):
        plt.plot(x(i), label='index {}'.format(i))

    plt.xlabel('index of snapshot')
    plt.legend()
    plt.ylabel('displacement along axis {} / um'.format(axis))
    plt.savefig('Displacement_ax{}.png'.format(axis), dpi=300)

    return


if __name__ == '__main__':
    params = Param("calculate displacement along different directions")
    h5file = 'Trajectory.hdf5'

    data = convert()

    nsteps = data.shape[0]
    npars = data.shape[1]

    plotTraj(data, 8, 0)

    plotVelDiff(data, params)

    for navg in params.navgs:
        amh5.newFile(params.case_foldername(navg)+os.sep+'Displacement')
        for i in range(params.start+navg,
                       nsteps if params.end <= 0 else params.end,
                       params.stride):
            process(data, i, navg, params)
