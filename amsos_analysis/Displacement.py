import os
import numba as nb
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt


import Util.AMSOS as am


parser = am.getDefaultArgParser(
    'Calculate net displacement along certain direction')
parser.add_argument('--calcmean', type=bool, dest='calcmean', default=False,
                    help='if calculate moving average')
parser.add_argument('--velmax', type=float, dest='velmax', default=1.0,
                    help='vel max (um/s) for histogram')
parser.add_argument('--avg', type=int, nargs='+', dest='avg', default=[10, 20, 50, 100],
                    help='moving average for every ... snapshots')
parser.add_argument('--axis', type=int, dest='axis', default=0,
                    help='0,1,2 -> x,y,z axis')
args = parser.parse_args()

config = am.parseConfig(args.config)
boxsize = np.array(config['simBoxHigh'])-np.array(config['simBoxLow'])
pbc = np.array(config['simBoxPBC'])
deltat = config['timeSnap']  # time between two snapshots


avgWindow = args.avg  # moving average for every ... snapshots
vel_max = args.velmax  # um/s for histogram

print('config: ', args.config)
print('box_size:', boxsize)
print('average window', avgWindow)
print('delta t', deltat)
print('vel_max', vel_max)
print('plot axis', args.axis)


def genPosHistory():
    '''posHistory, nMT*nSteps*7: [gid, center, orientation]'''
    SylinderFileList = am.getFileListSorted('./result*/SylinderAscii_*.dat')

    frame = am.FrameAscii(SylinderFileList[0])
    nMT = frame.TList.shape[0]
    nSteps = len(SylinderFileList)
    posHistory = np.zeros(shape=[nMT, nSteps, 7])

    for i in nb.prange(0, len(SylinderFileList)):
        mt = am.FrameAscii(SylinderFileList[i]).TList
        posHistory[:, i, 0] = mt[:, 0]  # gid
        posHistory[:, i, 1:4] = 0.5*(mt[:, 2:5]+mt[:, 5:8])  # center
        posHistory[:, i, 4:7] = mt[:, 5:8]-mt[:, 2:5]  # orientation

    return posHistory


@nb.jit
def genDispHistory(pos):
    '''dispHistory, nMT*nSteps*3: [displacement]'''
    nMT = pos.shape[0]
    steps = pos.shape[1]
    disp = np.zeros((nMT, steps, 3))
    for i in nb.prange(nMT):
        for j in range(1, steps):
            for k in range(3):
                x0 = pos[i, j-1, 1+k]
                x1 = pos[i, j, 1+k]
                # filter PBC jump
                dx = am.findMove(x0, x1, boxsize[k]) if pbc[k] else x1-x0
                disp[i, j, k] = disp[i, j-1, k]+dx

    return disp


@nb.jit
def genMovingAverage(disp, avgsteps):
    '''meanVelHistory, nMT*(nSteps-avgsteps)*3, [meanVel] '''
    nMT = disp.shape[0]
    nsteps = disp.shape[1]
    meanVel = np.zeros((nMT, nsteps-avgsteps, 3))
    avgTime = deltat * avgsteps
    for i in nb.prange(nMT):
        for j in range(nsteps-avgsteps):
            meanVel[i, j] = (disp[i, j+avgsteps, :]-disp[i, j, :])/avgTime

    return meanVel


def plot_fit(disp, dt, title):
    ntraj = disp.shape[0]
    nsteps = disp.shape[1]
    frame_transient = int(nsteps/2)  # find fit using the last half of data
    if ntraj < 1:
        plt.clf()
        plt.xlabel('time s')
        plt.ylabel('displacement um')
        plt.legend(title='nsamples {:6d}\n'.format(ntraj))
        plt.title(title)
        filename = title.replace(' ', '_')
        plt.savefig(filename+'.png', dpi=300)
        return
    # mean
    mean = np.mean(disp, axis=0)
    time = np.linspace(0, dt*nsteps, nsteps)
    # std
    stddev = np.std(disp, axis=0)
    # fit
    f = lambda t, *v: v[0] * t + v[1]
    popt, pcov = so.curve_fit(
        f, time[frame_transient:], mean[frame_transient:], [1, 0])
    vel = popt[0]
    offset = popt[1]
    print(popt, pcov)
    f = lambda t, *d: 2*d[0]*(t-dt*frame_transient)+d[1]
    popt, pcov = so.curve_fit(
        f, time[frame_transient:], stddev[frame_transient:]**2, [1, 0])
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
    plt.legend(title='nsamples {:6d}\n'.format(
        ntraj)+'V = {:6e}'.format(vel)+' um/s\n'+'D = {:6e}'.format(diff)+' um^2/s')
    plt.title(title)
    filename = title.replace(' ', '_')
    plt.savefig(filename+'.png', dpi=300)


def plotVelDiff(pos, disp, axis=0):
    '''axis=0,1,2 means x,y,z axis'''
    # initial toward right, moving right
    rr = np.logical_and(pos[:, 0, 4+axis] > 0,
                        disp[:, -1, axis] > disp[:, 0, axis])
    rl = np.logical_and(pos[:, 0, 4+axis] > 0,
                        disp[:, -1, axis] < disp[:, 0, axis])
    lr = np.logical_and(pos[:, 0, 4+axis] < 0,
                        disp[:, -1, axis] > disp[:, 0, axis])
    ll = np.logical_and(pos[:, 0, 4+axis] < 0,
                        disp[:, -1, axis] < disp[:, 0, axis])

    dispRR = disp[rr, :, axis]
    dispRL = disp[rl, :, axis]
    dispLR = disp[lr, :, axis]
    dispLL = disp[ll, :, axis]
    plot_fit(dispRR, deltat,
             'Orient right moving right axis {:d}'.format(axis))
    plot_fit(dispRL, deltat,
             'Orient right moving left axis {:d}'.format(axis))
    plot_fit(dispLR, deltat,
             'Orient left moving right axis {:d}'.format(axis))
    plot_fit(dispLL, deltat,
             'Orient left moving left axis {:d}'.format(axis))


def plotMovingAverage(avgsteps, axis=0):
    folderName = 'velx_'+str(avgsteps)
    try:
        os.mkdir(folderName)
    except FileExistsError:
        pass

    pos = np.load('posHistory.npy')
    meanVelHistory = np.load('meanVelHistory_'+str(avgsteps)+'.npy')
    # plot vx history of 5 rods
    plt.clf()
    for i in range(5):
        plt.plot(meanVelHistory[i, :, axis], label="Gid"+str(i))
    plt.xlabel('frame index')
    plt.ylabel('Vel {:d}, um/s'.format(axis))
    plt.savefig(folderName+os.sep +
                'Vel_ax{:d}_History.png'.format(axis), dpi=300)

    # plot vx histogram for each frame
    bins = np.linspace(-vel_max, vel_max, num=41)
    ortr = pos[:, 0, 4+axis] > 0  # orient right
    ortl = pos[:, 0, 4+axis] < 0  # orient left

    nsteps = meanVelHistory.shape[1]
    for i in range(nsteps):
        plt.clf()
        plt.hist(meanVelHistory[ortr][:, i, axis], bins=bins,
                 label='orient right', density=True, alpha=0.4)
        plt.hist(meanVelHistory[ortl][:, i, axis], bins=bins,
                 label='orient left', density=True, alpha=0.4)
        plt.legend(loc='upper left')
        plt.xlabel('Vel axis {:d}, um/s'.format(axis))
        plt.ylabel('PDF')
        plt.ylim(0, 4/vel_max)
        plt.xlim(-vel_max, vel_max)
        plt.tick_params(axis="x", direction="in")
        plt.tick_params(axis="y", direction="in")
        plt.annotate('Snapshot {:06d}'.format(
            i), (0.7, 0.9), xycoords='axes fraction')
        plt.savefig(folderName+os.sep +
                    'vel_ax{:d}_hist_{:06d}'.format(axis, i)+'.png', dpi=300)
    # generate movies
    path = os.path.abspath(os.path.curdir)
    names = path.split('/')
    casename = names[-2]
    cmd = "ffmpeg -y -framerate 30 -pattern_type glob -i '"+folderName+os.sep+'vel_ax{:d}_hist_'.format(axis) + \
        "*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 " + \
        casename+'_vel_ax{:d}_hist_'.format(axis)+str(avgsteps)+".mp4"
    print(cmd)
    os.system(cmd)


def main():
    if not os.path.isfile('posHistory.npy'):
        pos = genPosHistory()
        np.save('posHistory', pos)
    else:
        pos = np.load('posHistory.npy')

    if not os.path.isfile('dispHistory.npy'):
        disp = genDispHistory(pos)
        np.save('dispHistory', disp)
    else:
        disp = np.load('dispHistory.npy')

    plotVelDiff(pos, disp, args.axis)

    if not args.calcmean:
        exit()

    for steps in avgWindow:
        avgfname = 'meanVelHistory_'+str(steps)
        if not os.path.isfile(avgfname+'.npy'):
            data = genMovingAverage(disp, steps)
            np.save(avgfname, data)
        plotMovingAverage(steps, args.axis)


if __name__ == '__main__':
    main()
