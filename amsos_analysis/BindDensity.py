import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special as ss
import os

import Util.AMSOS as am
import Util.HDF5_Wrapper as h5


parser = am.getDefaultArgParser('Calculate G(r) and S(q)')
parser.add_argument('--pconfig', type=str, dest='pconfig',
                    default='../ProteinConfig.yaml', help='ProteinConfig.yaml file location')
parser.add_argument('--plot_each', type=bool, dest='plot_each',
                    default=False, help='plot hist for each frame')


args = parser.parse_args()

config = am.parseConfig(args.config)
pconfig = am.parseConfig(args.pconfig)

foldername = 'BindDensity'

h5FileName = 'BindDensity'
h5.newFile(h5FileName)


def point_line_proj(point, p0, p1):
    '''find projection of point on p0-p1'''
    u = point-p0
    v = p1-p0
    v_norm = np.sqrt(v.dot(v))
    proj_of_u_on_v = (np.dot(u, v)/v_norm**2)*v
    proj = p0+proj_of_u_on_v  # projection of point to p line
    dist = np.sqrt(proj.dot(proj))
    return proj, dist


def process_frame(frame):
    path = am.get_basename(frame.filename)

    minus_pts = frame.TList[:, 2:5]
    plus_pts = frame.TList[:, 5:8]

    xlinkers = frame.PList
    loc = []
    for xl in xlinkers:
        idBind = xl[-2:]
        if int(idBind[0]) == -1 or int(idBind[1]) == -1:
            continue
        assert int(idBind[0]) == 0 and int(idBind[1]) == 1
        xl0 = xl[2:5]  # position end0
        xl1 = xl[5:8]  # position end1
        proj, dist = point_line_proj(xl0, minus_pts[1], plus_pts[1])
        normT1 = plus_pts[1]-minus_pts[1]
        normT1 = normT1/np.linalg.norm(normT1)
        loc.append((xl1-proj).dot(normT1))
    loc = np.array(loc)
    h5.saveData(h5FileName, loc, path+'/loc', 'loc', float)
    # print(loc)

    name = am.get_basename(frame.filename)

    if not args.plot_each:
        return

    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(111)
    ax.hist(loc, bins=50, range=[-0.075, 0.075],
            density=True)
    plt.savefig(foldername+'/'+name+'.png', dpi=150)
    plt.close(fig)
    plt.clf()

    return


class FrameAscii:
    '''Load Ascii.dat data'''

    def __init__(self, filename):
        self.filename = filename
        self.TList = am.parseSylinderAscii(filename, False, True)
        filename = filename.replace('Sylinder', 'Protein')
        self.PList = np.loadtxt(filename, skiprows=2,
                                usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), dtype=np.float64)


try:
    os.mkdir(foldername)
except FileExistsError:
    pass

files = am.getFileListSorted('result*-*/SylinderAscii_*.dat')

for f in files:
    frame = FrameAscii(f)
    process_frame(frame)
