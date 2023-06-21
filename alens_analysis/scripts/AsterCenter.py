import numpy as np
import matplotlib.pyplot as plt
import os
import itertools as itt

import sklearn.cluster as skc  # density based method
import networkx as nx   # connectivity based method

import Util.aLENS as am
import Util.HDF5_Wrapper as h5
import point_cloud.PointCloud as pc


parser = am.getDefaultArgParser('detect aster centers ')
parser.add_argument('--rcut', type=float, dest='rcut', default=-1,
                    help='rcut for gr, default to half box size')
parser.add_argument('--min', type=int, dest='min', default=5,
                    help='minimal number of pts for a cluster')
parser.add_argument('--dbs_eps', type=float, dest='dbs_eps', default=0.1,
                    help='eps for dbscan')
parser.add_argument('--nframe', type=int, dest='nframe', default=10,
                    help='use the last nframe s of data')

args = parser.parse_args()

config = am.parseConfig(args.config)
boxsize = np.array(config['simBoxHigh'])-np.array(config['simBoxLow'])

if args.rcut < 0:
    args.rcut = np.min(0.5*boxsize)

print('rcut: ', args.rcut)
print('config: ', args.config)
print('box_size:', boxsize)

h5FileName = 'AsterCenter_min{:d}'.format(args.min)
h5.newFile(h5FileName)

figpath = 'AsterCenterFigs_min{:d}'.format(args.min)

try:
    os.mkdir(figpath)
except FileExistsError:
    pass


def pbc_replicate(pts=None, boxsize=None,
                  pbc=[True, True, True]):
    px = [-1, 0, 1] if pbc[0] else [0]
    py = [-1, 0, 1] if pbc[1] else [0]
    pz = [-1, 0, 1] if pbc[2] else [0]
    print(px, py, pz)
    rep = list(itt.product(px, py, pz))
    npts = np.shape(pts)[0]
    pbcpts = np.zeros((len(rep)*npts, 3))
    for i in range(len(rep)):
        shift = np.zeros(pts.shape)
        shift[:, 0] = rep[i][0]*boxsize[0]
        shift[:, 1] = rep[i][1]*boxsize[1]
        shift[:, 2] = rep[i][2]*boxsize[2]
        # print(shift)
        pbcpts[i*npts:(i+1)*npts] = shift+pts
    # keep only part
    xcon = np.logical_and(pbcpts[:, 0] > -boxsize[0]
                          * 0.5, pbcpts[:, 0] < boxsize[0]*1.5)
    ycon = np.logical_and(pbcpts[:, 1] > -boxsize[1]
                          * 0.5, pbcpts[:, 1] < boxsize[1]*1.5)
    zcon = np.logical_and(pbcpts[:, 2] > -boxsize[2]
                          * 0.5, pbcpts[:, 2] < boxsize[2]*1.5)

    return pbcpts[np.logical_and(np.logical_and(xcon, ycon), zcon)]


def calc_gr_sq(pts, path, label):
    pc.impose_pbc(pts, boxsize)
    rcut = args.rcut
    npar = len(pts)
    rho = pc.rho(npar, boxsize)
    pairs = pc.get_pair(pts, boxsize, rcut)
    rvec = pc.get_rvec(pts, boxsize, rcut, pairs)
    rnorm = np.linalg.norm(rvec, axis=1)
    r, rdf = pc.gen_rdf(rvec, npar, rho, rcut=rcut, nbins=int(rcut/0.05+1))
    q = np.array([(2*np.pi/boxsize[0])*(j+1) for j in range(500)])
    Sq = np.zeros(q.shape)
    for j in range(len(q)):
        Sq[j] = (pc.Sint3D(q[j], r, rdf, rho))

    RDF = np.vstack([r, rdf]).T
    SQ = np.vstack([q, Sq]).T

    h5.saveData(h5FileName, pts, path+'/'+label, 'pts', float)
    h5.saveData(h5FileName, rnorm, path+'/'+label, 'rnorm', float)
    h5.saveData(h5FileName, RDF, path+'/'+label, 'rdf', float)
    h5.saveData(h5FileName, SQ, path+'/'+label, 'sq', float)

    return RDF, SQ


def plot_gr_sq(grsq1, grsq2, labels, foldername, path):

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for data, label in zip((grsq1, grsq2), labels):
        r = data[0][:, 0]
        rdf = data[0][:, 1]
        q = data[1][:, 0]
        Sq = data[1][:, 1]
        ax1.plot(r, rdf, label=label)
        ax1.set_xlabel('$r / um$')
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlim(0, args.rcut)
        ax2.plot(q, Sq, label=label)
        ax2.set_xlabel(r'$q/ \frac{2\pi k}{L}$')
        ax2.set_ylabel('$S(q)$')
        ax2.set_xlim(0, np.max(q))

    plt.tight_layout()
    name = am.get_basename(path)
    plt.savefig(foldername+os.path.sep+'GrSq_'+name+'.jpg', dpi=300)
    plt.close(fig)
    plt.clf()

    return


def plot_centers(pts, ccs1, ccs2, label1, label2, path):
    '''plot all pts and centers of pts'''
    plt.cla()
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.view_init(15, 30)
    # ax.scatter3D(pts[:, 0], pts[:, 1], pts[:, 2],
    #              marker='+', c='#AAAAAA', alpha=0.1)
    ax.scatter3D(ccs1[:, 0], ccs1[:, 1], ccs1[:, 2],
                 marker='o', c='r', label=label1+' N='+str(len(ccs1)))
    ax.scatter3D(ccs2[:, 0], ccs2[:, 1], ccs2[:, 2],
                 marker='x', c='b', label=label2+' N='+str(len(ccs2)))
    ax.axes.set_xlim3d(left=0, right=boxsize[0])
    ax.axes.set_ylim3d(bottom=0, top=boxsize[1])
    ax.axes.set_zlim3d(bottom=0, top=boxsize[2])
    plt.legend()
    plt.tight_layout()
    name = am.get_basename(path)
    plt.savefig(figpath+os.sep+'AsterCenters'+name+'.jpg', dpi=300)
    return


def ac_dbscan(frame):
    minus_pts = frame.TList[:, 2:5]
    pts = pbc_replicate(minus_pts, boxsize, config['simBoxPBC'])
    clustering = skc.DBSCAN(
        eps=args.dbs_eps, min_samples=args.min).fit(pts)
    result = clustering.labels_
    nc = np.max(result)
    centers = []
    for i in range(0, nc+1):
        idx = (result == i)
        pts_cluster = pts[idx]
        # print for debug
        # print(pts_cluster)
        center = np.mean(pts_cluster, axis=0)
        if center[0] > 0 and center[0] < boxsize[0] \
                and center[1] > 0 and center[1] < boxsize[1] \
                and center[2] > 0 and center[2] < boxsize[2]:
            centers.append(np.mean(pts_cluster, axis=0))

    centers = np.array(centers)
    pc.impose_pbc(coords=centers, boxsize=boxsize)
    print('centers detected: ', len(centers))
    return centers


def ac_nx(frame):
    minus_pts = frame.TList[:, 2:5]
    gid = frame.TList[:, 0]
    pairs = frame.PList
    edges = pairs[np.logical_and(pairs[:, 0] >= 0, pairs[:, 1] >= 0)]
    G = nx.Graph()
    G.add_edges_from(edges)
    subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    centers = []
    for s in subgraphs:
        if len(s) < args.min:
            continue
        gid_cluster = s.nodes()
        idx = np.isin(gid, gid_cluster)
        pts_cluster = minus_pts[idx]
        trg = pts_cluster[0]
        for i in range(1, pts_cluster.shape[1]):
            nbpts = pts_cluster[i]
            for k in range(3):
                if config['simBoxPBC'][k]:
                    m = am.findMove(trg[k], nbpts[k], boxsize[k])
                    nbpts[k] = trg[k]+m
        center = np.mean(pts_cluster, axis=0)
        centers.append(center)

    centers = np.array(centers)
    pc.impose_pbc(coords=centers, boxsize=boxsize)
    print('centers detected: ', len(centers))
    return centers


def process_frame(frame):
    minus_pts = frame.TList[:, 2:5]
    path = am.get_basename(frame.filename)
    ccs1 = ac_dbscan(frame)
    ccs2 = ac_nx(frame)

    ccs1 = ccs1[ccs1[:, 0].argsort()]
    ccs2 = ccs2[ccs2[:, 0].argsort()]
    plot_centers(minus_pts, ccs1, ccs2, 'DBSCAN',
                 'Graph', path)
    data1 = calc_gr_sq(ccs1, path, 'DBSCAN')
    data2 = calc_gr_sq(ccs2, path, 'Graph')
    plot_gr_sq(data1, data2, ['DBSCAN', 'Graph'], figpath, path)
    # np.set_printoptions(precision=6, suppress=True)
    # print(ccs1)
    # print(ccs2)

    return


files = am.getFileListSorted('result*-*/SylinderAscii_*.dat')

for f in files[-args.nframe:]:
    frame = am.FrameAscii(f, readProtein=True, info=False)
    process_frame(frame)
