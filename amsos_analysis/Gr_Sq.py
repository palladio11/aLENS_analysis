import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special as ss
import h5py

import os
import argparse

import AMSOS as am
import point_cloud.PointCloud as pc

boxsize = np.array([10, 10, 10])

data_gr_sq = h5py.File('gr_sq.hdf5', 'w')
data_gr_sq.close()

parser = argparse.ArgumentParser(
    description='Calculate sylinder gr/sq from saved Ascii files.')
parser.add_argument('--rcut', type=float, dest='rcut', default=1.5,
                    help='rcut for gr')

args = parser.parse_args()
# meanWindow = args.window
print('rcut: ', args.rcut)


def calc_gr_sq(points, foldername, filename):
    try:
        os.mkdir(foldername)
    except FileExistsError:
        pass

    pc.impose_pbc(points, boxsize)
    rcut = args.rcut
    npar = len(points)
    rho = pc.rho(npar, boxsize)
    pairs = pc.get_pair(points, boxsize, rcut)
    rvec = pc.get_rvec(points, boxsize, rcut, pairs)
    r, rdf = pc.gen_rdf(rvec, npar, rho, rcut=rcut, nbins=int(rcut/0.005+1))
    # min q = 2*pi/boxsize[0]
    # max q = 500*(min q)
    q = np.array([(2*np.pi/boxsize[0])*(j+1) for j in range(1000)])
    Sq = np.zeros(q.shape)
    for j in range(len(q)):
        Sq[j] = (pc.Sint3D(q[j], r, rdf, rho))

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.semilogy(r/0.025, rdf)
    ax1.set_xlabel('$r/D_{MT}$')
    ax1.set_ylabel('$g(r)$')
    ax1.set_xlim(0, rcut/0.025)
    ax1.set_ylim(0.1, 200)
    ax2.semilogy(q, Sq)
    ax2.set_xlabel(r'$q/ \frac{2\pi k}{L}$')
    ax2.set_ylabel('$S(q)$')
    ax2.set_xlim(0, np.max(q))
    # ax2.set_ylim(0.5, 50)
    plt.tight_layout()
    name = am.get_basename(filename)
    plt.savefig(foldername+os.path.sep+name+'.png', dpi=300)
    plt.close(fig)
    plt.clf()

    return r, rdf, q, Sq


def save_data(grp, name, r, rdf, q, Sq):
    subgrp = grp.create_group(name)
    data1 = np.vstack([r, rdf]).T
    data2 = np.vstack([q, Sq]).T

    dset = subgrp.create_dataset(
        "rdf", data1.shape, dtype='float64', chunks=True)
    dset[...] = data1[...]

    dset = subgrp.create_dataset(
        "Sq", data2.shape, dtype='float64', chunks=True)
    dset[...] = data2[...]
    return


def process_frame(frame):
    minus_pts = frame.TList[:, 2:5]
    plus_pts = frame.TList[:, 5:8]
    center_pts = 0.5*(minus_pts+plus_pts)

    r, rdf, q, Sq = calc_gr_sq(minus_pts, 'rdf_sq_minus', frame.filename)
    data_gr_sq = h5py.File('gr_sq.hdf5', 'a')
    grp = data_gr_sq.create_group(am.get_basename(frame.filename))
    save_data(grp, 'minus', r, rdf, q, Sq)
    data_gr_sq.close()

    r, rdf, q, Sq = calc_gr_sq(plus_pts, 'rdf_sq_plus', frame.filename)
    data_gr_sq = h5py.File('gr_sq.hdf5', 'a')
    grp = data_gr_sq.require_group(am.get_basename(frame.filename))
    save_data(grp, 'plus', r, rdf, q, Sq)
    data_gr_sq.close()

    r, rdf, q, Sq = calc_gr_sq(center_pts, 'rdf_sq_center', frame.filename)
    data_gr_sq = h5py.File('gr_sq.hdf5', 'a')
    grp = data_gr_sq.require_group(am.get_basename(frame.filename))
    save_data(grp, 'center', r, rdf, q, Sq)
    data_gr_sq.close()

    return


files = am.getFileListSorted('result*-*/SylinderAscii_*.dat')

for f in files[-10:]:
    frame = am.FrameAscii(f, readProtein=False, info=False)
    process_frame(frame)
