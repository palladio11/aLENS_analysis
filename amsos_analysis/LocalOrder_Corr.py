import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import structured_to_unstructured
import pyvista as pv

import Util.AMSOS as am

parser = am.getDefaultArgParser('compute correlation')
parser.add_argument('--nbins', type=int, dest='nbins', default=40,
                    help='number of hist bins')

args = parser.parse_args()


def plot_corr(xdata, ydata, xlabel, ylabel, filename, xlim=None, ylim=None, nbins=20):
    assert xdata.shape == ydata.shape

    hist_norm = mpl.colors.LogNorm(vmin=1e-4, vmax=1)
    ax = plt.subplot(111)
    h = ax.hist2d(xdata, ydata,
                  density=True, norm=hist_norm, bins=nbins)
    plt.colorbar(h[3], label='Probability Density', ax=ax)
    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim != None:
        ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(filename)
    plt.savefig(am.get_basename(filename)+'.jpg', dpi=300)
    plt.clf()
    plt.cla()
    return


def calc_corr(file):
    print(file)
    mesh = pv.read(file)
    # print(mesh)
    pgrad = mesh.compute_derivative(scalars="polarity")['gradient']
    assert pgrad.shape[1] == 9
    # print(pgrad)
    # # project pgrad to et and ep directions
    # er, et, ep = am.e_sph(
    #     mesh.points-np.array([100.0, 100.0, 100.0])[np.newaxis, :])
    # print(er.shape, et.shape, ep.shape)
    # pt = np.vstack([np.einsum('ij, ij->i', pgrad[:, :3], et),
    #                 np.einsum('ij, ij->i', pgrad[:, 3:6], et),
    #                 np.einsum('ij, ij->i', pgrad[:, 6:], et)
    #                 ])
    # pp = np.vstack([np.einsum('ij, ij->i', pgrad[:, :3], ep),
    #                 np.einsum('ij, ij->i', pgrad[:, 3:6], ep),
    #                 np.einsum('ij, ij->i', pgrad[:, 6:], ep)
    #                 ])
    # pgrad_s = (pt+pp).T
    pdiv = pgrad[:, 0]+pgrad[:, 4]+pgrad[:, 8]

    n_ref = 200000/(4.0*np.pi/3*(5.102**3-5.0**3))
    plot_corr(xdata=mesh['xlinker_n_db']/n_ref,
              ydata=pdiv,
              xlabel=r'number density of doubly bound $n/n_{\rm ave}$',
              ylabel='div p',
              filename=file,
              xlim=[0, 4],
              ylim=[-5, 5],
              nbins=args.nbins
              )
    return


files = am.getFileListSorted('./*.vtu', info=False)
for f in files[::50]:
    calc_corr(f)
