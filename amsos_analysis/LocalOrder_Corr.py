import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from numpy.lib.recfunctions import structured_to_unstructured
import pyvista as pv

import Util.AMSOS as am


def plot_corr(xdata, ydata, xlabel, ylabel, filename, xlim=None, ylim=None, nbins=20):
    assert xdata.shape == ydata.shape

    hist_norm = mpl.colors.LogNorm()
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
    mesh = pv.read(file)
    print(mesh)
    pgrad = mesh.compute_derivative(scalars="polarity")[
        'gradient']   # 9D per point
    print(pgrad)
    # project pgrad to et and ep directions
    er, et, ep = am.e_sph(
        mesh.points-np.array([100.0, 100.0, 100.0])[np.newaxis, :])
    print(er.shape, et.shape, ep.shape)
    pt = np.vstack([np.einsum('ij, ij->i', pgrad[:, :3], et),
                    np.einsum('ij, ij->i', pgrad[:, 3:6], et),
                    np.einsum('ij, ij->i', pgrad[:, 6:], et)
                    ])
    pp = np.vstack([np.einsum('ij, ij->i', pgrad[:, :3], ep),
                    np.einsum('ij, ij->i', pgrad[:, 3:6], ep),
                    np.einsum('ij, ij->i', pgrad[:, 6:], ep)
                    ])
    pgrad_s = (pt+pp).T

    plot_corr(xdata=mesh['xlinker_n_db'],
              #   ydata=np.linalg.norm(mesh['polarity'], axis=1),
              ydata=np.linalg.norm(pgrad_s, axis=1),
              xlabel='number density of doubly bound $1/um^3$',
              ylabel='grad_s p',
              filename=file,
              )
    return


files = am.getFileListSorted('./*.vtu')
for f in files[900:905]:
    calc_corr(f)
