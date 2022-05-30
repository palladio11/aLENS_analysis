#!/usr/bin/env python

"""@package docstring
File: hic_animation.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path
from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation, FFMpegWriter
import math
import yaml
from numba import jit, vectorize
import argparse

from ..objects import filament
from ..read_func import (read_dat_sylinder,
                         get_file_number,
                         get_png_number,
                         count_fils)


SQRT2 = np.sqrt(2)


def make_separation_mat(com_arr, downsample=1):
    nbeads = com_arr.shape[0]
    reduc_com_arr = com_arr[::downsample]
    x = np.arange(nbeads + 1)[::int((nbeads) / reduc_com_arr.shape[0])]
    X, Y = np.meshgrid(x, x)
    dist_mat = np.linalg.norm(
        reduc_com_arr[:, np.newaxis, :] - reduc_com_arr[np.newaxis, :, :],
        axis=-1)
    return dist_mat, X, Y


def gauss_weighted_contact(sep_mat, sigma=.020):
    return np.exp(-np.power(sep_mat, 2) / (2. * (sigma * sigma)))


def create_hic_frame(fil_dat_path, style='sep',
                     downsample=1, bead_range=None, **kwargs):
    # Get filament data
    fils = read_dat_sylinder(fil_dat_path)
    com_arr = np.asarray([fil.get_com()
                          for fil in fils if (fil.fil_type != 'L')])
    if not bead_range is None:
        if len(bead_range) == 1:
            com_arr = com_arr[bead_range[0]:]
        else:
            com_arr = com_arr[bead_range[0]:bead_range[1]]

    sep_mat, X, Y = make_separation_mat(com_arr, downsample)

    if style == 'sep':
        return sep_mat, X, Y
    if style == 'contact':
        return gauss_weighted_contact(sep_mat), X, Y
    if style == 'log_contact':
        return -.5 * np.power(sep_mat, 2) / (.02 * .02 * np.log(10)), X, Y
    else:
        raise RuntimeError(f' The style "{style}" is not supported currently.')


def animate(i, fig, axarr, fil_dat_paths, png_paths, init_mutable, opts):
    for ax in axarr:
        ax.clear()

    png = plt.imread(str(png_paths[i]))
    img = axarr[0].imshow(png, resample=False)
    axarr[0].set_axis_off()

    print(f'Making frame {i}')

    t0 = time()
    frame, X, Y = create_hic_frame(fil_dat_paths[i], **opts.params)
    t1 = time()
    print(f"Frame {i} created in: {t1-t0:.2g} sec")
    # c = axarr[1].pcolorfast(frames[i], vmax=vmax, vmin=0)

    c = axarr[1].pcolorfast(X, Y, frame, vmin=opts.params['vmin'])
    t2 = time()
    print(f"Frame {i} drawn in: {t2 - t1:.2g} sec")
    if init_mutable[0] == True:  # Cludge, there is a better way to do thisdf cdf
        fig.colorbar(
            c,
            ax=axarr[1],
            # label=r"$\log$(Inferred contact map) $\sim$
            # ($r_{ij}^2/2\sigma^2$)")
            label=r"Log contact probability $\sim$ ($-r_{ij}^2$)")
        init_mutable[0] = False

    axarr[1].set_xlabel(r"Bead index")
    axarr[1].set_ylabel(r"Bead index")

    # pcm = ax.pcolorfast(frames[i], cmap='gray', vmax=vmax)
    axarr[0].set_title("Time {:.2f} sec".format(
        float(i * opts.params['time_step'] * opts.params['n_graph'])))
    return [img, c]


def hic_animation(opts):

    alens_stl = {
        "axes.titlesize": 20,
        "axes.labelsize": 24,
        "lines.linewidth": 3,
        "lines.markersize": 10,
        "xtick.labelsize": 24,
        "ytick.labelsize": 24,
        "font.size": 20,
        "font.sans-serif": 'Helvetica',
        "text.usetex": False,
        'mathtext.fontset': 'cm',
    }
    plt.style.use(alens_stl)

    with open(opts.path / 'RunConfig.yaml', 'r') as yf:
        run_params = yaml.safe_load(yf)
        opts.params['time_step'] = run_params['timeSnap']

    result_dir = opts.result_dir
    fil_dat_paths = sorted(result_dir.glob("**/SylinderAscii*.dat"),
                           key=get_file_number)[::opts.params['n_graph']]
    png_paths = sorted(result_dir.glob("PNG/*.png"),
                       key=get_png_number)[::opts.params['n_graph']]

    # print(png_paths)
    init_mutable = [True]
    nframes = len(png_paths)
    print(nframes)
    fig, axarr = plt.subplots(1, 2, figsize=(20, 8))
    # fig = plt.figure(figsize=(20, 8))
    # gs = GridSpec(1, 3)
    # axarr = [fig.add_subplot(gs[:-1])]
    # axarr += [fig.add_subplot(gs[-1])]

    axarr[1].set_aspect('equal')
    writer = FFMpegWriter(
        fps=opts.params['fps'],
        codec='libx264',
        bitrate=-1,
        extra_args=[
            '-pix_fmt',
            'yuv420p'])
    vmax = 0
    ani = FuncAnimation(fig, animate, nframes,
                        fargs=(
                            fig,
                            axarr,
                            fil_dat_paths,
                            png_paths,
                            init_mutable, opts),
                        blit=True)
    ani.save(
        opts.analysis_dir /
        f'{opts.params["style"]}_mat_vid.mp4',
        writer=writer)


##########################################
if __name__ == "__main__":
    pass
    # opts = parse_args()
    # hic_animation(opts)
