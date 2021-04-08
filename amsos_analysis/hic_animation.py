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
from matplotlib.animation import FuncAnimation, FFMpegWriter
import math
import yaml
from numba import jit, vectorize
import argparse

from .objects import filament
from .read_func import (read_dat_sylinder,
                        get_file_number,
                        get_png_number,
                        count_fils)


SQRT2 = np.sqrt(2)


def parse_args():

    parser = argparse.ArgumentParser(prog='hic_animation.py')
    parser.add_argument('-i', "--input", default=None,
                        help="Image parameter yaml file")
    opts = parser.parse_args()

    if opts.input is None:
        opts.params = {
            'n_graph': 10,
            'fps': 25,
            'time_step': 1.,
            'style': "log_contact"
        }
    else:
        param_path = Path(opts.input)
        if not param_path.exists():
            raise IOError(
                " {} does not exist. Put in valid path.".format(param_path))

        with param_path.open('r') as pf:
            opts.params = yaml.safe_load(pf)

    return opts


def make_separation_mat(com_arr, downsample=1):
    reduc_com_arr = com_arr[::downsample]
    dist_mat = np.linalg.norm(
        reduc_com_arr[:, np.newaxis, :] - reduc_com_arr[np.newaxis, :, :],
        axis=-1)
    return dist_mat


def gauss_weighted_contact(sep_mat, sigma=.020):
    return np.exp(-np.power(sep_mat, 2) / (2. * (sigma * sigma)))


def create_hic_frame(fil_dat_path, style='sep', downsample=1, **kwargs):
    # Get filament data
    fils = read_dat_sylinder(fil_dat_path)
    com_arr = np.asarray([fil.get_com() for fil in fils])
    sep_mat = make_separation_mat(com_arr, downsample)

    if style == 'sep':
        return sep_mat
    if style == 'contact':
        return gauss_weighted_contact(sep_mat)
    if style == 'log_contact':
        return -.5 * np.power(sep_mat, 2) / (.02 * .02 * np.log(10))
    else:
        raise RuntimeError(f' The style "{style}" is not supported currently.')


def animate(i, fig, axarr, fil_dat_paths, png_paths, init_mutable, opts):
    for ax in axarr:
        ax.clear()

    png = plt.imread(str(png_paths[i]))
    img = axarr[0].imshow(png)

    print(f'Making frame {i}')

    t0 = time()
    frame = create_hic_frame(fil_dat_paths[i], **opts.params)
    t1 = time()
    print(f"Frame {i} created in: {t1-t0:.2g} sec")
    # c = axarr[1].pcolorfast(frames[i], vmax=vmax, vmin=0)

    c = axarr[1].pcolorfast(frame, vmin=opts.params['vmin'])
    t2 = time()
    print(f"Frame {i} drawn in: {t2 - t1:.2g} sec")
    if init_mutable[0] == True:  # Cludge, there is a better way to do thisdf cdf
        fig.colorbar(c, ax=axarr[1], label=r"Log prob contact")
        init_mutable[0] = False

    axarr[1].set_xlabel("Loci bin")
    axarr[1].set_ylabel("Loci bin")

    # pcm = ax.pcolorfast(frames[i], cmap='gray', vmax=vmax)
    axarr[0].set_title("Time {:.2f} sec".format(float(i *
                                                      opts.params['time_step'] *
                                                      opts.params['n_graph'])))
    return [img, c]


def hic_animation(opts):

    amsos_stl = {
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
    plt.style.use(amsos_stl)

    with open(opts.path / 'RunConfig.yaml', 'r') as yf:
        run_params = yaml.safe_load(yf)
        opts.params['time_step'] = run_params['timeSnap']

    result_dir = opts.data_dir
    # print(result_dir)
    fil_dat_paths = sorted(result_dir.glob("**/SylinderAscii*.dat"),
                           key=get_file_number)[::opts.params['n_graph']]
    # print(fil_dat_paths)
    png_paths = sorted(result_dir.glob("PNG/*.png"),
                       key=get_png_number)[::opts.params['n_graph']]

    # print(png_paths)
    init_mutable = [True]
    nframes = len(fil_dat_paths)
    print(nframes)
    # for i, fdp in enumerate(fil_dat_paths[::opts.params['n_graph']]):
    #     t0 = time()
    #     frames += [create_hic_frame(fdp, **opts.params)]
    #     print("Frame {} created in: {:.2g} sec".format(i, time() - t0))
    fig, axarr = plt.subplots(1, 2, figsize=(20, 8))
    # axarr[0].set_aspect('equal')
    axarr[1].set_aspect('equal')
    # vmax = np.amax(frames)
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
    ani.save(result_dir / f'{opts.params["style"]}_mat_vid.mp4', writer=writer)


##########################################
if __name__ == "__main__":
    opts = parse_args()
    hic_animation(opts)
