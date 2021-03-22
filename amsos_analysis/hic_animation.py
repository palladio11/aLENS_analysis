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

rng = np.random.default_rng()  # initialize generator instance

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


def make_separation_mat(com_arr):
    row_extend_arr = np.tile(np.expand_dims(
        com_arr, axis=0), (com_arr.shape[0], 1, 1))
    col_extend_arr = np.transpose(row_extend_arr, (1, 0, 2))
    dist_mat = np.linalg.norm(row_extend_arr - col_extend_arr, axis=-1)

    return dist_mat


def gauss_weighted_contact(sep_mat, sigma=.020):
    return np.exp(-np.power(sep_mat, 2) / (2. * (sigma * sigma)))


def create_hic_frame(fil_dat_path, style='sep', **kwargs):
    # Get filament data
    fils = read_dat_sylinder(fil_dat_path)
    com_arr = []
    for fil in fils:
        fil.parse()
        com_arr += [fil.get_com()]

    com_arr = np.asarray(com_arr)
    # com_arr = np.asarray([fil.get_com() for fil in fils])
    sep_mat = make_separation_mat(com_arr)
    if style == 'sep':
        return sep_mat
    if style == 'contact':
        return gauss_weighted_contact(sep_mat)
    if style == 'log_contact':
        print
        return -.5 * np.power(sep_mat, 2) / (.02 * .02 * np.log(10))
    else:
        raise RuntimeError(f' The style "{style}" is not supported currently.')


def animate(i, fig, axarr, frames, png_paths, vmax, init_mutable, opts):
    for ax in axarr:
        ax.clear()

    png = plt.imread(png_paths[opts.params['n_graph'] * i])
    img = axarr[0].imshow(png)

    print(f'Making frame {i}')

    # c = axarr[1].pcolorfast(frames[i], vmax=vmax, vmin=0)
    c = axarr[1].pcolorfast(frames[i], vmin=-20)
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

    with open('../RunConfig.yaml', 'r') as yf:
        run_params = yaml.safe_load(yf)

    result_dir = Path(".")
    fil_dat_paths = sorted(result_dir.glob("**/SylinderAscii*.dat"),
                           key=get_file_number)
    png_paths = sorted(result_dir.glob("PNG/*.png"),
                       key=get_png_number)

    # rng = np.random.default_rng()

    init_mutable = [True]
    frames = []
    for i, fdp in enumerate(fil_dat_paths[::opts.params['n_graph']]):
        t0 = time()
        frames += [create_hic_frame(fdp, **opts.params)]
        print("Frame {} created in: {:.2g} sec".format(i, time() - t0))
    fig, axarr = plt.subplots(1, 2, figsize=(18, 8))
    # axarr[0].set_aspect('equal')
    axarr[1].set_aspect('equal')
    vmax = np.amax(frames)
    ani = FuncAnimation(fig, animate, len(frames),
                        fargs=(
        fig,
        axarr,
        frames,
        png_paths,
        vmax,
        init_mutable, opts),
        blit=True)
    writer = FFMpegWriter(fps=opts.params['fps'], bitrate=1800, )
    ani.save("log_contact_mat_vid.mp4", writer=writer)


##########################################
if __name__ == "__main__":
    opts = parse_args()
    hic_animation(opts)
