#!/usr/bin/env python

"""@package docstring
File: hic_animation.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path
from time import time
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import FuncAnimation, FFMpegWriter
import yaml

from .objects import filament
from .read_func import (read_dat_sylinder,
                        get_file_number,
                        get_png_number,
                        count_fils)


def animate(i, fig, ax, png_paths, opts):
    ax.clear()

    png = plt.imread(str(png_paths[i]))
    img = ax.imshow(png, resample=False)
    ax.set_axis_off()

    print(f'Making frame {i}')

    ax.set_title("Time {:.2f} sec".format(
        float(i * opts.params['time_step'] * opts.params['n_graph'])))
    return [img]


def min_animation(opts):

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
    png_paths = sorted(result_dir.glob("PNG/*.png"),
                       key=get_png_number)[::opts.params['n_graph']]

    nframes = len(png_paths)
    print(nframes)
    fig, ax = plt.subplots(figsize=(20, 8))
    writer = FFMpegWriter(
        fps=opts.params['fps'],
        codec='libx264',
        # bitrate=20000,
        extra_args=[
            # '-vcodec', 'libx264',
            '-pix_fmt', 'yuv420p',
        ]
    )
    vmax = 0
    ani = FuncAnimation(fig, animate, nframes,
                        fargs=(
                            fig,
                            ax,
                            png_paths,
                            opts),
                        blit=True)
    ani.save(opts.analysis_dir / f'min_vid.mp4', writer=writer)


##########################################
if __name__ == "__main__":
    pass
