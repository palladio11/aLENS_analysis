#!/usr/bin/env python

"""@package docstring
File: aa_controller.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import argparse
from pathlib import Path

import h5py
import yaml
import time
import re
import numpy as np
from datetime import datetime

import matplotlib.pyplot as plt

from .read_func import convert_dat_to_hdf
from .hic_animation import hic_animation
from .min_animation import min_animation
from .colormaps import register_cmaps
from .chrom_analysis import get_pos_kymo_data, get_pos_cond_data
from .chrom_graph_funcs import make_all_condensate_graphs

MOVIE_DICT = {'hic': hic_animation,
              'min': min_animation,
              }


def parse_args():
    parser = argparse.ArgumentParser(
        prog='aa_controller.py',
        formatter_class=argparse.RawTextHelpFormatter)

    # TODO: Implement subparser if necessary <26-02-21, ARL> #
    # subparsers = parser.add_subparsers(title="Simulation types")
    # chrom_parser = subparsers.add_parser

    parser.add_argument("-p", "--path", default=".",
                        help="Path used in AMSOS Analysis functions.")
    parser.add_argument('-i', "--image_input", default=None,
                        help="Image parameter yaml file")

    parser.add_argument("-t", "--type", default=None,
                        help="Type of analysis")

    parser.add_argument("-A ", "--analysis",
                        choices=[
                            None,
                            'collect',
                            'read',
                            'load',
                            'analyze',
                            'overwrite'],
                        default=None,
                        help=" Specify analysis type to determine if data will"
                        " be overwritten. Options include "
                        "(overwrite, None(default), or load.")
    parser.add_argument("-M", "--movie", choices=[None, "hic", "min"], default=None,
                        help=("Create an animation from a seed. "
                              "hic: movie with instantaneous Hi-C map"
                              "min: images only"))
    parser.add_argument("-G", "--graph", choices=[None, "condense"], default=False,
                        help=("Create graph of a seed's end state. "
                              "condense: graphs related to condensates"))

    parser.add_argument("-cm", "--colormap", default=None,
                        help=("Specify a colormap to use in graphs"))

    opts = parser.parse_args()

    if opts.colormap:
        register_cmaps()
        plt.rcParams['image.cmap'] = opts.colormap

    opts.params = {
        'n_graph': 5,
        'fps': 20,
        # 'n_graph': 1,
        # 'fps': 2,
        'style': "log_contact",
        'downsample': 1,
        'vmin': -20,
        # 'bead_range': [0, 520],
        'bead_range': None,
    }

    # Post parsing changes to options
    opts.path = Path(opts.path).resolve()
    opts.data_dir = opts.path / 'result'
    opts.analysis_dir = opts.data_dir / 'analysis'
    print(opts.path)

    return opts


def get_walltime(log_path):
    """Uses the log file to calculate the total time the simulation took.
    This will not work for restarted simulations. Might want to fix that.

    @param log_path TODO
    @return: TODO

    """
    with open(log_path, 'r') as rlf:
        pattern = re.compile(r'\[(\d+-\d+-\d+\s\d+:\d+:\d+\.\d+)\]')
        line = rlf.readline()
        while not pattern.search(line):
            line = rlf.readline()
        start_wtime = pattern.search(line).group(0)

        for line in reversed(rlf.readlines()):
            if not pattern.search(line):
                continue
            end_wtime = pattern.search(line).group(0)
            break

    stripstr = '[%Y-%m-%d %H:%M:%S.%f]'
    end_dt = datetime.strptime(end_wtime, stripstr)
    start_dt = datetime.strptime(start_wtime, stripstr)
    return end_dt - start_dt


def make_graphs(opts):
    """TODO: Docstring for make_graphs.

    @param opts TODO
    @return: TODO

    """
    h5_path = opts.path / f'{opts.path.stem}.h5'
    opts.analysis_dir.mkdir(exist_ok=True)

    with h5py.File(h5_path, 'r+') as h5_data:
        make_all_condensate_graphs(h5_data, opts, opts.analysis)


def main():
    """!Main function for AMSOS analysis controller

    Parameters
    ----------
    Will decipher command line arguments

    Returns
    -------
    Void

    """
    opts = parse_args()
    if opts.analysis:
        if opts. analysis == 'collect':
            t0 = time.time()
            print(f'{opts.path.stem}')
            h5_data = convert_dat_to_hdf(f'{opts.path.stem}.h5', opts.path)
            print(f" HDF5 created in {time.time() - t0}")
            # Check to see if run.log is present in current simulation
            # Wall time analysis
            if (opts.path / 'run.log').exists():
                dwtime = get_walltime(opts.path / 'run.log')
                h5_data.attrs['total_seconds'] = dwtime.total_seconds()
                h5_data.attrs['walltime'] = str(dwtime)

    if opts.movie:
        MOVIE_DICT[opts.movie](opts)
        return

    if opts.graph:
        t0 = time.time()
        make_graphs(opts)
        print(f" Graphs created in {time.time() - t0}")

    # if opts.graph:
    try:
        h5_data.flush()
        h5_data.close()
    except BaseException:
        pass


##########################################
if __name__ == "__main__":
    main()
