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
import numpy as np

import matplotlib.pyplot as plt

from .colormaps import register_cmaps
from .controller_funcs import TYPE_FUNC_DICT
# from .chrom_analysis import get_pos_kymo_data, get_pos_cond_data


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

    parser.add_argument("-t", "--type",
                        choices=[
                            "seed",
                            "seed_scan",
                            "param_scan",
                            "param_seed_scan"],
                        default="seed",
                        help="Type of analysis. Options"
                        " seed: single simulation analysis (default)\n"
                        " seed_scan: analysis of data from multiple seeds\n"
                        " param_scan: comparison of data from single seed\n"
                        "             simulations with different parameters\n"
                        " param_seed_scan: comparison of data of different\n"
                        "                  parameter with multiple seeds\n"
                        )

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
    parser.add_argument("-G", "--graph", choices=[None, "condense"], default=None,
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
    print(opts.path)

    opts.result_dir = opts.path / 'result'
    opts.analysis_dir = opts.path / 'analysis'
    opts.analysis_dir.mkdir(exist_ok=True)

    return opts


def main():
    """!Main function for aLENs analysis controller

    Parameters
    ----------
    Will decipher command line arguments

    Returns
    -------
    Void

    """
    opts = parse_args()
    TYPE_FUNC_DICT[opts.type](opts)


##########################################
if __name__ == "__main__":
    main()
