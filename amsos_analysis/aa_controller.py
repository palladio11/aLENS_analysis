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

from .read_func import convert_dat_to_hdf
from .hic_animation import hic_animation


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
                        "(overwrite, analyze(default), or load.")
    parser.add_argument("-M", "--movie", action="store_true", default=False,
                        help=("Create an animation from a seed."))
    parser.add_argument("-G", "--graph", action="store_true", default=False,
                        help=("Create graph of a seed's end state."))

    opts = parser.parse_args()

    opts.params = {
        'n_graph': 1,
        'fps': 25,
        'style': "log_contact"
    }

    # Post parsing changes to options
    opts.path = Path(opts.path).resolve()
    opts.data_dir = opts.path / 'result'

    return opts


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
        # run_analysis(opts)

    if opts.movie:
        hic_animation(opts)
        return

    # if opts.graph:
        # make_graphs(opts)
    h5_data.flush()
    h5_data.close()


##########################################
if __name__ == "__main__":
    main()
