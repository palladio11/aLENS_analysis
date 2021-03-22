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


def parse_args():
    parser = argparse.ArgumentParser(
        prog='aa_controller.py',
        formatter_class=argparse.RawTextHelpFormatter)

    # TODO: Implement subparser if necessary <26-02-21, ARL> #
    # subparsers = parser.add_subparsers(title="Simulation types")
    # chrom_parser = subparsers.add_parser

    parser.add_argument("-p", "--path", default=".",
                        help="Path used in AMSOS Analysis functions.")
    parser.add_argument('-i', "--input", default=None,
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

    # parser.add_argument('-i', "--input", default=None,
    #                     help="Image parameter yaml file")
    # if opts.input is None:
    #     opts.params = {
    #         'n_graph': 10,
    #         'fps': 25,
    #         'time_step': 1.,
    #         'style': "log_contact"
    #     }
    # else:
    #     param_path = Path(opts.input)
    #     if not param_path.exists():
    #         raise IOError(
    #             " {} does not exist. Put in valid path.".format(param_path))

    #     with param_path.open('r') as pf:
    #         opts.params = yaml.safe_load(pf)

    # parser.add_argument(
    #     "-r", "--run_type", type=str,
    #     choices=['single_seed',
    #              'multi_seed',
    #              'param_scan',
    #              'param_single_scan'],
    #     default="single_seed",
    #     help=(
    #         "AMSOS can analyze multiple simulations and aggregate "
    #         "data according to various schemes. The 'run_type' argument "
    #         "specifies how to collect the data from nested data directories.\n"
    #     ))

    # parser.add_argument(
    #     "-p", "--param", type=str, default=None,
    #     help=("Parameter to use in analysis or graphing functions."))
    # parser.add_argument("--spec", type=str, default='',
    #                     help=("Specify if parameter used in param scan or "
    #                           "full run anlaysis is a specific species "
    #                           "parameter. e.g. 'crosslink' "))
    opts = parser.parse_args()

    # Post parsing changes to options
    opts.path = Path(opts.path).resolve()
    opts.data_dir = opts.path / 'result'

    return opts


def main():
    """!Main function for simcore_analysis

    Parameters
    ----------
    param_file: Yaml Parameter file

    Returns
    -------
    TODO

    """
    opts = parse_args()
    if opts.analysis:
        if opts. analysis == 'collect':
            t0 = time.time()
            print(f'{opts.path.stem}')
            h5_data = convert_dat_to_hdf(f'{opts.path.stem}.h5', opts.path)
            print(f" HDF5 created in {time.time() - t0}")
        # run_analysis(opts)

    # if opts.graph:
        # make_graphs(opts)
    h5_data.flush()
    h5_data.close()


##########################################
if __name__ == "__main__":
    main()
