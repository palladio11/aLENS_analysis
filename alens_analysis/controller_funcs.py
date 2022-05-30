#!/usr/bin/env python

"""@package docstring
File: chrom_controller_funcs.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import re
import h5py
import time
import shutil

from datetime import datetime
from .chromatin.chrom_graph_funcs import (make_all_condensate_graphs)
from .chromatin.chrom_seed_scan_graph_funcs import (
    make_all_seed_scan_condensate_graphs)
from .read_func import convert_dat_to_hdf
from .chromatin.hic_animation import hic_animation
from .min_animation import min_animation
from .result_to_pvd import make_pvd_files
from .runlog_funcs import get_walltime

MOVIE_DICT = {'hic': hic_animation,
              'min': min_animation,
              }


def make_seed_graphs(opts):
    """TODO: Docstring for make_graphs.

    @param opts TODO
    @return: TODO

    """
    opts.analysis_dir.mkdir(exist_ok=True)
    h5_path = opts.analysis_dir / f'{opts.path.stem}.h5'

    with h5py.File(h5_path, 'a') as h5_data:
        overwrite = True if opts.analysis == 'overwrite' else False
        make_all_condensate_graphs(h5_data, opts, overwrite=overwrite)


def make_seed_scan_graphs(opts):
    """TODO: Docstring for make_graphs.

    @param opts TODO
    @return: TODO

    """
    opts.analysis_dir.mkdir(exist_ok=True)
    h5_path = opts.analysis_dir / f'{opts.path.stem}.h5'

    with h5py.File(h5_path, 'a') as h5_scan_data:
        overwrite = True if opts.analysis == 'overwrite' else False
        try:
            sd_h5_data_lst = [
                h5py.File(sdh5, 'r')
                for sdh5 in opts.result_dir.glob('s*/analysis/*.h5')]
            make_all_seed_scan_condensate_graphs(
                h5_scan_data, sd_h5_data_lst, opts, overwrite=overwrite)
        except BaseException:
            raise
        finally:
            for sdh5 in sd_h5_data_lst:
                sdh5.close()


def seed_analysis(opts):
    """ All subprocesses that go into analyzing a single seed

    @param opts Parsed options object
    @return: void, hdf5 files, graphs, and movies maybe saved in analysis
             directory

    """
    if opts.analysis == 'collect':
        # make_pvd_files(opts.result_dir)
        t0 = time.time()
        h5_path = opts.analysis_dir / f'{opts.path.stem}.h5'
        print(f'{opts.path.stem}')
        h5_data = convert_dat_to_hdf(h5_path, opts.path)
        print(f" HDF5 created in {time.time() - t0}")
        # Check to see if run.log is present in current simulation
        # Wall time analysis
        if (opts.path / 'run.log').exists():
            dwtime = get_walltime(opts.path / 'run.log')
            h5_data.attrs['total_seconds'] = dwtime.total_seconds()
            h5_data.attrs['walltime'] = str(dwtime)
        try:
            h5_data.flush()
            h5_data.close()
        except BaseException:
            print("Could not close h5_data file")
            pass

    if opts.movie:
        MOVIE_DICT[opts.movie](opts)
        return

    if opts.graph:
        t0 = time.time()
        make_seed_graphs(opts)
        print(f" Graphs created in {time.time() - t0}")


def seed_scan_analysis(opts):
    """ All subprocesses that go into analyzing a single parameter set with
    multiple random number seeds

    @param opts Parsed options object
    @return: void, hdf5 files, graphs, and movies maybe saved in analysis
             directory

    """
    # Store a simulation directory
    opts.result_dir = opts.path / 'simulations'
    if opts.analysis == 'collect':
        for sd in opts.result_dir.glob('s*'):
            sd_copy = opts.analysis_dir / sd.stem
            # Collect all graphs and images and put in new directories
            shutil.copytree(sd / 'analysis', sd_copy,
                            ignore=shutil.ignore_patterns('*.h5'),
                            dirs_exist_ok=True)

        # TODO Create a master hdf5 file for all of these for easier

    elif opts.graph:
        t0 = time.time()
        make_seed_scan_graphs(opts)
        print(f" Graphs created in {time.time() - t0}")


def param_scan_analysis(opts):
    """ All subprocesses that go into analyzing a single parameter set with
    multiple random number seeds

    @param opts Parsed options object
    @return: void, hdf5 files, graphs, and movies maybe saved in analysis
             directory

    """
    opts.result_dir = opts.path / 'simulations'
    if opts.analysis:
        if opts.analysis == 'collect':
            for pd in opts.result_dir.glob('*/'):
                sd_copy = opts.analysis_dir / pd.stem
                # Collect all graphs and images and put in new directories
                shutil.copytree(pd / 'result/analysis', sd_copy,
                                ignore=shutil.ignore_patterns('*.h5'),
                                dirs_exist_ok=True)
            # TODO Create a master hdf5 file for all of these for easier
            # Collect all hdf5 files and put into new directory
            pass
    pass


def param_seed_scan_analysis(opts):
    """ All subprocesses that go into analyzing a single parameter set with
    multiple random number seeds

    @param opts Parsed options object
    @return: void, hdf5 files, graphs, and movies maybe saved in analysis
             directory

    """
    if opts.analysis:
        if opts.analysis == 'collect':
            # Collect all hdf5 files and put into new directory
            pass
    pass


TYPE_FUNC_DICT = {'seed': seed_analysis,
                  'seed_scan': seed_scan_analysis,
                  'param_scan': param_scan_analysis,
                  'param_seed_scan': param_seed_scan_analysis,
                  }

##########################################
if __name__ == "__main__":
    print("Not implemented yet")
