#!/usr/bin/env python

"""@package docstring
File: time_testing.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import os
from shutil import copy, rmtree
from pathlib import Path
from subprocess import run
from .runlog_funcs import get_walltime, get_wt_timestep, calc_timestep_stats
import toml
import yaml
from pprint import pprint


def run_time_testing(n_time_steps, opts):
    """TODO: Docstring for run_time_testing.
    @return: TODO

    """
    start_dir = Path.cwd()
    # Create mock directory
    tt_path = opts.path / 'time_test'
    if tt_path.exists():
        # Start with a clean directory
        rmtree(tt_path)
    tt_path.mkdir()
    files = [p for p in opts.path.glob('*')
             if p.is_file() and p.suffix != '.zip']
    for f in files:
        copy(f, tt_path)

    try:
        # cd into mock directory to run simulations
        os.chdir(tt_path)

        # Change parameter files
        rc_path = tt_path / 'RunConfig.yaml'
        with rc_path.open('r') as rcf:
            params = yaml.safe_load(rcf)
            omp_num_threads = params.get('omp_num_threads', 1)
            if opts.omp_num_threads:
                print("   ! Overwriting OMP_NUM_THREADS in RunConfig.yaml. !")
                omp_num_threads = opts.omp_num_threads
            print("Number of OMP threads =", omp_num_threads)
            dt = float(params['dt'])
            params['timeTotal'] = dt * float(n_time_steps)
            # Don't write out during this process
            params['timeSnap'] = 10 * dt * float(n_time_steps)
        with rc_path.open('w') as rcf:
            yaml.dump(params, rcf)

        # Run alens (with specified number of
        out_path = (tt_path / 'runlog.out')
        err_path = (tt_path / 'runlog.err')
        run(['./aLENS.X'],
            stdout=out_path.open('w'),
            stderr=err_path.open('w'),
            env=dict(OMP_NUM_THREADS=str(omp_num_threads), **os.environ),
            )

        # Analyze runlog.out for run information
        tot_walltime = get_walltime(out_path)

        stats = {'total walltime': float(tot_walltime.total_seconds())}
        stats['mean step walltime'], stats['median step walltime'], stats['std of step walltime'], stats['max step walltime'] = calc_timestep_stats(
            out_path)
        # Put time analysis file in the analysis directory
        # (create if necessary)
        analysis_dir = start_dir / 'analysis'
        analysis_dir.mkdir(exist_ok=True)
        with (analysis_dir / f'timing_threads_{omp_num_threads}.toml').open('w') as tf:
            toml.dump(stats, tf)
        print("Run time stats")
        pprint(stats, sort_dicts=False)

    except BaseException:
        raise
    finally:
        os.chdir(start_dir)
        # Remove the time testing directory?

    pass


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
