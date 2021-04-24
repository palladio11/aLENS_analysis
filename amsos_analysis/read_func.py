#!/usr/bin/env python

"""@package docstring
File: read_func.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
import h5py
import re
import yaml
import numpy as np
from pathlib import Path
from .objects import filament, protein


def get_file_number(path):
    name = path.stem
    num = name.split("_")[-1]
    return int(num)


def get_png_number(path):
    name = path.stem
    num = name.split(".")[1]
    return int(num)


def count_fils(path):
    with path.open('r') as pf:
        for i, l in enumerate(pf, -1):  # Don't count first two lines
            pass
    return i


def read_dat_sylinder(fpath):
    # Read a SylinderAscii_X.dat file

    # open the file and read the lines
    with fpath.open('r') as file1:
        filecontent = file1.readlines()

        # Delete the first two lines because they dont have any data
        filecontent[0:2] = []

        # Create list of filaments
        filaments = sorted([filament(line)
                            for line in filecontent], key=lambda x: int(x.gid))
    return filaments


def read_dat_xlp(fpath):
    # Read a ProteinAscii_X.dat file

    # open the file and read the lines
    with fpath.open('r') as file1:
        filecontent = file1.readlines()

        # Delete the first two lines because they dont have any data
        filecontent[0:2] = []

        # Create list of proteins
        filaments = sorted([protein(line)
                            for line in filecontent], key=lambda x: int(x.gid))
    return filaments


def read_time(fnames, h5_data):
    """!Read in data from all protein files

    @param fnames: List posit file names
    @param h5_data: HDF5 position data gropu
    @return: HDF5 data set containing protein data

    """
    t = []
    for fn in fnames:
        p = Path(fn).resolve()
        with p.open(mode='r') as f:
            t += [float(f.readlines()[1])]
    time_dset = h5_data.create_dataset('time', data=t)
    return time_dset


def read_protein_data(protein_fnames, posit_grp):
    """!Read in data from all protein files

    @param protein_fnames: List of protein posit file names
    @param posit_grp: HDF5 position data group
    @return: HDF5 data set containing protein data

    """
    nframes = len(protein_fnames)
    p = Path(protein_fnames[0]).resolve()
    with p.open(mode='r') as fn:
        nproteins = int(fn.readline())

    # Create dataset for protein info (input shape)
    protein_dset = posit_grp.create_dataset('proteins',
                                            shape=(nproteins, 10, nframes),
                                            )
    # Set protein attribute: names of columns
    protein_dset.attrs['nproteins'] = nproteins
    protein_dset.attrs['axis dimensions'] = ['protein', 'state', 'frame']
    protein_dset.attrs['column labels'] = ['gid', 'tag',
                                           'end1 pos x', 'end1 pos y', 'end1 pos z',
                                           'end2 pos x', 'end2 pos y', 'end2 pos z',
                                           'end1 bindID', 'end2 bindID']
    # Loop over files adding to h5_data
    for frame, pfname in enumerate(protein_fnames):
        ppath = Path(pfname).resolve()
        xlps = read_dat_xlp(ppath)
        data_arr = [p.get_dat() for p in xlps]
        # with ppath.open(mode='r') as pf:
        #     lines = pf.readlines()
        #     data_arr = []
        #     for line in lines[2:]:
        #         data_arr += [line.split()[1:]]
        #     protein_dset[:, :, frame] = np.asarray(data_arr, dtype='f8')
        protein_dset[:, :, frame] = np.asarray(data_arr, dtype='f8')

    return protein_dset


# @profile
def read_sylinder_data(tubule_fnames, posit_grp):
    """!Read in data from all tubule files

    @param tubule_fnames: List of tubule posit file names
    @param posit_grp: HDF5 position data gropu
    @return: HDF5 data set containing tubule data

    """
    nframes = len(tubule_fnames)
    # Get number of tubules
    with open(tubule_fnames[0], 'r') as fn:
        ntubules = int(fn.readline())
    # Create dataset for MT info
    sy_dset = posit_grp.create_dataset('sylinders',
                                       shape=(ntubules, 9, nframes))
    sy_dset.attrs['ntubules'] = ntubules
    sy_dset.attrs['axis dimensions'] = ['sylinders', 'state', 'frame']
    # Set tubule attribute: names of columns
    sy_dset.attrs['column labels'] = ['gid', 'radius',
                                      'minus pos x', 'minus pos y', 'minus pos z',
                                      'plus pos x', 'plus pos y', 'plus pos z',
                                      'group', ]
    for frame, tfname in enumerate(tubule_fnames):
        fpath = Path(tfname).resolve()
        filaments = read_dat_sylinder(fpath)

        data_arr = [fil.get_dat()
                    for fil in filaments if (fil.fil_type is not 'L')]

        sy_dset[:, :, frame] = np.asarray(data_arr, dtype='f8')
    return sy_dset


def convert_dat_to_hdf(fname="TS_data.h5", path='.'):
    """!TODO: Docstring for convert_dat_to_hdf.

    @param fname: TODO
    @return: TODO

    """
    # TODO: Check to make sure path is a Path object <01-03-21, ARL> #
    data_dir = path / 'result'
    # Create h5 data set
    h5_data = h5py.File(fname, 'w')

    # Get list of tubule files, sort according to time
    sy_dat_paths = sorted(data_dir.glob("**/SylinderAscii*.dat"),
                          key=get_file_number)
    # Get list of all protein files, sort according to time
    xlp_dat_paths = sorted(data_dir.glob("**/ProteinAscii*.dat"),
                           key=get_file_number)
    # tubule_fnames = list(dpath.glob('SylinderAscii_*.dat'))
    # tubule_fnames.sort(key=lambda x: int(re.findall(r"[0-9]+", x.name)[0]))

    # assert(len(protein_fnames) == len(tubule_fnames))
    with open(path / 'RunConfig.yaml', 'r') as rc_file:
        rc_params = yaml.safe_load(rc_file)
        h5_data.attrs['RunConfig'] = yaml.dump(rc_params)

    with open(path / 'ProteinConfig.yaml', 'r') as xlp_file:
        xlp_params = yaml.safe_load(xlp_file)
        h5_data.attrs['ProteinConfig'] = yaml.dump(xlp_params)

    time_dset = read_time(sy_dat_paths, h5_data)
    # Create group of position data
    posit_grp = h5_data.create_group('raw_data')

    # protein_dset = read_protein_data(protein_fnames, posit_grp)

    sy_dset = read_sylinder_data(sy_dat_paths, posit_grp)
    xlp_dset = read_protein_data(xlp_dat_paths, posit_grp)

    return h5_data


##########################################
if __name__ == "__main__":
    print("Not implemented.")
