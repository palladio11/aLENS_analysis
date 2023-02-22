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
import vtk
from vtk.util import numpy_support as vn
import time
import numpy as np
from pathlib import Path
from tqdm import tqdm
from .objects import filament, protein, con_block
from .runlog_funcs import get_walltime
import zipfile


def get_file_number(path):
    name = Path(path).stem
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
        proteins = sorted([protein(line)
                           for line in filecontent], key=lambda x: int(x.gid))
    return proteins


def read_dat_constraint(fpath):

    con_blocks = []
    # open the file and read the lines
    # with fpath.open('r') as file1:
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName(str(fpath))
    reader.Update()
    data = reader.GetOutput()
    # filecontent = file1.readlines()

    nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
    # print("parsing data for ", nObj, " object(s)")
    for i in range(nObj):
        cb = con_block()
        cb.end0 = data.GetPoints().GetPoint(2 * i)
        cb.end1 = data.GetPoints().GetPoint(2 * i + 1)
        con_blocks += [cb]

    # step 2, member cell data
    numCellData = data.GetCellData().GetNumberOfArrays()
    # print("Number of CellDataArrays: ", numCellData)
    for i in range(numCellData):
        cdata = data.GetCellData().GetArray(i)
        dataName = cdata.GetName()
        # print("Parsing Cell Data", dataName)
        for j in range(len(con_blocks)):
            setattr(con_blocks[j], dataName, cdata.GetTuple(j))

    # step 3, member point data
    numPointData = data.GetPointData().GetNumberOfArrays()
    # print("Number of PointDataArrays: ", numPointData)
    for i in range(numPointData):
        pdata = data.GetPointData().GetArray(i)
        dataName = pdata.GetName()
        # print("Parsing Point Data", dataName)
        for j in range(len(con_blocks)):
            setattr(con_blocks[j], dataName + "0", pdata.GetTuple(2 * j))
            setattr(con_blocks[j], dataName + "1", pdata.GetTuple(2 * j + 1))

    return con_blocks


def read_stress_from_con(fpath):
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName(str(fpath))
    reader.Update()
    data = reader.GetOutput()

    stress_arr = vn.vtk_to_numpy(data.GetCellData().GetVectors('Stress'))
    bilat_flag_arr = vn.vtk_to_numpy(
        data.GetCellData().GetVectors('bilateral'))
    collision_stress = stress_arr[bilat_flag_arr ==
                                  0, :].sum(axis=0).reshape((3, 3))
    bilat_stress = stress_arr[bilat_flag_arr ==
                              1, :].sum(axis=0).reshape((3, 3))
    return bilat_stress, collision_stress

    # return stress_arr.sum(axis=0).reshape((3, 3))


def read_time(fpaths, h5_data):
    """!Read in data from all protein files

    @param fnames: List posit file names
    @param h5_data: HDF5 position data gropu
    @return: HDF5 data set containing protein data

    """
    t = []
    for fp in tqdm(fpaths, disable=True):
        # p = Path(fn).resolve()
        with fp.open(mode='r') as f:
            t += [float(f.readlines()[1])]
    time_dset = h5_data.create_dataset('time', data=t)
    return time_dset


# @profile
def read_sylinder_data(syl_paths, posit_grp):
    """!Read in data from all tubule files

    @param tubule_fnames: List of tubule posit file names
    @param posit_grp: HDF5 position data gropu
    @return: HDF5 data set containing tubule data

    """
    nframes = len(syl_paths)
    # Get number of tubules
    with syl_paths[0].open('r') as sp:
        n_syl = int(sp.readline())
    # Create dataset for MT info
    sy_dset = posit_grp.create_dataset('sylinders',
                                       shape=(n_syl, 9, nframes))
    sy_dset.attrs['n_yslinders'] = n_syl
    sy_dset.attrs['axis dimensions'] = ['sylinders', 'state', 'frame']
    # Set tubule attribute: names of columns
    sy_dset.attrs['column labels'] = ['gid', 'radius',
                                      'minus pos x', 'minus pos y', 'minus pos z',
                                      'plus pos x', 'plus pos y', 'plus pos z',
                                      'group', ]
    for frame, syl_path in tqdm(enumerate(syl_paths), total=len(syl_paths), disable=True):
        filaments = read_dat_sylinder(syl_path)
        data_arr = [fil.get_dat()
                    for fil in filaments if (fil.fil_type != 'L')]
        sy_dset[:, :, frame] = np.asarray(data_arr, dtype='f8')
    return sy_dset


def read_protein_data(xlp_paths, posit_grp):
    """!Read in data from all protein files

    @param protein_fnames: List of protein posit file names
    @param posit_grp: HDF5 position data group
    @return: HDF5 data set containing protein data

    """
    nframes = len(xlp_paths)
    with xlp_paths[0].open(mode='r') as xp:
        nproteins = int(xp.readline())

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
    for frame, xlp_path in tqdm(enumerate(xlp_paths), total=len(xlp_paths), disable=True):
        # for frame, xlp_path in enumerate(xlp_paths):
        xlps = read_dat_xlp(xlp_path)
        data_arr = [p.get_dat() for p in xlps]
        protein_dset[:, :, frame] = np.asarray(data_arr, dtype='f8')

    return protein_dset


def read_constraint_data(cons_fnames, h5_data):
    """!Read in data from constraint files

    @param cons_fnames: List constraint file names
    @param h5_data: HDF5 data file to add stress
    @return: HDF5 data set containing tubule data

    """
    nframes = len(cons_fnames)
    bi_dset = h5_data.create_dataset('bilateral_stress',
                                     shape=(3, 3, nframes))
    col_dset = h5_data.create_dataset('collision_stress',
                                      shape=(3, 3, nframes))
    bi_dset.attrs['axis labels'] = ['dim', 'dim', 'frame']
    col_dset.attrs['axis labels'] = ['dim', 'dim', 'frame']
    for frame, tfname in enumerate(cons_fnames):
        collision_stress, bilateral_stress = read_stress_from_con(tfname)
        # print(f'Step {frame}: {collision_stress.flatten()}')
        col_dset[:, :, frame] = collision_stress
        bi_dset[:, :, frame] = bilateral_stress

    return bi_dset, col_dset


def collect_stress_from_con_pvtp(fname="stress.h5", path=Path('.')):
    result_dir = path / 'result'
    if not result_dir.exists():
        raise FileNotFoundError(
            f'Result directory {str(result_dir)} does not exist.')

    with h5py.File(fname, 'w') as h5_data:
        con_dat_paths = sorted(result_dir.glob("**/ConBlock*.pvtp"),
                               key=get_file_number)
        bi_dset, col_dset = read_constraint_data(con_dat_paths, h5_data)
        t4 = time.time()


def convert_dat_to_hdf(fname="raw_data.h5", path=Path('.'), store_stress=False):
    """Convert separate ascii and vtk data files into a single hdf5 file

    Parameters
    ----------
    fname : str, optional
        Name of the file to save, by default "raw_data.h5"
    path : Path object, optional
        The seed directory of the simulation, by default Path('.')
    store_stress : bool, optional
        Should you spend the space to store the stress calculated in the system, by default False

    Raises
    ------
    OSError
        If neither result.zip file or result directory do not exist, raise error.
    """
    if (path / 'result').exists():
        result_dir = path / 'result'
        is_zip = False
    elif (path / 'result.zip').exists():
        result_zip = zipfile.ZipFile(path / 'result.zip')
        result_path = zipfile.Path(path / 'result.zip')
        is_zip = True
    else:
        raise OSError(f'Could not find result directory or zipfile in {path}.')

    # Open raw h5 data objec to write to
    with h5py.File(fname, 'w') as h5_data:
        # Get paths (depends on if you are using zip archive or not)
        if is_zip:
            sy_reg = re.compile(r'.*SylinderAscii.*.dat')
            sy_dat_paths = sorted(list(filter(sy_reg.search, result_zip.namelist())),
                                  key=get_file_number)
            sy_dat_paths = [result_path / f for f in sy_dat_paths]

            xlp_reg = re.compile(r'.*ProteinAscii.*.dat')
            xlp_dat_paths = sorted(list(filter(xlp_reg.search, result_zip.namelist())),
                                   key=get_file_number)
            xlp_dat_paths = [result_path / f for f in xlp_dat_paths]
        else:
            sy_dat_paths = sorted(result_dir.glob("**/SylinderAscii*.dat"),
                                  key=get_file_number)
            xlp_dat_paths = sorted(result_dir.glob("**/ProteinAscii*.dat"),
                                   key=get_file_number)

        # assert(len(protein_fnames) == len(tubule_fnames))
        with (path / 'RunConfig.yaml').open('r') as rc_file:
            rc_params = yaml.safe_load(rc_file)
            h5_data.attrs['RunConfig'] = yaml.dump(rc_params)

        with (path / 'ProteinConfig.yaml').open('r') as xlp_file:
            xlp_params = yaml.safe_load(xlp_file)
            h5_data.attrs['ProteinConfig'] = yaml.dump(xlp_params)

        # Make time array
        t0 = time.time()
        time_dset = read_time(sy_dat_paths, h5_data)
        t1 = time.time()
        print(f"Made time data set in {t1-t0} seconds.")

        # Create group of position data
        posit_grp = h5_data.create_group('raw_data')

        # Make sylinder data
        sy_dset = read_sylinder_data(sy_dat_paths, posit_grp)
        t2 = time.time()
        print(f"Made sylinder data set in {t2-t1} seconds.")

        # Make protein data
        xlp_dset = read_protein_data(xlp_dat_paths, posit_grp)
        t3 = time.time()
        print(f"Made protin data set in {t3-t2} seconds.")

        # Make stress data
        if not is_zip:
            # Get list of all constraint files, sort according to time
            try:
                h5_stress_path = fname.parent / f'stress_{path.stem}.h5'
                con_dat_paths = sorted(result_dir.glob("**/ConBlock*.pvtp"),
                                       key=get_file_number)
                bi_dset, col_dset = collect_stress_from_con_pvtp(
                    h5_stress_path, path)
                t4 = time.time()
                print(f"Made stress data set in {t4-t3} seconds.")
            except:
                print("Could not make stress data.")

        # Wall time analysis
        # Check to see if run.log is present in current simulation
        log_path = list(path.glob('*run*.(log|out)'))
        if log_path:
            dwtime = get_walltime(log_path[0])  # TODO make this more robust
            h5_data.attrs['total_seconds'] = dwtime.total_seconds()
            h5_data.attrs['walltime'] = str(dwtime)
        t = time.time()
        print(f"Made raw data file in a total of {t-t0} seconds.")


##########################################
if __name__ == "__main__":
    print("Not implemented.")
