import os
import numba as nb
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import structured_to_unstructured
import h5py
import glob
import pyvista as pv
import vtk
import Util.AMSOS as am

h5name = 'Trajectory'
foldername = 'Traj2PolyLines'

# one poly line file is written for each gid.
# avoid too many small files
# a few hundred is ok

rng = np.random.default_rng()
gids = rng.integers(low=0, high=100000, size=100)


def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    return poly


def traj2Polyline(h5name, gid_choice=[]):
    '''convert hdf5 trajectory to polylines'''
    file = h5py.File(h5name+'.hdf5', 'r')
    steps = list(file.keys())
    nTraj = np.array(file[steps[0]]['traj'].shape[0])
    nSteps = len(steps)
    # access data as this np.array(file[keys[0]]['traj'])
    if len(gid_choice) <= 0:
        gid_choice = range(nTraj)
    print('gid_choice: ', gid_choice)

    for gid in gid_choice:
        if gid < 0 or gid >= nTraj:
            print('invalid gid: ', gid)
            continue
        pts = np.zeros([nSteps, 3])
        for k in range(nSteps):
            pts[k, :] = np.array(file[steps[k]]['traj'])[gid, :]
        pl = polyline_from_points(pts)
        pl.save(foldername+'/TrajVTK_{:08d}.vtp'.format(gid), binary=True)

    file.close()
    return


def mergePolyline():
    '''merge all vtk polylines into a single file'''
    reader = vtk.vtkXMLPolyDataReader()
    append = vtk.vtkAppendPolyData()
    filenames = glob.glob(foldername+'/TrajVTK_*.vtp')
    for file in filenames:
        reader.SetFileName(file)
        reader.Update()
        polydata = vtk.vtkPolyData()
        polydata.ShallowCopy(reader.GetOutput())
        append.AddInputData(polydata)
    append.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('TrajAll.vtp')
    writer.SetInputData(append.GetOutput())
    writer.Write()


am.mkdir(foldername)
# traj2Polyline(h5name, gids)
mergePolyline()
