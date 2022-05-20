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
from codetiming import Timer

h5name = 'Trajectory'
foldername = 'Traj2PolyLines'

# one poly line file is written for each gid.
# avoid too many small files
# a few hundred is ok

parser = am.getDefaultArgParser(
    'Calculate net displacement along certain direction')
parser.add_argument('--ntraj', type=int, dest='ntraj', default=100,
                    help='number of trajs')
parser.add_argument('--start', type=int, dest='start', default=0,
                    help='start frame number of traj')
parser.add_argument('--end', type=int, dest='end', default=-1,
                    help='end frame number of traj')

args = parser.parse_args()

rng = np.random.default_rng(seed=0)
gids = rng.integers(low=0, high=100000, size=args.ntraj)


def polyline_from_points(points):
    poly = pv.PolyData()
    poly.points = points
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))
    poly.lines = the_cell
    poly["step"] = np.arange(poly.n_points)
    return poly


@Timer("One Polyline")
def traj2Polyline(h5name, gid_choice=[], start=0, end=-1):
    '''convert hdf5 trajectory to polylines'''
    file = h5py.File(h5name+'.hdf5', 'r')
    steps = list(file.keys())
    nTraj = np.array(file[steps[0]]['traj'].shape[0])
    nSteps = len(steps)
    # access data as this np.array(file[keys[0]]['traj'])
    if len(gid_choice) <= 0:
        gid_choice = range(nTraj)
    print('gid_choice: ', gid_choice)

    end = nSteps+end+1 if end < 0 else end
    pts = np.zeros((len(gid_choice), end-start, 3))

    for k in range(start, end):
        pts[:, k-start,
            :] = np.array(file[steps[k]]['traj'])[gid_choice, :]

    for i in range(len(gid_choice)):
        gid = gid_choice[i]
        pl = polyline_from_points(pts[i])
        pl.save(
            foldername+'/TrajVTK_{:08d}_{:d}_{:d}.vtp'.format(gid, start, end), binary=True)

    file.close()
    return


def mergePolyline():
    '''merge all vtk polylines into a single file'''
    reader = vtk.vtkXMLPolyDataReader()
    append = vtk.vtkAppendPolyData()
    filenames = glob.glob(
        foldername+'/TrajVTK_*_{:d}_{:d}.vtp'.format(args.start, args.end))
    for file in filenames:
        reader.SetFileName(file)
        reader.Update()
        polydata = vtk.vtkPolyData()
        polydata.ShallowCopy(reader.GetOutput())
        append.AddInputData(polydata)
    append.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('TrajAll_{:d}_{:d}.vtp'.format(args.start, args.end))
    writer.SetInputData(append.GetOutput())
    writer.Write()


am.mkdir(foldername)
traj2Polyline(h5name, gids, args.start, args.end)
mergePolyline()
