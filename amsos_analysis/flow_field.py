# %%
import vtk
import numpy as np
import yaml
import re
from pathlib import Path
import time
import Util.AMSOS as am
from mpi4py import MPI
from amsos_analysis.PySTKFMM import Stk3DFMM, DArray, KERNEL
from amsos_analysis.PrintVTKData import Frame
# %%

# cubic mesh has (mesh+1) points in each direction
# example: mesh=5, box=1:
# grid points along each coord axis: [0,0.2,0.4,0.6,0.8,1.0]
# grid spacing = box/mesh
mesh = 128
write_vti = True


def read_params(configyaml):
    global boxL
    global box_origin
    global mult_order
    global viscosity
    # global dipole_offset
    # global dipole_strength

    # read params into global variables
    with open(configyaml, 'r') as cyf:
        config = yaml.load(cyf, Loader=yaml.FullLoader)
    # with open(cellyaml, 'r') as clf:
    #     cell = yaml.load(clf, Loader=yaml.FullLoader)

    boxLow = config['simBoxLow']
    boxHigh = config['simBoxHigh']
    boxX = boxHigh[0] - boxLow[0]  # um
    boxY = boxHigh[1] - boxLow[1]  # um
    boxZ = boxHigh[2] - boxLow[2]  # um
    assert np.abs(boxX - boxY) < 1e-6
    assert np.abs(boxX - boxZ) < 1e-6
    # box_origin = np.array([0, 0, 0])
    # box_origin = np.array([-3, -3, -3])
    box_origin = np.asarray(boxLow)
    # box_origin[0] = boxLow[0]
    # box_origin[1] = boxLow[1]
    # box_origin[2] = boxLow[2]
    # print(box_origin)
    viscosity = config['viscosity']
    boxL = boxX

    mult_order = 8
    # mult_order = cell['fmmMultOrder']
    # dipole_offset = cell['dipoleOffset']
    # dipole_strength = cell['dipoleStrength']
    # TODO Read in fluid force data file

    print('Params from yaml files:')
    print('mult_order: ', mult_order)
    print('boxL: ', boxL)
    print('viscosity: ', viscosity)
    # print('dipole_offset: ', dipole_offset)
    # print('dipole_strength: ', dipole_strength)

    return


def print_rank0(*args, **kwargs):
    # Convenience wrapper to print only on MPI rank 0
    if rank == 0:
        print(*args, **kwargs)


def create_cubicgrid(mesh, dx, origin=[0., 0., 0.]):
    # a grid matching vtkImageData order
    npts = (mesh)**3
    coord = np.zeros((npts, 3))
    index = 0
    for z in range(mesh):
        for y in range(mesh):
            for x in range(mesh):
                # x,y,z are integers
                coord[index, 0] = x * dx + origin[0]
                coord[index, 1] = y * dx + origin[1]
                coord[index, 2] = z * dx + origin[2]
                index += 1
    return coord


def write_cubicgrid(mesh, dx, data, filename):
    if write_vti:
        filename = str(filename) + '.vti'
        imageData = vtk.vtkImageData()
        imageData.SetDimensions(mesh, mesh, mesh)
        imageData.SetSpacing(dx, dx, dx)
        imageData.SetOrigin(box_origin[0], box_origin[1], box_origin[2])
        datadim = data.shape[1]
        imageData.AllocateScalars(vtk.VTK_DOUBLE, datadim)

        dims = imageData.GetDimensions()
        npts = imageData.GetNumberOfPoints()

        # Fill every entry of the image data with "2.0"
        index = 0
        for z in range(dims[2]):
            for y in range(dims[1]):
                for x in range(dims[0]):
                    for k in range(datadim):
                        imageData.SetScalarComponentFromDouble(
                            x, y, z, k, data[index, k] * (1 / viscosity))
                    index += 1

        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(imageData)
        writer.Write()
        writer.Write()
    else:
        data *= (1 / viscosity)
        np.save(filename, data)


def calc_flow_write(frame, output):
    # setup src_coord/src_value
    # columns: gid, radius, end0(dim3), end1(dim3)
    if rank:
        Sdata = DArray(None)
    else:
        syl_list = []
        for syl in frame.sylinders:
            center = .5 * (np.asarray(syl.end0) + np.asarray(syl.end1))
            syl_list += [center.tolist() + list(syl.force) + list(syl.radius)]

        Sdata = DArray(np.asarray(syl_list))

    # Sdata = DArray(None if rank else np.loadtxt(
    #     file, skiprows=2, usecols=(1, 2, 3, 4, 5, 6, 7, 8)))
    Sdata.scatter()
    # print("rank: ",rank, "Sdata local: ", Sdata.chunk.shape)

    nrod = Sdata.chunk.shape[0]
    src_coord = np.zeros((nrod, 3))
    src_value = np.zeros((nrod, 4))
    for i in range(nrod):
        src_coord[i, :] = Sdata.chunk[i, :3]
        src_value[i, :] = Sdata.chunk[i, 3:]
        # src_coord[2*i, :] = center+orient*dipole_offset
        # src_coord[2*i+1, :] = center-orient*dipole_offset
        # src_value[2*i, :3] = orient*dipole_strength
        # src_value[2*i+1, :3] = -orient*dipole_strength
        # src_value[2*i, 3] = radius
        # src_value[2*i+1, 3] = radius

    # call fmm
    fmm.set_box(box_origin, boxL)
    fmm.set_points(src_coord, trg_coord.chunk, np.empty((0, 3)))
    fmm.setup_tree(KERNEL.RPY)
    fmm.clear_fmm(KERNEL.RPY)

    before = time.time()
    trg_value.scatter()
    trg_value.chunk.fill(0)
    fmm.evaluate_fmm(KERNEL.RPY, src_value, trg_value.chunk, np.empty((0, 0)))
    trg_value.gather()
    comm.barrier()
    after = time.time()
    print_rank0("fmm time: ", after - before)

    if rank == 0:
        before = time.time()
        write_cubicgrid(mesh, dx, trg_value.data, output)
        after = time.time()
        print_rank0("writing time: ", after - before)
    print(trg_value.data)


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', str(filename)).group(1))


# %%
pbc = 0  # Number of dimensions to periodize (0, 1, 2, 3)
# wrk_dir = Path(
# '/home/alamson/DATA/Chromatin/21-08-26_aLchr1_LEF_acute_angle_hydro_test/result')
wrk_dir = Path.cwd()
read_params(wrk_dir / '../RunConfig.yaml')

# initialize FMM
dx = boxL * 1.0 / mesh
max_pts = 2000  # Max points per OctTree leaf
# Get MPI parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print(box_origin)


fmm = Stk3DFMM(mult_order, max_pts, pbc, KERNEL.RPY)

kdimSL = 4
kdimDL = 0
kdimTrg = 6
# initialize trg coord and values
trg_coord = DArray(None if rank else create_cubicgrid(mesh, dx, box_origin))
print(trg_coord.data)
trg_coord.scatter()
trg_value = DArray(None if rank else np.zeros(((mesh)**3, kdimTrg)))
trg_value.scatter()

print(trg_coord.chunk.shape)
print(trg_value.chunk.shape)


syl_vtk_files = list(wrk_dir.glob('./result*-*/Sylinder_*.pvtp'))
prot_vtk_files = list(wrk_dir.glob('./result*-*/Protein_*.pvtp'))
con_vtk_files = list(wrk_dir.glob('./result*-*/ConBlock_*.pvtp'))
syl_vtk_files.sort(key=getFrameNumber_lambda)
prot_vtk_files.sort(key=getFrameNumber_lambda)
con_vtk_files.sort(key=getFrameNumber_lambda)

folder = wrk_dir / f'flow{mesh}'
try:
    folder.mkdir()
except FileExistsError:
    pass
# SylinderFileList = am.getFileListSorted('./result*-*/SylinderAscii_*.dat')

index = 0
for syl_file, prot_file, con_file in zip(
        syl_vtk_files, prot_vtk_files, con_vtk_files):
    frame = Frame(str(syl_file), str(prot_file),
                  str(con_file), print_flag=False)
    outfile = folder / f"flow_{getFrameNumber_lambda(syl_file)}"
    calc_flow_write(frame, outfile)
    index += 1
#     comm.barrier()

# %%
