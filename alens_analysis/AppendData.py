import Util.aLENS as am
import h5py as h5
import numpy as np
import pyvista as pv

# open data
foldername = 'Sylinder_SND'

h5FileName = 'OrderLinkS.hdf5'
orderFile = h5.File(h5FileName, 'r')


# append data and export new vtk file
def process_frame(syFile):
    frame = am.FrameVTK(syFile)
    points = frame.data["points"]
    GID = frame.data["gid"]
    basename = am.get_basename(syFile)
    asciiname = basename.replace('Sylinder', 'SylinderAscii')

    # ordered by increasing gid
    order_data = orderFile[asciiname]['OrderLinkS'][:]
    num_sy = len(GID)
    S = np.zeros(num_sy)
    ND = np.zeros(num_sy)

    for i in range(num_sy):
        gid = GID[i]
        S[i] = order_data[gid, 0]
        ND[i] = order_data[gid, 1]

    pointsPerLine = np.zeros(num_sy, dtype=int)
    pointsPerLine.fill(2)
    x = points[:, 0].copy(order='C')
    y = points[:, 1].copy(order='C')
    z = points[:, 2].copy(order='C')
    polyLinesToVTK(path=foldername+'/'+basename, x=x, y=y, z=z,
                   pointsPerLine=pointsPerLine,
                   cellData={"gid": GID, "S": S, "ND": ND})


# open vtk file
SylinderFileList = am.getFileListSorted("./result*-*/Sylinder_*.pvtp")
for f in SylinderFileList:
    process_frame(f)
