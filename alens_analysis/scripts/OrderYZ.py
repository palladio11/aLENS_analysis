import yaml
import sys
import vtk
import glob
import re
import os
import argparse

import numba as nb
import numpy as np
import scipy as sp
import scipy.spatial as ss
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py

parser = argparse.ArgumentParser()
# parser.add_argument('pbcX', type=float, help='periodic bc length along X')
# parser.add_argument('pbcY', type=float, help='periodic bc length along Y')
# parser.add_argument('pbcZ', type=float, help='periodic bc length along Z')
parser.add_argument('ngrid', type=int, help='number of samples along X axis')
parser.add_argument('--rcut', type=float,
                    help='cut-off radius of g(r), default 0.1um', default=0.1)

args = parser.parse_args()

# read box from RunConfig.yaml
config_file = open('../RunConfig.yaml', 'r')
config = yaml.load(config_file, Loader=yaml.FullLoader)
boxLow = config['simBoxLow']
boxHigh = config['simBoxHigh']
pbcX = boxHigh[0]-boxLow[0]  # um
pbcY = boxHigh[1]-boxLow[1]  # um
pbcZ = boxHigh[2]-boxLow[2]  # um
config_file.close()

gridX = np.linspace(0, pbcX, num=args.ngrid)
density_max = 0.9069/(0.25*np.pi*(0.025**2))  # hex packing number density
rcut = args.rcut

print('pbc X:', pbcX)
print('pbc Y:', pbcY)
print('pbc Z:', pbcZ)
print('ngrid:', args.ngrid)
print('rcut :', args.rcut)


def find_intersection(end0, end1, xloc, xpbc):
    p0 = end0
    p1 = end1
    assert xpbc > 0
    if (end0[0]-xloc)*(end1[0]-xloc) < 0:  # intersect
        pass
    elif (end0[0]-xpbc-xloc)*(end1[0]-xpbc-xloc) < 0:  # intersect with image -xpbc
        p0[0] -= xpbc
        p1[0] -= xpbc
    elif (end0[0]+xpbc-xloc)*(end1[0]+xpbc-xloc) < 0:  # intersect with image +xpbc
        p0[0] += xpbc
        p1[0] += xpbc
    else:
        return None
    vec = p1-p0
    frac = (xloc-p0[0])/vec[0]
    point = p0+frac*vec
    return point


# find closest point
def closest_point(target, points):
    target = np.array(target)
    points = np.array(points)
    distance = []
    for p in points:
        distance.append(np.linalg.norm(p-target))
    distance = np.array(distance)
    ind = np.argmin(distance)
    return points[ind], ind


def closest_point1d(target, points):
    distance = []
    for p in points:
        distance.append(np.abs(p-target))
    distance = np.array(distance)
    ind = np.argmin(distance)
    return points[ind], ind


def get_closetimage(target, source, boxsize):
    dim = target.shape[0]
    assert source.shape[0] == dim
    image = np.zeros(dim)
    for i in range(dim):
        pos, ind = closest_point1d(
            target[i], [source[i], source[i]-boxsize[i], source[i]+boxsize[i]])
        image[i] = pos
    return image


def applypbc(points, boxsize):
    # reset coord to [0,boxsize)
    dim = len(boxsize)
    for d in range(dim):
        x = points[:, d]
        boxL = boxsize[d]
        newx = x-np.floor(x/boxL)*boxL
        points[:, d] = newx
    return points


def get_pair(coords, boxsize, rcut=None):
    tree = ss.cKDTree(data=coords, boxsize=boxsize)
    boxsize = np.array(boxsize)
    if rcut == None:
        rcut = np.sqrt(boxsize.dot(boxsize))/10
    pairs = tree.query_pairs(r=rcut)  # this returns only pairs (i<j)
    pairs2 = set()
    for p in pairs:
        pairs2.add((p[1], p[0]))
    for p in pairs2:
        pairs.add(p)
    return pairs


def get_rvec(coords, boxsize, rcut, pairs):
    rvec = []
#     print(pairs)
    for pair in pairs:
        id0 = pair[0]
        id1 = pair[1]
        pos0 = coords[id0]
        pos1 = coords[id1]
        vec01 = pos1-pos0
        if np.linalg.norm(vec01) < rcut:
            rvec.append(vec01)
        else:  # fix periodic image
            image = get_closetimage(pos0, pos1, boxsize)
            rvec.append(image-pos0)
    return np.array(rvec)


def gen_rdf(rvec, npar, density, rcut=0, nbins=20):
    rnorm = np.linalg.norm(rvec, axis=1)
    lb = 0
    if rcut == 0:
        ub = np.max(rnorm)
    else:
        ub = rcut
    dim = rvec.shape[1]
    bins = np.linspace(lb, ub, nbins)
    count, bins = np.histogram(rnorm, bins)
    print(count)
    # print(bins)
    # scale with vol and density
    vol = np.zeros(count.shape)
    if dim == 2:    # area = pi(r1^2-r0^2)
        for i in range(nbins-1):
            vol[i] = np.pi*(bins[i+1]**2-bins[i]**2)
    elif dim == 3:  # area = 4pi/3(r1^3-r0^3)
        for i in range(nbins-1):
            vol[i] = (4.0/3.0)*np.pi*(bins[i+1]**3-bins[i]**3)
    rdf = count/(npar*vol*density)
    r = 0.5*(bins[:-1]+bins[1:])
    return r, rdf


def gen_rdf2d(rvec, npar, density, rcut=0, nbins=20):
    if rcut == 0:
        lbx = np.min(rvec[:, 0])
        ubx = np.max(rvec[:, 0])
        lby = np.min(rvec[:, 1])
        uby = np.max(rvec[:, 1])
    else:
        lbx = -rcut*0.7
        ubx = rcut*0.7
        lby = -rcut*0.7
        uby = rcut*0.7
    dim = rvec.shape[1]
    assert dim == 2
    xedges = np.linspace(lbx, ubx, nbins)
    yedges = np.linspace(lby, uby, nbins)
    H, xedges, yedges = np.histogram2d(
        rvec[:, 0], rvec[:, 1], bins=(xedges, yedges))
    H = H/(npar*density*(xedges[1]-xedges[0])*(yedges[1]-yedges[0]))
    return H, xedges, yedges


def plot_rdf(r, rdf, gmax=0):
    fig, ax0 = plt.subplots(nrows=1)
    im = ax0.plot(r, rdf)
    ax0.set_xlim(0, rcut)
    if gmax > 0:
        ax0.set_ylim(0, gmax)
    ax0.set_xlabel('r/um')
    ax0.set_ylabel('$g(r)/n_{cp}$')
    # ax0.set_title('')
    fig.tight_layout()
    return fig


def plot_rdf2d(H, xedges, yedges, gmax=0):
    X, Y = np.meshgrid(xedges, yedges)
    fig, ax0 = plt.subplots(nrows=1)
    if gmax > 0:
        im = ax0.pcolormesh(X, Y, H, cmap='jet', vmin=0, vmax=gmax)
    else:
        im = ax0.pcolormesh(X, Y, H, cmap='jet', vmin=0)
    ax0.set_xlabel('r/um')
    ax0.set_ylabel('r/um')
    fig.colorbar(im, ax=ax0)
    # ax0.set_title('')
    fig.tight_layout()
    return fig


def plot_rdf_combine(r, g, H, x, y, gmax=0):
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(7, 3))
    rmax = np.max(r)
    # g(r), convert to nm
    ax0.plot(1000*r, g)
    if gmax > 0:
        ax0.set_ylim(0, gmax)
    else:
        ax0.set_ylim(0, np.max(g)*1.2)
    ax0.set_xlim(0, rmax*1000)
    ax0.set_xlabel('r/nm')
    ax0.set_ylabel('$g(r)/n_{cp}$')
    # g2d
    X, Y = np.meshgrid(x, y)
    # convert to nm
    X *= 1000
    Y *= 1000
    if gmax > 0:
        im = ax1.pcolormesh(X, Y, H, cmap='jet', vmin=0, vmax=gmax)
    else:
        im = ax1.pcolormesh(X, Y, H, cmap='jet', vmin=0)
    fig.colorbar(im, ax=ax1, label='$g(r)/n_{cp}$')
    ax1.set_aspect(aspect=1)
    ax1.set_xlabel('r/nm')
    ax1.set_ylabel('r/nm')
    for axis in [ax1.xaxis, ax1.yaxis]:
        axis.set_major_locator(mpl.ticker.MaxNLocator(7))
#     ax1.set_title('pcolormesh with levels')
    fig.tight_layout()
    return fig


class Sylinder(object):
    end0 = None
    end1 = None
    pass


class Frame:

    def __init__(self, sylinderFile=None):
        self.sylinders = []
        self.parseSylinderFile(sylinderFile)

    def parseFile(self, dataFile, objType, objList):
        print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        print("parsing data for ", nObj, " sylinders")
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

        print("-------------------------------------")
        self.sylinders.sort(key=lambda x: x.gid, reverse=False)

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)

    def printData(self):
        # output all data for debug
        for s in self.sylinders[:10]:
            # print(s.end0, s.end1)
            attrs = vars(s)
            print('*************************************')
            print('\n'.join("%s: %s" % item for item in attrs.items()))
            print('*************************************')


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def genEndHistory():
    SylinderFileList = glob.glob('./result*/Sylinder_*.pvtp')
    # sort as numerical order
    SylinderFileList.sort(key=getFrameNumber_lambda)
    print(SylinderFileList)

    endHistory = []
    frame = Frame(SylinderFileList[0])
    for i in range(1, len(SylinderFileList[:])):
        # for i in range(1, 10):
        pos = np.zeros((len(frame.sylinders), 7))
        # get frame
        frame = Frame(SylinderFileList[i])
        for j in range(len(frame.sylinders)):
            sy = frame.sylinders[j]
            pos[j] = np.array(
                [sy.gid[0],
                 sy.end0[0], sy.end0[1], sy.end0[2],
                 sy.end1[0], sy.end1[1], sy.end1[2]])
        endHistory.append(pos)

    # generate file
    endHistory = np.array(endHistory)
    endHistory = np.swapaxes(endHistory, 0, 1)
    np.save('endHistory', endHistory)
    pass


def genPairVec():
    data = np.load('endHistory.npy')
    nsteps = data.shape[1]
    pairVec = h5py.File('pairVec.hdf5', 'w')
    pairVec.attrs["pbcX"] = pbcX
    pairVec.attrs["pbcY"] = pbcY
    pairVec.attrs["pbcZ"] = pbcZ
    pairVec.attrs["rcut"] = rcut
    pairVec.attrs["density_max"] = density_max
    pairVec.attrs["ngridX"] = len(gridX)

    # hdf5 file hierarchy:
    # step_i/xloc_xindex/points
    # step_i/xloc_xindex/rvec
    for i in range(nsteps):
        rods = data[:, i, :]

        grp = pairVec.create_group("step_"+str(i))
        xindex = 0
        for xloc in gridX:
            points = []
            for r in rods:
                end0 = r[1:4]
                end1 = r[4:7]
                point = find_intersection(end0, end1, xloc, pbcX)
                if point is not None:
                    points.append(point[1:])  # Y-Z coordinate
            boxsize = [pbcY, pbcZ]
            points = np.array(points)
            applypbc(points, boxsize)
            pairs = get_pair(np.array(points), boxsize, rcut)
            rvec = get_rvec(points, boxsize, rcut, pairs)
            grp_step = grp.create_group('xloc_'+str(xindex))
            dset1 = grp_step.create_dataset(
                "points", points.shape, dtype='float64', compression="lzf", chunks=True)
            dset1[...] = points[...]
            dset2 = grp_step.create_dataset(
                "rvec", rvec.shape, dtype='float64', compression="lzf", chunks=True)
            dset2[...] = rvec[...]
            xindex += 1
    pairVec.close()


def genOrderYZ():
    data = h5py.File('pairVec.hdf5', 'r')
    nsteps = len(data)
    print(nsteps)

    folderName = 'order_yz'
    try:
        os.mkdir(folderName)
    except FileExistsError:
        pass

    for step in range(nsteps):
        npartotal = 0
        rvectotal = np.array([], dtype=np.float).reshape(0, 2)
        frame = data['step_'+str(step)]
        # each step has ngrid crosssection samples
        assert len(frame) == args.ngrid
        for xloc in frame:
            points = frame[xloc]['points']
            rvec = frame[xloc]['rvec']
            npartotal += len(points)
            rvectotal = np.concatenate((rvectotal, rvec), axis=0)
        print(npartotal)
        print(rvectotal.shape)
        # generate rdf and rdf_2d
        meandensity = density_max
        r, rdf = gen_rdf(rvectotal, npartotal, meandensity, rcut, 200)
        H, X, Y = gen_rdf2d(rvectotal, npartotal, meandensity, rcut, 50)
        fig = plot_rdf_combine(r, rdf, H, X, Y, gmax=1.2)
        fig.savefig(folderName+os.sep +
                    'rdf2d_{:06d}'.format(step)+'.png', dpi=300)
        plt.close(fig)
    # generate movie
    path = os.path.abspath(os.path.curdir)
    names = path.split('/')
    casename = names[-2]
    cmd = "ffmpeg -y -framerate 30 -pattern_type glob -i '"+folderName+os.sep+'rdf2d_' + \
        "*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 " + \
        casename + '_rdf2d.mp4'
    print(cmd)
    os.system(cmd)

    return


def main():
    if not os.path.isfile('endHistory.npy'):
        genEndHistory()
    if not os.path.isfile('pairVec.hdf5'):
        genPairVec()
    genOrderYZ()


if __name__ == '__main__':
    main()
