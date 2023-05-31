import h5py
import numpy as np


def newFile(name):
    '''create a file, clear if existing'''
    file = h5py.File(name+'.hdf5', 'w')
    file.close()
    return


def saveData(fname, data, path, dname, dtype):
    h5file = h5py.File(fname+'.hdf5', 'a')
    grp = h5file.require_group(path)
    dset = grp.create_dataset(dname, data.shape, dtype, chunks=True)
    dset[...] = data[...]
    h5file.close()

    return


if __name__ == '__main__':
    newFile('test')
    data = np.random.uniform(low=0, high=1, size=(100, 3))
    saveData('test', data, '/', 'u01', dtype=float)
    data = np.random.normal(0, 1, size=(100, 3))
    saveData('test', data, 'normal/', 'n01', dtype=float)
