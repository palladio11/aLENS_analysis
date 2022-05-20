import re
import os
from pathlib import Path


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', str(filename)).group(1))


def result_to_pvd(path, name, ext):
    FileList = list(path.glob(f'./result*/{name}_*.{ext}'))
    # sort as numerical order
    FileList.sort(key=getFrameNumber_lambda)
    with open(path / f'{name}{ext}.pvd', "w") as fpvd:
        fpvd.write(
            '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64"> \n')
        fpvd.write('<Collection>\n')
        for i, f in enumerate(FileList):
            fpvd.write(f'<DataSet timestep="{i}"  file="{str(f)}"/>\n')
        fpvd.write('</Collection>\n')
        fpvd.write('</VTKFile> \n')


def make_pvd_files(result_dir):
    result_to_pvd(result_dir, 'Sylinder', 'pvtp')
    result_to_pvd(result_dir, 'Protein', 'pvtp')
    result_to_pvd(result_dir, 'ConBlock', 'pvtp')
