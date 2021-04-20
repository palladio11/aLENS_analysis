import os
import scipy.special as ss
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re
import os
import sys

currentdir = os.path.dirname(os.path.realpath(__file__))  # nopep8
parentdir = os.path.dirname(currentdir)  # nopep8
sys.path.append(parentdir)  # nopep8

import Util.AMSOS as am
import Util.HDF5_Wrapper as h5
import point_cloud.PointCloud as pc

parser = am.getDefaultArgParser('Parse log file')
parser.add_argument('--logfile', type=str, dest='logfile',
                    help='logfile to be parsed')
parser.add_argument('--maxline', type=int, dest='maxline', default=0,
                    help='max number of lines to parse, default unlimited')

args = parser.parse_args()

print('logfile: ', args.logfile)
print('maxline: ', args.maxline)


def parseOneStep(lines):
    '''parse lines starts from "CurrentTime xxx.xxx" '''
    timestamp = np.nan
    wtime_mainloop = np.nan
    wtime_bcqp = np.nan
    bcqp_steps = np.nan
    bcqp_residue = np.nan

    for line in lines:
        if line.startswith('CurrentTime'):
            timestamp = float(line.split(' ')[1])
        if line.startswith('RECORD: BCQP residue'):
            data = line.split(',')
            bcqp_steps = float(data[-1])
            bcqp_residue = float(data[-2])
        if line.startswith('AMSOS main loop'):
            data = re.split(' +', line)
            wtime_mainloop = float(data[5])
        if line.startswith('SylinderSystem::SolveConstraints'):
            data = re.split(' +', line)
            wtime_bcqp = float(data[3])

    return [timestamp, wtime_bcqp, wtime_mainloop, bcqp_steps, bcqp_residue]


file = open(args.logfile, 'r')
msg_onestep = []
data = []
count = 0
append = False  # skip until first CurrentTime tag
for line in file:
    if line.startswith('CurrentTime'):
        if append:
            data.append(parseOneStep(msg_onestep))
        msg_onestep.clear()
        append = True
    if append:
        msg_onestep.append(line)
        count += 1
    if count > args.maxline and args.maxline > 0:
        break
    pass

np.savetxt('logdata.txt', np.array(data), fmt='%.14g',
           header='timestamp, wtime_bcqp, wtime_mainloop, bcqp_steps, bcqp_residue')

file.close()
