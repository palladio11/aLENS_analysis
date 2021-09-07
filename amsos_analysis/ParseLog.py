import numpy as np
import re

import Util.AMSOS as am

parser = am.getDefaultArgParser('Parse log file')
parser.add_argument('--logfile', type=str, dest='logfile',
                    help='logfile to be parsed')
parser.add_argument('--maxline', type=int, dest='maxline', default=0,
                    help='max number of lines to parse, default unlimited')

args = parser.parse_args()

print('logfile: ', args.logfile)
print('maxline: ', args.maxline)


# spdlong line example:
# [2021-09-06 16:10:51.194] [rank 0] [warning] CurrentStep 36255
# or
# [2021-09-06 16:12:00.532] [rank 0] [info] RECORD:

def parseOneStep(lines):
    '''parse lines starts from "CurrentTime xxx.xxx" '''
    timestamp = np.nan
    wtime_mainloop = np.nan
    wtime_bcqp = np.nan
    bcqp_steps = np.nan
    bcqp_residue = np.nan

    # preprocess, remove spdlog headers
    msg_levels = ['[trace]',
                  '[debug]',
                  '[info]',
                  '[warn]',
                  '[err]',
                  '[critical]']
    newlines = []
    for line in lines:
        for level in msg_levels:
            pos = line.find(level)
            if pos != -1:
                newlines.append(line[pos+len(level)+1:])

    for line in newlines:
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
    pos = line.find('CurrentTime')
    if pos != -1:
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
