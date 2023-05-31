import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Process XF from log file.')
parser.add_argument('-w', type=int, dest='window', default=2000,
                    help='moving average mean window')
parser.add_argument('file', type=str,
                    help='log file name')

args = parser.parse_args()
meanWindow = args.window
file = args.file

print('mean window: ', meanWindow)
print('log file: ', file)


def analyzeXF(file, meanWindow, record):
    grepString = record+'XF'
    recordFile = 'record_'+grepString+'.csv'
    cmd = 'grep '+grepString+' '+file+' > '+recordFile
    print(cmd)
    os.system(cmd)
    XFHistory = np.genfromtxt(recordFile, usecols=(
        1, 2, 3, 4, 5, 6, 7, 8, 9), delimiter=',')
    np.savetxt(recordFile, XFHistory, delimiter=',', fmt='%6g',
               header='Pxx,Pxy,Pxz,Pyx,Pyy,Pyz,Pzx,Pzy,Pzz')
    # tensor stat
    SigmaHistory = np.reshape(XFHistory, (-1, 3, 3))
    print("stress / nkbT:")
    print(np.mean(SigmaHistory[-meanWindow:], axis=0))
    # 3 eigen values for each tensor
    EigSigma = np.linalg.eigvals(np.mean(SigmaHistory[-meanWindow:], axis=0))
    print("stress eigen / nkbT: ", EigSigma)
    print("stress eigen max/min ratio: ", np.max(EigSigma)/np.min(EigSigma))
    # Trace stat
    PHistory = np.mean(XFHistory[:, [0, 4, 8]], axis=1)
    print("pressure / nkbT: ",
          np.mean(PHistory[-meanWindow:]), "+ - std: ", np.std(PHistory[-meanWindow:]))
    plt.plot(PHistory)
    plt.savefig('record_'+grepString+'.png')
    plt.show()


analyzeXF(file, meanWindow, 'Col')
analyzeXF(file, meanWindow, 'Bi')
