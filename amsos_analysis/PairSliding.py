import Util.AMSOS as am
import numpy as np
import scipy as sp
import scipy.optimize as so
from numpy.lib.recfunctions import structured_to_unstructured

import matplotlib as mpl
import matplotlib.pyplot as plt

param = am.ParamBase(text='fit velocity for a pair of sliding mts')


def process_frame(file):
    frame = am.FrameAscii(file, sort=False, info=True)
    TList = frame.TList
    Tm = structured_to_unstructured(TList[['mx', 'my', 'mz']])
    Tp = structured_to_unstructured(TList[['px', 'py', 'pz']])
    Tc = (Tm+Tp)*0.5
    distm01 = np.linalg.norm(Tm[0]-am.point_line_proj(Tm[0], Tm[1], Tp[1]))
    distp01 = np.linalg.norm(Tp[0]-am.point_line_proj(Tp[0], Tm[1], Tp[1]))
    distm10 = np.linalg.norm(Tm[1]-am.point_line_proj(Tm[1], Tm[0], Tp[0]))
    distp10 = np.linalg.norm(Tp[1]-am.point_line_proj(Tp[1], Tm[0], Tp[0]))
    distc01 = np.linalg.norm(Tc[1]-am.point_line_proj(Tc[0], Tm[1], Tp[1]))
    distc10 = np.linalg.norm(Tc[0]-am.point_line_proj(Tc[1], Tm[0], Tp[0]))

    return np.array([distm01, distp01, distm10, distp10, distc01, distc10])


nframes = len(param.syfiles)
data = np.zeros((nframes, 6))
for i, f in enumerate(param.syfiles):
    data[i, :] = process_frame(f)

np.savetxt(fname='pairsliding_dist.txt', X=data,
           header='distm01, distp01, distm10, distp10, distc01, distc10')

dt = param.config['timeSnap']*param.stride
t = dt*np.linspace(0, nframes, nframes)
x = (data[:, 4]+data[:, 5])/2
z = (data[:, 0]+data[:, 1]+data[:, 2]+data[:, 3])/4
# fit and plot
plt.figure()
plt.plot(t, data[:, 0], label='distm01')
plt.plot(t, data[:, 1], label='distp01')
plt.plot(t, data[:, 2], label='distm10')
plt.plot(t, data[:, 3], label='distp10')
plt.plot(t, data[:, 4], label='distc01')
plt.plot(t, data[:, 5], label='distc10')


def f_fit(time, a, b):
    return a*time+b


start_fit = nframes//2
sol = so.curve_fit(f_fit, t[start_fit:], x[start_fit:], p0=[1, 1])
print(sol)
z_dist = np.mean(z[start_fit:])
print(z_dist)

ufit = sol[0][0]
plt.plot(t[start_fit:], f_fit(t[start_fit:],
         ufit, sol[0][1]), label='u_fit={:g} um/s'.format(ufit))
plt.plot(t, z, label='distance = {:g} um'.format(z_dist))

plt.legend()
plt.savefig('PairSligin.png', dpi=300)
