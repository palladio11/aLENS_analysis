import Util.AMSOS as am
import glob
from numpy.lib.recfunctions import structured_to_unstructured
import numpy as np
import matplotlib.pyplot as plt
import json as js
import dask.distributed as dd
import matplotlib
matplotlib.use('Agg')

files = am.getFileListSorted('BindDensity/*.json')
nfiles = len(files)
data = []
u = []
sa = []
sb = []
db = []

for i in range(nfiles//2):
    with open(files[i], 'r') as f:
        item = js.load(f)
        u.append(item['uCount'])
        sa.append(item['saCount'])
        sb.append(item['sbCount'])
        db.append(item['dCount'])

for i in range(nfiles//2, nfiles):
    with open(files[i], 'r') as f:
        item = js.load(f)
        data.extend(item['head1_to_foot_of0_dimensionless'])
        u.append(item['uCount'])
        sa.append(item['saCount'])
        sb.append(item['sbCount'])
        db.append(item['dCount'])

fig = plt.figure()
ax = fig.add_subplot(111)
xlim = 10  # dimensionless
ax.hist(data, bins=1000, range=[-xlim, xlim],
        density=True, alpha=0.5)
plt.savefig('BindDensity_Avg.png', dpi=600)
plt.close(fig)

fig = plt.figure()
ax = fig.add_subplot(111)
t = range(len(u))
ax.semilogy(t, u, label='unbound')
ax.semilogy(t, sa, label='only head a bound')
ax.semilogy(t, sb, label='only head b bound')
ax.semilogy(t, db, label='doubly bound')
ax.legend()
plt.savefig('BindNumber.png', dpi=600)
plt.close(fig)
