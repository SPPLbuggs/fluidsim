import numpy as np
from scipy.io import FortranFile
import glob
import sys

res = ['1e7', '5e6', '3e6', '2e6', '1e6', '5e5', '3e5', '2e5', '1e5', '5e4', '3e4', '2e4', '1e4', '5e3', '3e3', '2e3', '1e3']

for i in range(len(res)):
    path = 'Output/2d_res_' + res[i] + '/'

    x = np.fromfile(path + 'meshx.dat',dtype=float)
    y = np.fromfile(path + 'meshy.dat',dtype=float)
    t = np.fromfile(path + 'time.dat', dtype=float)
    dt = np.fromfile(path + 'dt.dat', dtype=float)
    Id = np.fromfile(path + 'id.dat', dtype=float)
    Vd = np.fromfile(path + 'vd.dat', dtype=float)

    nx = len(x)
    ny = max(len(y),1)
    nt = len(t)

    # Load Data
    temp = np.fromfile(path + 'f1.dat',dtype=float)
    f1 = temp.reshape([nt, ny, nx])

    temp = np.fromfile(path + 'f2.dat',dtype=float)
    f2 = temp.reshape([nt, ny, nx])

    temp = np.fromfile(path + 'f3.dat',dtype=float)
    f3 = temp.reshape([nt, ny, nx])

    temp = np.fromfile(path + 'f4.dat',dtype=float)
    f4 = temp.reshape([nt, ny, nx])
    f4 = f4/f2

    temp = np.fromfile(path + 'f5.dat',dtype=float)
    f5 = temp.reshape([nt, ny, nx])

    path = 'Data/2d_res_' + res[i] + '_100x100.npz'
    np.savez(path, x=x, y=y, t=t, dt=dt, Id=Id, Vd=Vd,
             phi=f1, ne=f2, ni=f3, nt=f4, nm=f5)
