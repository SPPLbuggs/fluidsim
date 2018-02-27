import numpy as np
from scipy.io import FortranFile
import glob
import sys

path = 'Output/2d_pulse_1500V/'

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

path = 'Data/data.npz'
np.savez(path, x=x, y=y, t=t, dt=dt, Id=Id, Vd=Vd,
         phi=f1, ne=f2, ni=f3, nt=f4, nm=f5)
