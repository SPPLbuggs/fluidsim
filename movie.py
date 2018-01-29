import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
from matplotlib import animation
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

size = 12
med_size = 13
big_size = 13

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)
#plt.rcParams['figure.autolayout'] = True

path = 'Output/2d_res_1e4/'
x = np.fromfile(path + 'meshx.dat',dtype=float)
y = np.fromfile(path + 'meshy.dat',dtype=float)
t = np.fromfile(path + 'time.dat', dtype=float)

nx = len(x)
ny = len(y)
nt = len(t)

temp = np.fromfile(path + 'f2.dat',dtype=float)
ne = temp.reshape([nt, ny, nx])
# ne = np.log10(ne)

fig = plt.figure()
im = plt.imshow(ne[0,:,:].T, cmap='viridis', animated=True, origin='lower')#, interpolation='bilinear')
plt.xlabel('R')
plt.ylabel('X')
tx = plt.title(r'Electron Density, t = 0.00 $\mu$s')

def anim(i):
    im.set_array(ne[i,:,:].T)
    tx.set_text(r'Electron Density : t = {:.2f} $\mu$s'.format(t[i]))
    im.autoscale()

ani = animation.FuncAnimation(fig, anim, frames=nt, interval=50)

plt.show()
#ani.save('figures/anim.gif', dpi = 100)
