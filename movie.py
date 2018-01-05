import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm
from scipy.io import FortranFile
import glob
from matplotlib import animation
import types
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

size = 12
med_size = 13
big_size = 14

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)
#plt.rcParams['figure.autolayout'] = True

cm_subsection = np.linspace(0.0, 1.0, 4)
colors = [ cm.viridis(x) for x in cm_subsection ]

path = 'output/'
x = np.fromfile('Output/meshx.dat',dtype=float)
y = np.fromfile('Output/meshy.dat',dtype=float)
t = np.fromfile('Output/time.dat', dtype=float)

nx = len(x)
ny = len(y)
ts = len(t)

Ex = np.zeros([nx,ny,ts])
Ey = np.zeros([nx,ny,ts])
ne = np.zeros([nx,ny,ts])
ni = np.zeros([nx,ny,ts])
nt = np.zeros([nx,ny,ts])
nm = np.zeros([nx,ny,ts])

temp = np.fromfile('Output/f2.dat',dtype=float)
ne = temp.reshape([ts, ny, nx])

interp = 'bilinear'
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))
im = NonUniformImage(ax, interpolation=interp, cmap='viridis',
                     extent=(y[0], y[-1], x[0], x[-1]), zorder= -20)
#ax.set_rasterization_zorder(-10)
im.set_data(y, x, ne[0,:,:].T)
ax.images.append(im)
ax.set_xlim(y[0]-0.1, y[-1]+0.1)
ax.set_ylim(x[0]-0.1, x[-1]+0.1)
tx = plt.title(r'$n_e$ : time = 0 $\mu$s')
ax.set_xlabel('r [mm]')
ax.set_ylabel('z [mm]')
plt.tight_layout()

def anim(i):
    im.set_data(y, x, ne[i,:,:].T)
    ax.set_xlim(y[0]-0.1, y[-1]+0.1)
    ax.set_ylim(x[0]-0.1, x[-1]+0.1)
    tx.set_text(r'$n_e$ : time = {:.2f} $\mu$s'.format(t[i]))
    im.autoscale()

ani = animation.FuncAnimation(fig, anim, frames = ts, interval = 50)

plt.show()
#ani.save('figures/anim.gif', dpi = 100)
