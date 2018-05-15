import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm
from scipy.io import FortranFile
import glob
from matplotlib import animation
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

path = 'Output/'

y = np.fromfile(path + 'meshy.dat',dtype=float) * 1e3

ny = max(len(y),1)
nx = ny * 2

temp = np.fromfile('Output/Hz.dat',dtype=float)
nt = len(temp) / ny / nx
Hz = temp.reshape([nt, ny, nx])

fig = plt.figure(figsize=(8,6))
im = plt.imshow(Hz[0,:,:], animated=True)
plt.xlabel('X')
plt.ylabel('Y')
tx = plt.title('Hz : time = 0')
plt.colorbar()

def animate(i):
    im.set_array(Hz[i,:,:])
    tx.set_text('Hz : time = {}'.format(i))
    im.autoscale()

ani = animation.FuncAnimation(fig, animate, frames = nt, interval = 75)

plt.show()
