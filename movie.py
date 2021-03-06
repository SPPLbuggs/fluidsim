import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
from matplotlib import animation
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

# Uncomment for LaTex fonts:
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

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
# plt.rcParams['figure.autolayout'] = True

path = 'Output/10Torr/1150V_60x40/'
# path = 'Output/'
x = np.fromfile(path + 'meshx.dat',dtype=float) *1e3
y = np.fromfile(path + 'meshy.dat',dtype=float) *1e3
t = np.fromfile(path + 'time.dat', dtype=float)

nx = len(x)
ny = len(y)
nt = len(t)

temp = np.fromfile(path + 'f2.dat',dtype=float)
ne = temp.reshape([nt, ny, nx])
# ne = np.log10(ne)

ne2 = np.zeros([nt, ny*2, nx])
ne2[:,:ny,:] = ne[:,::-1,:]
ne2[:,ny:,:] = ne

# ne2 = np.log10(ne2)
y2 = np.concatenate([-1*y[::-1],y])

fig = plt.figure(figsize=(7,3.5))
gs = gridspec.GridSpec(1,1)
ax0 = fig.add_subplot(gs[0,0])

im = ax0.imshow(ne2[0,:,:].T, cmap='viridis', animated=True, origin='lower') #interpolation='bicubic')#, vmin = 13.8, vmax = 17.2) vmin=14.7, vmax=18.3,
ax0.set_xlabel(r'R [$mm$]')
ax0.set_ylabel(r'X [$mm$]')
tx = ax0.set_title(r'Electron Density, t = 0.00 $\mu$s')
# tck = [14, 15, 16, 17]
# tck_lbl = [r'$10^{14}$',r'$10^{15}$',r'$10^{16}$',r'$10^{17}$']
# cbar = fig.colorbar(im, ticks=tck)
# cbar.ax.set_yticklabels(tck_lbl)
# cbar.ax.set_title(r'[$m^{-3}$]')
# ax0.colorbar()

y_tck = np.linspace(0,len(x),5,dtype='int')
y_lbl = np.linspace(0,10,5,dtype='int')
ax0.set_yticks(y_tck)
ax0.set_yticklabels(y_lbl)

x_tck = np.linspace(0,2*len(y),9,dtype='int')
x_lbl = np.abs(np.linspace(-10,10,9,dtype='int'))
ax0.set_xticks(x_tck)
ax0.set_xticklabels(x_lbl)
fig.colorbar(im, ax=ax0)

gs.tight_layout(fig, rect=[0, 0, 1, 1])

nt = 36
# tgt = np.logspace(np.log10(3e-2), np.log10(t[-1]), nt)
tgt = np.linspace(0.15, 0.5, nt)
ii = np.unique([np.argmax(t - i > 0) for i in tgt])

def anim(i):
    # j = 5*(i+100)
    j = ii[i]
    im.set_array(ne2[j,:,:].T)
    tx.set_text(r'Electron Density : t = {:.2f} $\mu$s'.format(t[j]))
    im.autoscale()
    return im

ani = animation.FuncAnimation(fig, anim, frames=nt, interval=150, blit=False)

# fig = plt.figure(figsize=(6, 3))
# gs = gridspec.GridSpec(1, 1)
# ax = fig.add_subplot(gs[0,0])
# quad = ax.pcolormesh(y2, x, ne2[0,:,:].T, vmin=14.7, vmax=18.3) #shading='gouraud',
# ax.set_xlabel('R')
# ax.set_ylabel('X')
# cb = fig.colorbar(quad,ax=ax)
# tx = plt.title(r'Electron Density, t = 0.00 $\mu$s')
#
# def init():
#     quad.set_array([])
#     return quad1
#
# nt = 36
# # tgt = np.logspace(np.log10(3e-2), np.log10(t[-1]), nt)
# tgt = np.linspace(0.15, 0.5, nt)
# ii = np.unique([np.argmax(t - i > 0) for i in tgt])
# def animate(i):
#     j = ii[i]
#     quad.set_array((ne2[j,:,:].T).ravel())
#     tx.set_text(r'Electron Density : t = {:.2f} $\mu$s'.format(t[j]))
#     # quad.autoscale()
#     return quad
#
# gs.tight_layout(fig)
#
# anim = animation.FuncAnimation(fig,animate,frames=nt,interval=150,blit=False,repeat=True)
plt.show()

# Set up formatting for the movie files
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=20, bitrate=900)
# ani.save('anim.mp4', writer=writer)
