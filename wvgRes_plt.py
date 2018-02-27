import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
import time
from scipy import integrate

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

cm_subsection = np.linspace(0.0, 1.0, 10)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]
color2= plt.rcParams['axes.prop_cycle'].by_key()['color']

# Uncomment for LaTex fonts:
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

data = np.load('Data/wvgResData.npz')

dw = data['dw']
t = data['t']
nt = len(t)

xl = 1.0
yl = 2.0
zl = 2.0

wr = np.pi * 2.9989e8 * (1e4 / yl**2 + 1e4 / zl**2)**0.5
Qa = 2000.
wT0 = wr * (((1 - (1j + 1.0) / Qa + 0)**0.5 - 1).real + 1)
nf = 50
wT = wT0 + np.linspace(0,60,nf) * 2e6 * np.pi
T = np.zeros([nt, nf])
Q = np.zeros(nt)

for n in range(nt):
    Q[n]= -(dw[n].real + 1.0) / dw[n].imag / 2.0
    w = wr * (dw[n].real + 1)
    for l in range(nf):
        T[n,l] = 1.0 / Qa**2 / ((wT[l]/w - w/wT[l])**2 + (1.0/Q[n])**2)

fig = plt.figure(figsize=(5.25,3))
gs = gridspec.GridSpec(1,1)
ax0 = fig.add_subplot(gs[0,0])
ax1 = ax0.twinx()

ax0.spines['top'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel(r'$\Delta f_\mathrm{res}$ [$MHz$]')
ax1.set_ylabel(r'Quality Factor')
ax0.set_xlabel(r'Time [$\mu \,s$]')
ax0.set_xscale('log')

ax0.plot(t, dw.real*wr/2e6/np.pi, color=color2[0], label=r'$\Delta f_\mathrm{res}$')
ax1.plot(t, Q, color=color2[1], label=r'Q')

ax0.legend(frameon=False, loc=2, bbox_to_anchor = (0.0, 1.075))
ax1.legend(frameon=False, loc=1, bbox_to_anchor = (1.0, 1.075))
# ax0.set_ylim([-10,275])
# ax1.set_ylim([-100, 2250])

plt.annotate('', xy=(0.17-0.13, 0.04+0.15), xycoords='axes fraction',
             xytext=(0.17, 0.03), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=0.38",
             color=color2[0]))

plt.annotate('', xy=(0.77+0.15, 0.56), xycoords='axes fraction',
             xytext=(0.77, 0.725), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=0.38",
             color=color2[1]))

plt.xlim([2e-2,35])
gs.tight_layout(fig, rect=[0, 0, 1, 1])

fig = plt.figure(figsize=(5.25,3))
gs = gridspec.GridSpec(1,2, hspace=0.0)
gs.set_width_ratios([1.0,0.03])
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[0,1])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=-10,vmax=210)
cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, ticks=np.arange(0,201,50))
cb.outline.set_visible(False)
cb.ax.set_title(r'$\Delta f$ [$MHz$]')

for l in range(0,nf,nf/5):
    ax0.plot(t, 10*np.log10(T[:,l]), color=colors[l/5],
             label=r'{}'.format(int((wT[l]-wr)/1e6/2.0/np.pi)))
ax0.set_xlabel(r'Time [$\mu$s]')
ax0.set_ylabel(r'Transmission [$dB$]')
# ax0.set_xticks(np.arange(0,22,5))
# ax0.set_xlim([-2,22])
ax0.set_xscale('log')
ax0.set_xlim([2e-2,5e-1])
gs.tight_layout(fig, rect=[0, 0, 1, 1])


### Plot 3
fig = plt.figure(figsize=(5.25,3))
gs = gridspec.GridSpec(1,1, hspace=0.0)
# gs.set_width_ratios([1.0,0.03])
ax0 = fig.add_subplot(gs[0,0])
# ax1 = fig.add_subplot(gs[0,1])

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

# cmap = mpl.cm.viridis
# norm = mpl.colors.Normalize(vmin=-10,vmax=210)
# cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, ticks=np.arange(0,201,50))
# cb.outline.set_visible(False)
# cb.ax.set_title(r'$\Delta f$ [$MHz$]')

tgt = np.logspace(np.log10(80e-3),np.log10(t[-1]),10)
ns = np.zeros(len(tgt), dtype='int')
for i in range(len(tgt)):
    ns[i] = np.argmin((tgt[i] - t)**2)

for j in range(len(ns)):
    n = ns[j]
    ax0.plot((wT - wT0)/2e6/np.pi, 10*np.log10(T[n,:]), color=colors[j])
ax0.set_xlabel(r'$\Delta$f [MHz]')
ax0.set_ylabel(r'Transmission [$dB$]')
# ax0.set_xticks(np.arange(0,22,5))
# ax0.set_xlim([-2,22])
# ax0.set_xscale('log')
# ax0.set_xlim([2e-2,5e-1])
gs.tight_layout(fig, rect=[0, 0, 1, 1])

plt.show()
