import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys

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

# Uncomment for LaTex fonts:
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.serif'] = ['Computer Modern']
# mpl.rc('text', usetex=True)

# Use log time scale
ylog = False

cm_subsection = np.linspace(0.0, 1.0, 4)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]

path = 'Output/'#1d_pulse_350V/'
if len(sys.argv) > 1:
    d = sys.argv[1]
    res = sys.argv[2]
    path = 'Output/{}d_res_'.format(d) + res + '/'

x = np.fromfile(path + 'meshx.dat',dtype=float) * 1e3
y = np.fromfile(path + 'meshy.dat',dtype=float) * 1e3
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

tt, xx = np.meshgrid(t,x)
yloc = 0

tloc = [nt/10, nt/2, -1]
tloc[0] = np.argmin(np.abs(t-0.2))

# Figure 1
fig = plt.figure(figsize=(5.25,5.25))
gs = gridspec.GridSpec(2,2)
gs.set_width_ratios([0.97,0.03])
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
ax2 = fig.add_subplot(gs[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
if ylog:
    ax1.set_yscale('log')

ax0.plot(x, f1[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f1[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f1[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, f1[:,yloc,:].T, 30, cmap='viridis')
fig.colorbar(im, cax = ax2)

#ax0.set_ylim([-1.05, 1.05])
#ax1.set_xlim([-0.05, 1.05])
#ax1.set_ylim([-0.5, 10.5])

fig.suptitle('Field')
gs.tight_layout(fig, rect=[0, 0, 1, 1])

# Figure 2
fig2 = plt.figure(figsize=(5.25,5.25))
gs2 = gridspec.GridSpec(2,2)
gs2.set_width_ratios([0.97,0.03])
ax0 = fig2.add_subplot(gs2[0,0])
ax1 = fig2.add_subplot(gs2[1,0], sharex=ax0)
ax2 = fig2.add_subplot(gs2[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax0.set_yscale('log')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
if ylog:
    ax1.set_yscale('log')

ax0.plot(x, f2[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f2[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f2[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, np.log10(f2[:,yloc,:]).T, 30, cmap='viridis')
fig2.colorbar(im, cax = ax2)

#ax0.set_ylim([-0.05, 1.05])
#ax1.set_xlim([-0.05, 1.05])
#ax1.set_ylim([-0.5, 10.5])

fig2.suptitle('Electron Density')
gs2.tight_layout(fig2, rect=[0, 0, 1, 1])

# Figure 3
fig3 = plt.figure(figsize=(5.25,5.25))
gs3 = gridspec.GridSpec(2,2)
gs3.set_width_ratios([0.97,0.03])
ax0 = fig3.add_subplot(gs3[0,0])
ax1 = fig3.add_subplot(gs3[1,0], sharex=ax0)
ax2 = fig3.add_subplot(gs3[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax0.set_yscale('log')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
if ylog:
    ax1.set_yscale('log')

ax0.plot(x, f3[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f3[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f3[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, np.log10(f3[:,yloc,:]).T, 30, cmap='viridis')
fig3.colorbar(im, cax = ax2)

#ax0.set_ylim([-0.05, 1.05])
#ax1.set_xlim([-0.05, 1.05])
#ax1.set_ylim([-0.5, 10.5])

fig3.suptitle('Ion Density')
gs3.tight_layout(fig3, rect=[0, 0, 1, 1])

# Figure 4
fig4 = plt.figure(figsize=(5.25,5.25))
gs4 = gridspec.GridSpec(2,2)
gs4.set_width_ratios([0.97,0.03])
ax0 = fig4.add_subplot(gs4[0,0])
ax1 = fig4.add_subplot(gs4[1,0], sharex=ax0)
ax2 = fig4.add_subplot(gs4[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax0.set_yscale('log')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
if ylog:
    ax1.set_yscale('log')

ax0.plot(x, f4[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f4[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f4[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, np.log10(f4[:,yloc,:]).T, 30, cmap='viridis')
fig4.colorbar(im, cax = ax2)

#ax0.set_ylim([-0.05, 1.05])
#ax1.set_xlim([-0.05, 1.05])
#ax0.set_ylim([1e-2, f4.max()*1.05])

fig4.suptitle('Electron Temperature')
gs4.tight_layout(fig4, rect=[0, 0, 1, 1])

# Figure 5
fig5 = plt.figure(figsize=(5.25,5.25))
gs5 = gridspec.GridSpec(2,2)
gs5.set_width_ratios([0.97,0.03])
ax0 = fig5.add_subplot(gs5[0,0])
ax1 = fig5.add_subplot(gs5[1,0], sharex=ax0)
ax2 = fig5.add_subplot(gs5[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Amplitude')
ax0.set_yscale('log')
ax1.set_xlabel('Position')
ax1.set_ylabel('Time')
if ylog:
    ax1.set_yscale('log')

ax0.plot(x, f5[tloc[0],yloc,:], label='{:.1f}'.format(t[tloc[0]]), color=colors[0])
ax0.plot(x, f5[tloc[1],yloc,:], label='{:.1f}'.format(t[tloc[1]]), color=colors[1])
ax0.plot(x, f5[tloc[2],yloc,:], label='{:.1f}'.format(t[tloc[2]]), color=colors[2])
ax0.legend(bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False, title='Time')

im = ax1.contourf(xx, tt, np.log10(f5[:,yloc,:]).T, 30, cmap='viridis')
fig5.colorbar(im, cax = ax2)

#ax0.set_ylim([-0.05, 1.05])
#ax1.set_xlim([-0.05, 1.05])
#ax1.set_ylim([-0.5, 10.5])

fig5.suptitle('Metastable Density')
gs5.tight_layout(fig5, rect=[0, 0, 1, 1])

# Figure 6
fig6 = plt.figure(figsize=(4.5,5.25))
gs6 = gridspec.GridSpec(2,1)
ax0 = fig6.add_subplot(gs6[0,0])
ax1 = fig6.add_subplot(gs6[1,0], sharex=ax0)

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_ylabel('Voltage')
# ax0.set_yscale('log')
ax1.set_ylabel('Current')
ax1.set_xlabel('Time')
ax1.set_yscale('log')
if ylog:
    ax0.set_xscale('log')
    ax1.set_xscale('log')

ax0.plot(t, Vd, color=colors[0], label='Vd')
ax1.plot(t, np.abs(Id)*1000., color=colors[1], label='Id')
# ax0.plot(t, dt*1000., color=colors[2], label=r'$\Delta$t')
# ax0.legend(loc='best')
ax1.set_ylim([3e-3,3e2])
# ax0.set_xlim([0.05,t[-1]])

gs6.tight_layout(fig6, rect=[0, 0, 1, 1])

if (ny > 1):
    # Phi
    fig = plt.figure(figsize = (4.5,3.5))
    gs = gridspec.GridSpec(1,2, width_ratios=[0.97,0.03])
    ax0 = fig.add_subplot(gs[0,0])
    cbax = plt.subplot(gs[0,1])

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)

    # im = ax0.contourf(y,x,f1[-1,:,:].T, 30)
    im = ax0.pcolormesh(y, x, f1[-1,:,:].T)#, shading='gouraud')
    fig.colorbar(im, cax = cbax)
    ax0.set_ylabel('X')
    ax0.set_xlabel('R')
    fig.suptitle('Field')
    gs.tight_layout(fig, rect=[0, 0, 1, 1])

    # Ne
    fig = plt.figure(figsize = (4.5,3.5))
    gs = gridspec.GridSpec(1,2, width_ratios=[0.97,0.03])
    ax0 = fig.add_subplot(gs[0,0])
    cbax = plt.subplot(gs[0,1])

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)

    # im = ax0.contourf(y,x,f2[-1,:,:].T, 30)
    # im = ax0.pcolormesh(y, x, f2[-1,:,:].T) #, shading='gouraud')
    im = ax0.contourf(y, x, np.log10(f2[-1,:,:].T), 30)
    fig.colorbar(im, cax = cbax)
    ax0.set_ylabel('X')
    ax0.set_xlabel('R')
    fig.suptitle('Density')
    gs.tight_layout(fig, rect=[0, 0, 1, 1])

    # Te
    fig = plt.figure(figsize = (4.5,3.5))
    gs = gridspec.GridSpec(1,2, width_ratios=[0.97,0.03])
    ax0 = fig.add_subplot(gs[0,0])
    cbax = plt.subplot(gs[0,1])

    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)

    # im = ax0.contourf(y,x,f4[-1,:,:].T, 30)
    im = ax0.pcolormesh(y, x, f4[-1,:,:].T)#, shading='gouraud')
    fig.colorbar(im, cax = cbax)
    ax0.set_ylabel('X')
    ax0.set_xlabel('R')
    fig.suptitle('Temperature')
    gs.tight_layout(fig, rect=[0, 0, 1, 1])

plt.show()
