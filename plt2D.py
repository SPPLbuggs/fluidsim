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

cm_subsection = np.linspace(0.0, 1.0, 4)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]
color2= plt.rcParams['axes.prop_cycle'].by_key()['color']

# Uncomment for LaTex fonts:
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

data = np.load('Data/2D_1250pulse_nx100_ny100_31us.npz')

x = data['x'] * 1e3
y = data['y'] * 1e3
t = data['t']
Vd = data['Vd']
Id = data['Id']

ne = np.log10(data['ne'])
nt = data['nt']
ni = data['ne']
ph = data['phi']

tgt = [0.111, 0.135, 0.501, 30.1]
ts = np.zeros(len(tgt), dtype='int')
j = 0
for i in range(len(t)-1):
    if t[i+1] > tgt[j]:
        ts[j] = i
        j = j+1

    if j == len(tgt):
        break

# ts = [240, 280, 350, 550, -1]

# Figure 1: Ne in time
fig = plt.figure(figsize=(8,8))
gs1 = gridspec.GridSpec(2,2, hspace=0.75, wspace=0.75)
gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs1[0,0],
      wspace=0.05, hspace=0.15, width_ratios=[10.0,0.5])
gs3 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs1[0,1],
      wspace=0.05, hspace=0.15, width_ratios=[10.0,0.5])
gs4 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs1[1,0],
      wspace=0.05, hspace=0.15, width_ratios=[10.0,0.5])
gs5 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs1[1,1],
      wspace=0.05, hspace=0.15, width_ratios=[10.0,0.5])

##First point in time
ax0  = fig.add_subplot(gs2[0,0])
ax0c = fig.add_subplot(gs2[0,1])
ax1  = fig.add_subplot(gs2[1,0])
ax1c = fig.add_subplot(gs2[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)

ax0.set_title(r'(a) {:.2f} $\mu$s'.format(t[ts[0]]))
ax0.set_ylabel(r'r $(mm)$')
ax1.set_ylabel(r'r $(mm)$')
ax1.set_xlabel(r'x $(mm)$')

vmin = 10.8
vmax = 17.2
v = np.linspace(vmin, vmax, 25)
tck_ne = [11, 13, 15, 17]
tck_lbl_ne = [r'$10^{11}$',r'$10^{13}$',r'$10^{15}$',r'$10^{17}$']
# ne
im = ax0.contourf(x, y, np.maximum(ne[ts[0],:,:],vmin), v)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
cbar = plt.colorbar(im, cax=ax0c, ticks=tck_ne)
cbar.ax.set_yticklabels(tck_lbl_ne)
cbar.set_label(r'n$_\mathrm{e}$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)

# phi
im = ax1.contourf(x, y, ph[ts[0],:,:], 25)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
tcks = np.arange(0,1201,300)
# tck_lbls = ['{}'.format(int(val)) for val in tcks]
cbar = plt.colorbar(im, cax=ax1c, ticks=tcks)
cbar.set_label(r'$\phi$')
# cbar.ax.set_yticklabels(tck_lbls)
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)

## Second point in time
ax0  = fig.add_subplot(gs3[0,0])
ax0c = fig.add_subplot(gs3[0,1])
ax1  = fig.add_subplot(gs3[1,0])
ax1c = fig.add_subplot(gs3[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)

ax0.set_title(r'(b) {:.2f} $\mu$s'.format(t[ts[1]]))
ax0.set_ylabel(r'r [$mm$]')
ax1.set_ylabel(r'r [$mm$]')
ax1.set_xlabel(r'x [$mm$]')

# ne
im = ax0.contourf(x, y, np.maximum(ne[ts[1],:,:],vmin), v)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
cbar = plt.colorbar(im, cax=ax0c, ticks=tck_ne)
cbar.ax.set_yticklabels(tck_lbl_ne)
cbar.set_label(r'n$_\mathrm{e}$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)
# phi
im = ax1.contourf(x, y, ph[ts[1],:,:], 25)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
tcks = np.arange(50,1251,300)
tck_lbls = ['{}'.format(int(val)) for val in tcks]
cbar = plt.colorbar(im, cax=ax1c, ticks=tcks)
# cbar.ax.set_yticklabels(tck_lbls)
cbar.set_label(r'$\phi$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)

## Third point in time
ax0  = fig.add_subplot(gs4[0,0])
ax0c = fig.add_subplot(gs4[0,1])
ax1  = fig.add_subplot(gs4[1,0])
ax1c = fig.add_subplot(gs4[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)

ax0.set_title(r'(c) {:.2f} $\mu$s'.format(t[ts[2]]))
ax0.set_ylabel(r'r [$mm$]')
ax1.set_ylabel(r'r [$mm$]')
ax1.set_xlabel(r'x [$mm$]')

# ne
im = ax0.contourf(x, y, np.maximum(ne[ts[2],:,:],vmin), v)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
cbar = plt.colorbar(im, cax=ax0c, ticks=tck_ne)
cbar.ax.set_yticklabels(tck_lbl_ne)
cbar.set_label(r'n$_\mathrm{e}$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)
# phi
im = ax1.contourf(x, y, ph[ts[2],:,:], 25)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
tcks = np.linspace(1.0,13,6)
tck_lbls = ['{}'.format(int(val)) for val in tcks]
cbar = plt.colorbar(im, cax=ax1c, ticks=tcks)
# cbar.ax.set_yticklabels(tck_lbls)
cbar.set_label(r'$\phi$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)

## Fourth point in time
ax0  = fig.add_subplot(gs5[0,0])
ax0c = fig.add_subplot(gs5[0,1])
ax1  = fig.add_subplot(gs5[1,0])
ax1c = fig.add_subplot(gs5[1,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax0.set_title(r'(d) {:.1f} $\mu$s'.format(t[ts[3]]))
ax0.set_ylabel(r'r [$mm$]')
ax1.set_ylabel(r'r [$mm$]')
ax1.set_xlabel(r'x [$mm$]')

# ne
im = ax0.contourf(x, y, np.maximum(ne[ts[3],:,:],vmin), v)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
cbar = plt.colorbar(im, cax=ax0c, ticks=tck_ne)
cbar.ax.set_yticklabels(tck_lbl_ne)
cbar.set_label(r'n$_\mathrm{e}$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)
# phi
im = ax1.contourf(x, y, ph[ts[3],:,:], 25)
for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)
tcks = np.arange(1,3.1,0.5)
tck_lbls = ['{}'.format(int(val)) for val in tcks]
cbar = plt.colorbar(im, cax=ax1c, ticks=tcks)
# cbar.ax.set_yticklabels(tck_lbls)
cbar.set_label(r'$\phi$')
cbar.solids.set_edgecolor("face")
cbar.outline.set_visible(False)

gs1.tight_layout(fig, rect=[0, 0, 1, 1])

plt.savefig('Figures/2dpulse.eps')

## Current/Voltage Figure ###
fig = plt.figure(figsize=(5,3))
gs = gridspec.GridSpec(1,1)
ax0 = fig.add_subplot(gs[0])

#ax.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

ax0.plot(t, Vd, color=color2[0], label=r'V$_\mathrm{d}$')
ax1 = ax0.twinx()
ax1.plot(t, Id * 1000, color=color2[1], label=r'I$_\mathrm{d}$')

ax0.set_xscale('log')
ax0.set_xlim([0.008,40])
ax0.set_ylim([-150,1350])
# ax0.set_yticks(np.arange(25,126,25))
# ax0.set_ylim([20,130])

# ax0.tick_params('y', colors=colors[0])
ax0.set_ylabel(r'Discharge Voltage [$V$]')
ax0.set_xlabel(r'Time [$\mu s$]')

ax1.set_ylabel(r'Discharge Current [$mA$]')
# ax1.tick_params('y', colors=colors[2])
# ax1.set_yscale('log')
#ax1.set_ylim([200,1000])
ax0.spines['top'].set_visible(False)
ax1.spines['top'].set_visible(False)

plt.annotate(r'(a)', fontsize=12, xy=(0.235,1.0), xycoords='axes fraction')
s = np.linspace(-150,1350,10)
v = np.ones(10)*0.105
ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)

plt.annotate(r'(b)', fontsize=12, xy=(0.34,1.0), xycoords='axes fraction')
v = np.ones(10)*0.13
ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)

plt.annotate(r'(c)', fontsize=12, xy=(0.51,1.0), xycoords='axes fraction')
v = np.ones(10)*0.5
ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)

plt.annotate(r'(d)', fontsize=12, xy=(0.9,1.0), xycoords='axes fraction')
v = np.ones(10)*30
ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)

plt.annotate('', xy=(0.095, 0.42), xycoords='axes fraction',
             xytext=(0.215, 0.27), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=0.42",
             color=color2[0]))

plt.annotate('', xy=(0.34+0.15, 0.42), xycoords='axes fraction',
             xytext=(0.34, 0.23), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.36",
             color=color2[1]))

ax1.legend(loc=1, bbox_to_anchor=(0.97,0.95), frameon=False)
ax0.legend(loc=2, bbox_to_anchor=(0.0,0.95), frameon=False)

gs.tight_layout(fig, rect=[0, 0, 1, 1])

plt.savefig('Figures/2dpulseIV.eps')

plt.show()
