import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
from scipy.interpolate import spline

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

data = np.load('Data/2D_350AC_nx80_ny80_80us.npz')

x = data['x'] * 1e3
y = data['y'] * 1e3
t = data['t']
Vd = data['Vd']
Id = data['Id']

ne = data['ne']
nt = data['nt']
ni = data['ni']
nm = data['nm']
ph = data['phi']

tgt = [78.55, 78.8, 79.02]
ts = np.zeros(len(tgt), dtype='int')
j = 0
for i in range(len(t)-1):
    if t[i+1] > tgt[j]:
        ts[j] = i
        j = j+1

    if j == len(tgt):
        break

fig = plt.figure(figsize=(5.25,5))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.0,0.05], height_ratios=[1.0,0.75])
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0])
axc = fig.add_subplot(gs[0,1])

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.setp(ax0.get_xticklabels(), visible=False)

j = 2

v = np.linspace(13.9,18.1,25)
im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
tck = [14.0, 15.0, 16.0, 17.0, 18]
tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
cb = plt.colorbar(im, cax=axc, ticks=tck)
cb.outline.set_visible(False)
cb.ax.set_yticklabels(tck_lbls)
cb.ax.set_title(r'[$m^{-3}$]')
cb.solids.set_edgecolor("face")

for c in im.collections:
    c.set_edgecolor("face")
    c.set_linewidth(0.000000000001)

ax0.set_ylabel(r'R [$mm$]')
ax0.set_title(r'Electron Density')

xnew = np.linspace(x[0],x[-1],300)
ne_s = spline(x,ne[ts[j],0,:],xnew)
ni_s = spline(x,ni[ts[j],0,:],xnew)
nm_s = spline(x,nm[ts[j],0,:],xnew)

ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
ax1.set_yscale('log')
ax1.set_xlabel(r'X [$mm$]')
ax1.set_ylabel(r'Density [$m^{-3}$]')
ax1.set_ylim([2e15,5e19])
ax1.legend(bbox_to_anchor=(1.0,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)

gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig('Figures/2dAC.eps', dpi = 300)


## Figure 2:: IV
fig = plt.figure(figsize=(5,3))
gs = gridspec.GridSpec(1,1)
ax0 = fig.add_subplot(gs[0])

#ax.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

tnew = np.linspace(t[870],t[-1],300)
Id_s = spline(t[870:],Id[870:]*1000.,tnew)
Vd_s = spline(t[870:],Vd[870:],tnew)

ax0.plot(tnew, Vd_s, label=r'V$_\mathrm{d}$')
ax1 = ax0.twinx()
next(ax1._get_lines.prop_cycler)
ax1.plot(tnew, Id_s, label=r'I$_\mathrm{d}$')
plt.legend()

# ax0.set_xscale('log')
# ax0.set_xlim([0.008,40])
# ax0.set_ylim([-150,1350])
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

ax1.legend(loc=1, bbox_to_anchor=(1.0,1.1), frameon=False)
ax0.legend(loc=2, bbox_to_anchor=(0.0,1.1), frameon=False)

plt.annotate('', xy=(0.02, 0.73+0.1), xycoords='axes fraction',
             xytext=(0.157, 0.73), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=0.27",
             color=color2[0]))

plt.annotate('', xy=(0.99, 0.44), xycoords='axes fraction',
             xytext=(0.9, 0.36), textcoords='axes fraction',
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.27",
             color=color2[1]))

gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig('Figures/2dAC_IV.eps', dpi = 300)

plt.show()
