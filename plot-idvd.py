import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from scipy.io import FortranFile

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
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

cm_subsection = np.linspace(0.0, 1.0, 22)
colors = [ mpl.cm.viridis_r(x) for x in cm_subsection ]

reses = ['1e7', '5e6', '3e6', '2e6', '1e6', '5e5', '3e5', '2e5', '1e5', '5e4', '3e4', '2e4', '1e4', '5e3', '3e3', '2e3', '1e3']#, '5e2']

nr = len(reses)

Id_f = np.zeros(nr)
Vd_f = np.zeros(nr)
Id_f2 = np.zeros(nr)
Vd_f2 = np.zeros(nr)

fig = plt.figure(figsize=(5.25,5.25))
gs = gridspec.GridSpec(2,2)
gs.set_width_ratios([0.97,0.03])
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
ax2 = fig.add_subplot(gs[:,1])

plt.setp(ax0.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=2.7,vmax=7.3)
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
cb.outline.set_visible(False)
cb.ax.set_title(r'log $\Omega$')

#ax0.set_ylim([50,550])
ax0.set_ylabel(r'Discharge Voltage [$V$]')

#ax1.set_xlim([5e-2,3e2])
#ax1.set_ylim([7e-6, 4e-1])
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'Time [$\mu s$]')
ax1.set_ylabel(r'Discharge Current [$A$]')

for i in range(nr):
    res = reses[i]
    path = 'Output/1d_res_' + res

    t = np.fromfile(path + '/time.dat', dtype=float)
    Id = np.fromfile(path + '/id.dat', dtype=float)
    Vd = np.fromfile(path + '/vd.dat', dtype=float)

    Id_f[i] = Id[-1]
    Vd_f[i] = Vd[-1]

    ax0.plot(t, Vd, color=colors[i], label=(res))
    ax1.plot(t, abs(Id), color=colors[i])

#ax0.annotate(r'$\Omega$', fontsize=14, color=(0.2,0.2,0.2),
#             xy=(0.8,0.48), xycoords='axes fraction',
#             xytext=(0.65,0.8), textcoords='axes fraction',
#             arrowprops=dict(arrowstyle='->', linestyle='dashed',
#                             connectionstyle='arc3,rad=-0.3',color=(0.3,0.3,0.3)))

gs.tight_layout(fig, rect=[0, 0, 1, 1])
fig.suptitle('(a) Standard Model', fontsize=big_size)


#*** Figure 2 ***#
fig2 = plt.figure()
gs2 = gridspec.GridSpec(1,1)
ax = fig2.add_subplot(gs2[0])

# Rafatov Data
rdata = np.genfromtxt('Data/rafatovData.csv',delimiter=',')
ax.plot(rdata[:,0], rdata[:,1], '--^', color=(0.9,0.3,0.3), label='Ref.',
        markerfacecolor=(1,1,1))

ax.plot(Id_f, Vd_f, '--o', color=(0.5,0.5,0.5), label='1d Std.')
for i in range(nr):
    ax.plot(Id_f[i], Vd_f[i], 'o', markerfacecolor=colors[i],
             markeredgecolor=colors[i], markersize = 4)

ax.set_ylim([0,500])
ax.set_xscale('log')
fig2.suptitle('(c) I-V Summary', fontsize=big_size)
ax.set_ylabel(r'Discharge Voltage [$V$]')
ax.set_xlabel(r'Discharge Current [$A$]')
ax.set_xscale('log')
plt.legend(frameon=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
gs2.tight_layout(fig2, rect=[0, 0, 1, 1])
plt.show()
