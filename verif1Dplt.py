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
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

cm_subsection = np.linspace(0.0, 1.0, 14)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]

data = np.load('Data/1D_resSwp_nx100_50per.npz')
rfdata =np.genfromtxt('Data/rafatov1D.csv',delimiter=',')

fig = plt.figure(figsize=(5.1,8.0))
gs1 = gridspec.GridSpec(2,2)
gs1.set_height_ratios([0.5,1.0])
gs1.set_width_ratios([0.97,0.03])
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs1[1,0],
                                       hspace=0.2)
gs3 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs1[:,1],
                                       height_ratios=[0.55,1.25,0.25])
ax0 = fig.add_subplot(gs1[0,0])
ax1 = fig.add_subplot(gs2[0,0])
ax2 = fig.add_subplot(gs2[1,0], sharex=ax1)
ax3 = fig.add_subplot(gs3[1])

plt.setp(ax1.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

Id = np.zeros(13)
Vd = np.zeros(13)
t = data['t']
nt = len(t)
t_idx = np.zeros(13, dtype='int')

j = 0
tgt = 49.99
for i in range(nt-1):
    if t[i+1] > tgt:
        t_idx[j] = i
        Id[j] = data['Id'][i]
        Vd[j] = data['Vd'][i]
        print('T={:.2f} Id={:.2e} Vd={:.2f}'.format(t[i], Id[j], Vd[j]))
        j+= 1
        tgt += 50
Id[-1] = data['Id'][-1]
Vd[-1] = data['Vd'][-1]

Id = Id / 1.5e-2

ax0.plot(rfdata[:,0], rfdata[:,1], '--o', color=(0.9,0.4,0.4), label='Ref.',
         markerfacecolor=(1,1,1))
ax0.plot(Id, Vd, color=colors[0], label='This work')
for i in range(13):
    ax0.plot(Id[i], Vd[i], 'o', markerfacecolor=colors[i],
            markeredgecolor=colors[i])

ax0.set_xscale('log')
ax0.set_xlabel(r'Current Density ($mA \, cm^{-2}$)')
ax0.set_ylabel(r'Discharge Voltage ($V$)')
ax0.set_title('(a)')
# ax0.legend(loc=2)
ax0.legend(bbox_to_anchor = (0.85, 0.6), frameon=False)

# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/verif1D_c.eps',dpi=300)

cmap = mpl.cm.viridis_r
norm = mpl.colors.Normalize(vmin=2.7,vmax=7.3)
cb = mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm)
cb.outline.set_visible(False)
cb.ax.set_title(r'log $\Omega$')

for i in range(13):
    if i == 0:
        ax2.plot(t[:t_idx[0]], data['Id'][:t_idx[0]]/ 1.5e-2, color=colors[i])
        ax1.plot(t[:t_idx[0]], data['Vd'][:t_idx[0]], color=colors[i])
    elif i == 12:
        ax2.plot(t[t_idx[i-1]:], data['Id'][t_idx[i-1]:]/ 1.5e-2, color=colors[i])
        ax1.plot(t[t_idx[i-1]:], data['Vd'][t_idx[i-1]:], color=colors[i])
    else:
        ax2.plot(t[t_idx[i-1]:t_idx[i]], data['Id'][t_idx[i-1]:t_idx[i]]/ 1.5e-2, color=colors[i])
        ax1.plot(t[t_idx[i-1]:t_idx[i]], data['Vd'][t_idx[i-1]:t_idx[i]], color=colors[i])

ax1.set_xlim([-20,670])
ax1.set_ylim([160,440])
ax2.set_ylim([3e-4,3e1])
ax2.set_yscale('log')
ax2.set_xlabel(r'Time ($\mu s$)')
ax1.set_ylabel(r'Discharge Voltage ($V$)')
ax2.set_ylabel(r'Current Density ($mA \, cm^{-2}$)')
ax1.set_title('(b)')
ax2.set_title('(c)')
gs1.tight_layout(fig, rect=[0, 0, 1, 1])

plt.savefig('Figures/verif1D.eps',dpi=300)

plt.show()
