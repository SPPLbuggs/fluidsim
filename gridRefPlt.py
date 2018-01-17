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

cm_subsection = np.linspace(0.0, 1.0, 5)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]

data_80 = np.load('Data/DC350V_nx80.npz')
data_160 = np.load('Data/DC350V_nx160.npz')
data_240 = np.load('Data/DC350V_nx240.npz')
data_120 = np.load('Data/DC350V_nx120_nonunif.npz')

# Figure 1: Ion and Electron Density
# fig = plt.figure(figsize=(4.5,5.25))
# gs = gridspec.GridSpec(2,1)
# ax0 = fig.add_subplot(gs[0,0])
# ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
#
# plt.setp(ax0.get_xticklabels(), visible=False)
# ax0.spines['right'].set_visible(False)
# ax0.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['top'].set_visible(False)
#
# ax0.plot( data_80['x']/ 1e2,  data_80['ne'][-1,0,:], label='80 pts', color=colors[0])
# ax0.plot(data_160['x']/ 1e2, data_160['ne'][-1,0,:], label='160 pts', color=colors[1])
# ax0.plot(data_240['x']/ 1e2, data_240['ne'][-1,0,:], label='240 pts', color=colors[2])
# ax0.plot(data_120['x']/ 1e2, data_120['ne'][-1,0,:], label=r'120 pts$^*$', color=colors[3])
#
# ax1.plot(data_80['x'] / 1e2, data_80['ni'][-1,0,:], label='80 pts', color=colors[0])
# ax1.plot(data_160['x']/ 1e2, data_160['ni'][-1,0,:], label='160 pts', color=colors[1])
# ax1.plot(data_240['x']/ 1e2, data_240['ni'][-1,0,:], label='240 pts', color=colors[2])
# ax1.plot(data_120['x']/ 1e2, data_120['ni'][-1,0,:], label=r'120 pts$^*$', color=colors[3])
# ax0.legend()
# ax0.set_yscale('log')
# ax1.set_yscale('log')
# ax0.set_ylim([8e15,4e18])
# ax1.set_ylim([8e15,4e18])
# ax1.set_xlabel(r'X [$mm$]')
# ax0.set_ylabel(r'Electron Density [$cm^{-3}$]')
# ax1.set_ylabel(r'Ion Density [$cm^{-3}$]')
# gs.tight_layout(fig, rect=[0, 0, 1, 1])


# Figure 2: Ion and Electron Density and potential
fig = plt.figure(figsize=(4.5,7.5))
gs = gridspec.GridSpec(3,1)
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0], sharex=ax0)
ax2 = fig.add_subplot(gs[2,0], sharex=ax0)

plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax0.plot( data_80['x']/ 1e2,  data_80['ne'][-1,0,:], label='80 pts', color=colors[0])
ax0.plot(data_160['x']/ 1e2, data_160['ne'][-1,0,:], label='160 pts', color=colors[1])
ax0.plot(data_240['x']/ 1e2, data_240['ne'][-1,0,:], label='240 pts', color=colors[2])
ax0.plot(data_120['x']/ 1e2, data_120['ne'][-1,0,:], label=r'120 pts$^*$', color=colors[3])

ax1.plot(data_80['x'] / 1e2, data_80['ni'][-1,0,:], label='80 pts', color=colors[0])
ax1.plot(data_160['x']/ 1e2, data_160['ni'][-1,0,:], label='160 pts', color=colors[1])
ax1.plot(data_240['x']/ 1e2, data_240['ni'][-1,0,:], label='240 pts', color=colors[2])
ax1.plot(data_120['x']/ 1e2, data_120['ni'][-1,0,:], label=r'120 pts$^*$', color=colors[3])

ax2.plot(data_80['x'] / 1e2, data_80['phi'][-1,0,:], label='80 pts', color=colors[0])
ax2.plot(data_160['x']/ 1e2, data_160['phi'][-1,0,:], label='160 pts', color=colors[1])
ax2.plot(data_240['x']/ 1e2, data_240['phi'][-1,0,:], label='240 pts', color=colors[2])
ax2.plot(data_120['x']/ 1e2, data_120['phi'][-1,0,:], label=r'120 pts$^*$', color=colors[3])
# ax2.plot(data_80['x'] / 1e2, data_80['ni'][-1,0,:] - data_80['ne'][-1,0,:],
#          label='80 pts', color=colors[0])
# ax2.plot(data_160['x'] / 1e2, data_160['ni'][-1,0,:] - data_160['ne'][-1,0,:],
#          label='160 pts', color=colors[1])
# ax2.plot(data_240['x'] / 1e2, data_240['ni'][-1,0,:] - data_240['ne'][-1,0,:],
#          label='240 pts', color=colors[2])
# ax2.plot(data_120['x'] / 1e2, data_120['ni'][-1,0,:] - data_120['ne'][-1,0,:],
#          label=r'120 pts$^*$', color=colors[3])

ax2.legend()
ax0.set_yscale('log')
ax1.set_yscale('log')
ax0.set_ylim([8e15,4e18])
ax1.set_ylim([8e15,4e18])
ax2.set_xlabel(r'Position [$cm$]')
ax0.set_ylabel(r'Electron Density [$m^{-3}$]')
ax1.set_ylabel(r'Ion Density [$m^{-3}$]')
ax2.set_ylabel(r'Potential [$V$]')

ax0.set_title('(a)')
ax1.set_title('(b)')
ax2.set_title('(c)')
gs.tight_layout(fig, rect=[0, 0, 1, 1])

plt.savefig('Figures/gridRef.eps',dpi=300)

plt.figure()
plt.contourf(data_120['x'], data_120['t'], np.log10(data_120['ne'][:,0,:]))

plt.show()
