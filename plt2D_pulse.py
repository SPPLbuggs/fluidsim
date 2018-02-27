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

data = np.load('Data/2D_1250VPulse_160x160_17us.npz')

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

tgt = [0.111, 0.135, 0.18, 0.501, 17.01]
tgt = [0.18, 0.19, 0.2, 0.21, 0.22]
ts = np.zeros(len(tgt), dtype='int')
j = 0
for i in range(len(t)-1):
    if t[i+1] > tgt[j]:
        ts[j] = i
        j = j+1

    if j == len(tgt):
        break

j = 2

plt.figure()
for i in range(len(tgt)):
    plt.plot(x, ne[ts[i], 0, :])#, color=colors[i])

plt.yscale('log')
plt.show()


# fig = plt.figure(figsize=(5.25,5.25))
# gs = gridspec.GridSpec(2, 1)#, width_ratios=[1.0,0.05])#, height_ratios=[1.0,0.75])
# gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# ax0 = fig.add_subplot(gs1[0,0])
# ax1 = fig.add_subplot(gs2[0,0])
# axc = fig.add_subplot(gs1[0,1])
# ax2 = ax1.twinx()
#
# ax1.spines['top'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# plt.setp(ax0.get_xticklabels(), visible=False)
#
# v = np.linspace(13.9,17.2,25)
# im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
# # im = ax0.pcolormesh(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]),vmin=v[0], vmax=v[-1], shading='gouraud')
# tck = [14.0, 15.0, 16.0, 17.0, 18.0]
# tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
# cb = plt.colorbar(im, cax=axc, ticks=tck)
# cb.outline.set_visible(False)
# cb.ax.set_yticklabels(tck_lbls)
# cb.ax.set_title(r'n$_\mathrm{e}$ [$m^{-3}$]')
# cb.solids.set_edgecolor("face")
#
# for c in im.collections:
#     c.set_edgecolor("face")
#     c.set_linewidth(0.000000000001)
#
# ax0.set_ylabel(r'R [$mm$]')
# ax0.set_title(r'Electron Density')
#
# xnew = np.linspace(x[0],x[-1],300)
# ne_s = spline(x,ne[ts[j],0,:],xnew)
# ni_s = spline(x,ni[ts[j],0,:],xnew)
# nm_s = spline(x,nm[ts[j],0,:],xnew)
# nt_s = spline(x,nt[ts[j],0,:],xnew)
#
# ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
# ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
# ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
# ax2.plot(xnew, nt_s,label=r'T$_\mathrm{e}$', color=color2[3])
# ax1.set_yscale('log')
# ax1.set_xlabel(r'X [$mm$]')
# ax1.set_ylabel(r'Density [$m^{-3}$]')
# ax2.set_ylabel(r'Temperature [$eV$]')
# ax1.set_ylim([2e13,6.3e17])
# ax2.set_yscale('log')
# ax2.set_ylim([5e-2, 5e1])
# # ax1.legend(bbox_to_anchor=(1.2,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)
#
# ax1.text(1.96+0.76, 3e14, r'n$_\mathrm{e}$', color=color2[0], va='center', size=12)
# ax1.annotate('', xy=(1.96, 1e14), xytext=(1.96+0.73, 3e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[0]))
#
# ax1.text(1.51+.76, 8e15, r'n$_\mathrm{i}$', color=color2[1], va='center', size=12)
# ax1.annotate('', xy=(1.51, 3e15), xytext=(1.51+.73, 8e15),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[1]))
#
# ax1.text(4.5+.76, 2e14, r'n$_\mathrm{m}$', color=color2[2], va='center', size=12)
# ax1.annotate('', xy=(4.5, 5e13), xytext=(4.5+.73, 2e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[2]))
#
# ax1.text(7.7+.76, 1e16, r'T$_\mathrm{e}$', color=color2[3], va='center', size=12)
# ax1.annotate('', xy=(7.7, 4e16), xytext=(7.7+.73, 1e16),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[3]))
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_t5.eps', dpi = 300)

# j = 0
# fig = plt.figure(figsize=(5.25,5.25))
# gs = gridspec.GridSpec(2, 1)#, width_ratios=[1.0,0.05])#, height_ratios=[1.0,0.75])
# gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# ax0 = fig.add_subplot(gs1[0,0])
# ax1 = fig.add_subplot(gs2[0,0])
# axc = fig.add_subplot(gs1[0,1])
# ax2 = ax1.twinx()
#
# ax1.spines['top'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# plt.setp(ax0.get_xticklabels(), visible=False)
#
# v = np.linspace(13.9,17.2,25)
# im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
# # im = ax0.pcolormesh(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]),vmin=v[0], vmax=v[-1], shading='gouraud')
# tck = [14.0, 15.0, 16.0, 17.0, 18.0]
# tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
# cb = plt.colorbar(im, cax=axc, ticks=tck)
# cb.outline.set_visible(False)
# cb.ax.set_yticklabels(tck_lbls)
# cb.ax.set_title(r'n$_\mathrm{e}$ [$m^{-3}$]')
# cb.solids.set_edgecolor("face")
#
# for c in im.collections:
#     c.set_edgecolor("face")
#     c.set_linewidth(0.000000000001)
#
# ax0.set_ylabel(r'R [$mm$]')
# ax0.set_title(r'Electron Density')
#
# xnew = np.linspace(x[0],x[-1],300)
# ne_s = spline(x,ne[ts[j],0,:],xnew)
# ni_s = spline(x,ni[ts[j],0,:],xnew)
# nm_s = spline(x,nm[ts[j],0,:],xnew)
# nt_s = spline(x,nt[ts[j],0,:],xnew)
#
# ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
# ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
# ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
# ax2.plot(xnew, nt_s,label=r'T$_\mathrm{e}$', color=color2[3])
# ax1.set_yscale('log')
# ax1.set_xlabel(r'X [$mm$]')
# ax1.set_ylabel(r'Density [$m^{-3}$]')
# ax2.set_ylabel(r'Temperature [$eV$]')
# ax1.set_ylim([2e13,6.3e17])
# ax2.set_yscale('log')
# ax2.set_ylim([5e-2, 5e1])
# # ax1.legend(bbox_to_anchor=(1.2,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)
#
# ax1.text(1.96+0.76, 3e14, r'n$_\mathrm{e}$', color=color2[0], va='center', size=12)
# ax1.annotate('', xy=(1.96, 1e14), xytext=(1.96+0.73, 3e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[0]))
#
# ax1.text(1.51+.76, 8e15, r'n$_\mathrm{i}$', color=color2[1], va='center', size=12)
# ax1.annotate('', xy=(1.51, 3e15), xytext=(1.51+.73, 8e15),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[1]))
#
# ax1.text(4.5+.76, 2e14, r'n$_\mathrm{m}$', color=color2[2], va='center', size=12)
# ax1.annotate('', xy=(4.5, 5e13), xytext=(4.5+.73, 2e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[2]))
#
# ax1.text(7.7+.76, 1e16, r'T$_\mathrm{e}$', color=color2[3], va='center', size=12)
# ax1.annotate('', xy=(7.7, 4e16), xytext=(7.7+.73, 1e16),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[3]))
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_t1.eps', dpi = 300)
#
# j = 1
# fig = plt.figure(figsize=(5.25,5.25))
# gs = gridspec.GridSpec(2, 1)#, width_ratios=[1.0,0.05])#, height_ratios=[1.0,0.75])
# gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# ax0 = fig.add_subplot(gs1[0,0])
# ax1 = fig.add_subplot(gs2[0,0])
# axc = fig.add_subplot(gs1[0,1])
# ax2 = ax1.twinx()
#
# ax1.spines['top'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# plt.setp(ax0.get_xticklabels(), visible=False)
#
# v = np.linspace(13.9,17.2,25)
# im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
# # im = ax0.pcolormesh(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]),vmin=v[0], vmax=v[-1], shading='gouraud')
# tck = [14.0, 15.0, 16.0, 17.0, 18.0]
# tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
# cb = plt.colorbar(im, cax=axc, ticks=tck)
# cb.outline.set_visible(False)
# cb.ax.set_yticklabels(tck_lbls)
# cb.ax.set_title(r'n$_\mathrm{e}$ [$m^{-3}$]')
# cb.solids.set_edgecolor("face")
#
# for c in im.collections:
#     c.set_edgecolor("face")
#     c.set_linewidth(0.000000000001)
#
# ax0.set_ylabel(r'R [$mm$]')
# ax0.set_title(r'Electron Density')
#
# xnew = np.linspace(x[0],x[-1],300)
# ne_s = spline(x,ne[ts[j],0,:],xnew)
# ni_s = spline(x,ni[ts[j],0,:],xnew)
# nm_s = spline(x,nm[ts[j],0,:],xnew)
# nt_s = spline(x,nt[ts[j],0,:],xnew)
#
# ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
# ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
# ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
# ax2.plot(xnew, nt_s,label=r'T$_\mathrm{e}$', color=color2[3])
# ax1.set_yscale('log')
# ax1.set_xlabel(r'X [$mm$]')
# ax1.set_ylabel(r'Density [$m^{-3}$]')
# ax2.set_ylabel(r'Temperature [$eV$]')
# ax1.set_ylim([2e13,6.3e17])
# ax2.set_yscale('log')
# ax2.set_ylim([5e-2, 5e1])
# # ax1.legend(bbox_to_anchor=(1.2,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)
#
# ax1.text(7.1-0.8, 5e14, r'n$_\mathrm{e}$', color=color2[0], va='center', ha='right', size=12)
# ax1.annotate('', xy=(7.1, 1e15), xytext=(7.1-0.73, 5e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[0]))
#
# ax1.text(1.5+.76, 3e16, r'n$_\mathrm{i}$', color=color2[1], va='center', size=12)
# ax1.annotate('', xy=(1.5, 1.5e17), xytext=(1.5+.73, 3e16),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[1]))
#
# ax1.text(3.9+.76, 5e17, r'n$_\mathrm{m}$', color=color2[2], va='center', size=12)
# ax1.annotate('', xy=(3.9, 2e17), xytext=(3.9+.73, 5e17),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[2]))
#
# ax1.text(8.5+.76, 7e16, r'T$_\mathrm{e}$', color=color2[3], va='center', size=12)
# ax1.annotate('', xy=(8.5, 1.5e17), xytext=(8.5+.73, 7e16),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[3]))
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_t2.eps', dpi = 300)
#
# j = 2
# fig = plt.figure(figsize=(5.25,5.25))
# gs = gridspec.GridSpec(2, 1)#, width_ratios=[1.0,0.05])#, height_ratios=[1.0,0.75])
# gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# ax0 = fig.add_subplot(gs1[0,0])
# ax1 = fig.add_subplot(gs2[0,0])
# axc = fig.add_subplot(gs1[0,1])
# ax2 = ax1.twinx()
#
# ax1.spines['top'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# plt.setp(ax0.get_xticklabels(), visible=False)
#
# v = np.linspace(13.9,17.2,25)
# im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
# # im = ax0.pcolormesh(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]),vmin=v[0], vmax=v[-1], shading='gouraud')
# tck = [14.0, 15.0, 16.0, 17.0, 18.0]
# tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
# cb = plt.colorbar(im, cax=axc, ticks=tck)
# cb.outline.set_visible(False)
# cb.ax.set_yticklabels(tck_lbls)
# cb.ax.set_title(r'n$_\mathrm{e}$ [$m^{-3}$]')
# cb.solids.set_edgecolor("face")
#
# for c in im.collections:
#     c.set_edgecolor("face")
#     c.set_linewidth(0.000000000001)
#
# ax0.set_ylabel(r'R [$mm$]')
# ax0.set_title(r'Electron Density')
#
# xnew = np.linspace(x[0],x[-1],300)
# ne_s = spline(x,ne[ts[j],0,:],xnew)
# ni_s = spline(x,ni[ts[j],0,:],xnew)
# nm_s = spline(x,nm[ts[j],0,:],xnew)
# nt_s = spline(x,nt[ts[j],0,:],xnew)
#
# ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
# ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
# ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
# ax2.plot(xnew, nt_s,label=r'T$_\mathrm{e}$', color=color2[3])
# ax1.set_yscale('log')
# ax1.set_xlabel(r'X [$mm$]')
# ax1.set_ylabel(r'Density [$m^{-3}$]')
# ax2.set_ylabel(r'Temperature [$eV$]')
# ax1.set_ylim([2e13,6.3e17])
# ax2.set_yscale('log')
# ax2.set_ylim([5e-2, 5e1])
# # ax1.legend(bbox_to_anchor=(1.2,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)
#
# ax1.text(0.2+0.76, 2e14, r'n$_\mathrm{e}$', color=color2[0], va='center', size=12)
# ax1.annotate('', xy=(0.2, 5e14), xytext=(0.2+0.73, 2e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[0]))
#
# ax1.text(1.5+.76, 3e16, r'n$_\mathrm{i}$', color=color2[1], va='center', size=12)
# ax1.annotate('', xy=(1.5, 1.2e17), xytext=(1.5+.73, 3e16),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[1]))
#
# ax1.text(3.9+.76, 6e17, r'n$_\mathrm{m}$', color=color2[2], va='center', size=12)
# ax1.annotate('', xy=(3.9, 2.5e17), xytext=(3.9+.73, 6e17),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[2]))
#
# ax1.text(7.6+.76, 5e14, r'T$_\mathrm{e}$', color=color2[3], va='center', size=12)
# ax1.annotate('', xy=(7.6, 2e15), xytext=(7.6+.73, 5e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[3]))
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_t3.eps', dpi = 300)
#
# j = 3
# fig = plt.figure(figsize=(5.25,5.25))
# gs = gridspec.GridSpec(2, 1)#, width_ratios=[1.0,0.05])#, height_ratios=[1.0,0.75])
# gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1,0],
#                                        hspace=0.0, width_ratios=[1.0,0.03])
# ax0 = fig.add_subplot(gs1[0,0])
# ax1 = fig.add_subplot(gs2[0,0])
# axc = fig.add_subplot(gs1[0,1])
# ax2 = ax1.twinx()
#
# ax1.spines['top'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# plt.setp(ax0.get_xticklabels(), visible=False)
#
# v = np.linspace(13.9,17.2,25)
# im = ax0.contourf(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]), v)
# # im = ax0.pcolormesh(x,y,np.maximum(np.log10(ne[ts[j],:,:]), v[0]),vmin=v[0], vmax=v[-1], shading='gouraud')
# tck = [14.0, 15.0, 16.0, 17.0, 18.0]
# tck_lbls = [r'$10^{14}$', r'$10^{15}$', r'$10^{16}$', r'$10^{17}$', r'$10^{18}$']
# cb = plt.colorbar(im, cax=axc, ticks=tck)
# cb.outline.set_visible(False)
# cb.ax.set_yticklabels(tck_lbls)
# cb.ax.set_title(r'n$_\mathrm{e}$ [$m^{-3}$]')
# cb.solids.set_edgecolor("face")
#
# for c in im.collections:
#     c.set_edgecolor("face")
#     c.set_linewidth(0.000000000001)
#
# ax0.set_ylabel(r'R [$mm$]')
# ax0.set_title(r'Electron Density')
#
# xnew = np.linspace(x[0],x[-1],300)
# ne_s = spline(x,ne[ts[j],0,:],xnew)
# ni_s = spline(x,ni[ts[j],0,:],xnew)
# nm_s = spline(x,nm[ts[j],0,:],xnew)
# nt_s = spline(x,nt[ts[j],0,:],xnew)
#
# ax1.plot(xnew, ne_s,label=r'n$_\mathrm{e}$')
# ax1.plot(xnew, ni_s,label=r'n$_\mathrm{i}$')
# ax1.plot(xnew, nm_s,label=r'n$_\mathrm{m}$')
# ax2.plot(xnew, nt_s,label=r'T$_\mathrm{e}$', color=color2[3])
# ax1.set_yscale('log')
# ax1.set_xlabel(r'X [$mm$]')
# ax1.set_ylabel(r'Density [$m^{-3}$]')
# ax2.set_ylabel(r'Temperature [$eV$]')
# ax1.set_ylim([2e13,6.3e17])
# ax2.set_yscale('log')
# ax2.set_ylim([5e-2, 5e1])
# # ax1.legend(bbox_to_anchor=(1.2,0.75), frameon=False)#bbox_to_anchor = (1.3, 0.55), loc = 5, frameon=False)
#
# ax1.text(0.21+0.76, 1e14, r'n$_\mathrm{e}$', color=color2[0], va='center', size=12)
# ax1.annotate('', xy=(0.21, 5e14), xytext=(0.21+0.73, 1e14),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[0]))
#
# ax1.text(1.5+.76, 6e15, r'n$_\mathrm{i}$', color=color2[1], va='center', size=12)
# ax1.annotate('', xy=(1.5, 2.2e16), xytext=(1.5+.73, 6e15),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[1]))
#
# ax1.text(3.9+.76, 6e17, r'n$_\mathrm{m}$', color=color2[2], va='center', size=12)
# ax1.annotate('', xy=(3.9, 2.3e17), xytext=(3.9+.73, 6e17),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.27",
#              color=color2[2]))
#
# ax1.text(7.6+.76, 8e13, r'T$_\mathrm{e}$', color=color2[3], va='center', size=12)
# ax1.annotate('', xy=(7.6, 3e14), xytext=(7.6+.73, 8e13),
#              arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=-0.27",
#              color=color2[3]))
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_t4.eps', dpi = 300)
#
# ## Figure 2:: IV
# ## Current/Voltage Figure ###
# fig = plt.figure(figsize=(5,3))
# gs = gridspec.GridSpec(1,1)
# ax0 = fig.add_subplot(gs[0])
#
# #ax.spines['right'].set_visible(False)
# ax0.spines['top'].set_visible(False)
#
# ax0.plot(t, Vd, color=color2[0], label=r'V$_\mathrm{d}$')
# ax1 = ax0.twinx()
# ax1.plot(t, Id * 1000, color=color2[1], label=r'I$_\mathrm{d}$')
#
# ax0.set_xscale('log')
# ax0.set_xlim([0.008,25])
# ax0.set_ylim([-150,1350])
# # ax0.set_yticks(np.arange(25,126,25))
# # ax0.set_ylim([20,130])
#
# # ax0.tick_params('y', colors=colors[0])
# ax0.set_ylabel(r'Discharge Voltage [$V$]')
# ax0.set_xlabel(r'Time [$\mu s$]')
#
# ax1.set_ylabel(r'Discharge Current [$mA$]')
# # ax1.tick_params('y', colors=colors[2])
# # ax1.set_yscale('log')
# #ax1.set_ylim([200,1000])
# ax0.spines['top'].set_visible(False)
# ax1.spines['top'].set_visible(False)
#
# # plt.annotate(r'(a)', fontsize=12, xy=(0.235,1.0), xycoords='axes fraction')
# # s = np.linspace(-150,1350,10)
# # v = np.ones(10)*0.105
# # ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)
# #
# # plt.annotate(r'(b)', fontsize=12, xy=(0.34,1.0), xycoords='axes fraction')
# # v = np.ones(10)*0.13
# # ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)
# #
# # plt.annotate(r'(c)', fontsize=12, xy=(0.51,1.0), xycoords='axes fraction')
# # v = np.ones(10)*0.5
# # ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)
# #
# # plt.annotate(r'(d)', fontsize=12, xy=(0.9,1.0), xycoords='axes fraction')
# # v = np.ones(10)*30
# # ax0.plot(v,s,'--',color=(0.5,0.5,0.5), linewidth=1)
#
# plt.annotate('', xy=(0.095, 0.42), xycoords='axes fraction',
#              xytext=(0.22, 0.27), textcoords='axes fraction',
#              arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=0.42",
#              color=color2[0]))
#
# plt.annotate('', xy=(0.4+0.15, 0.42), xycoords='axes fraction',
#              xytext=(0.4, 0.23), textcoords='axes fraction',
#              arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.36",
#              color=color2[1]))
#
# ax1.legend(loc=1, bbox_to_anchor=(0.97,0.95), frameon=False)
# ax0.legend(loc=2, bbox_to_anchor=(0.0,0.95), frameon=False)
#
# gs.tight_layout(fig, rect=[0, 0, 1, 1])
# plt.savefig('Figures/2dpulse_IV_t1.eps', dpi = 300)

plt.show()
