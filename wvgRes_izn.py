import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
import time
from scipy import integrate
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

cm_subsection = np.linspace(0.0, 1.0, 7)
colors = [ mpl.cm.viridis(x) for x in cm_subsection ]
color2= plt.rcParams['axes.prop_cycle'].by_key()['color']

# Uncomment for LaTex fonts:
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
mpl.rc('text', usetex=True)

t1 = time.time()

data = np.load('Data/2D_1250pulse_100x100_200ns.npz')

x = data['x']
y = data['y']
t = data['t']
Vd = data['Vd']
Id = data['Id']
ne = data['ne']
Te = data['nt']

eps = 8.8542e-12
me = 9.109e-31
e = 1.6022e-19
kb = 1.3806503e-23
Tg = 300.0
p = 3.0
ninf = p * 101325.0 / 760.0 / kb / Tg

Ej = lambda y, z: (np.sin(np.pi * y / 2.0) * np.sin(np.pi * z / 2.0))**2
iEj, err = integrate.dblquad(Ej, 0, 2.0, lambda x: 0, lambda x: 2.0)

def int2D(A):
    return np.trapz(np.trapz(A, axis=1), axis=0)

def int3D(A):
    return np.trapz(np.trapz(np.trapz(A, axis=2), axis=1), axis=0)

def nu(x):
    y = np.log(np.maximum(2.34e-1,np.minimum(1.57e2, x)))

    a = -32.275912575
    b =   1.45173283977
    c =   0.00936933121094
    d =   0.129397015353
    f =  -0.0414865809044
    g =  -0.0582934303409
    h =   0.0309832277826
    i =  -0.00542014733763
    j =   0.000325615321708

    return np.exp(a+b*y + c*y**2 + d*y**3 + f*y**4 + g*y**5 + h*y**6
                  + i*y**7 + j*y**8) * ninf
nx = 60
ny = 80
nz = 80
Mat = np.zeros([nx, ny, nz], dtype='complex')
Ex = np.zeros([ny, nz])

X = np.linspace(0,1,nx+2)
Y = np.linspace(0,2,ny+2)
Z = np.linspace(0,2,nz+2)

dX = ((X[2:] - X[1:-1]) + (X[1:-1] - X[:-2])) / 2.0
dY = ((Y[2:] - Y[1:-1]) + (Y[1:-1] - Y[:-2])) / 2.0
dZ = ((Z[2:] - Z[1:-1]) + (Z[1:-1] - Z[:-2])) / 2.0

wp = np.zeros([len(x), len(y)], dtype='complex')
wr = np.pi * 2.9989e8 * (1.0 / 2e-2**2 + 1.0 / 2e-2**2)**0.5
a0 = wr / 2.0 / 2000.

# nt_mi = 30
# t = np.concatenate([np.linspace(-1,0,nt_mi),t])

tt = 1
jmin = np.argmin((t-0.09)**2)
jmax = np.argmin((t-0.2)**2)
nt = len(t[jmin:jmax:tt])
dw = np.zeros(nt, dtype='complex')
Q = np.zeros(nt)
Qa = 2000.

wT0 = wr * (((1 - (1j + 1.0) / Qa + 0)**0.5 - 1).real + 1)
wT = np.array([wT0, wT0+2*np.pi*50e6, wT0+2*np.pi*75e6, wT0+2*np.pi*100e6,
               wT0+2*np.pi*125e6, wT0+2*np.pi*150e6, wT0+2*np.pi*200e6])
T = np.zeros([nt, len(wT)])

n = 0
for m in np.arange(jmin,jmax,tt):
    val = 0
    print 'Running at t = {:.3f} us and Vd = {:.2f} V'.format(t[m], Vd[m])

    for i in range(len(x)):
        for j in range(len(y)):
            # Plasm freq Term: wp**2 / w**2 / (w**2 + nu**2)
            wp[i,j] = e**2 * ne[m,j,i] / me / eps / (wr**2 + 1j * wr * nu(Te[m,j,i]))

    for j in range(len(dY)):
        for k in range(len(dZ)):
            Ex[j,k] = Ej(Y[j], Z[k])

    print 'Assembling Matrix...'
    for i in range(len(dX)):
        for j in range(len(dY)):
            for k in range(len(dZ)):
                r = np.sqrt((Y[j+1] - 1.0)**2 + (Z[k+1] - 1.0)**2)
                rj = (np.abs(y*1000 - r)).argmin()
                Mat[i,j,k] = wp[i,rj] * Ex[j,k] * dX[i] * dY[j] * dZ[k]

    print 'Integrating Matrix...'
    val = int3D(Mat)

    dw[n] = ((1 - (1j + 1.0) / Qa + val)**0.5 - 1)
    w = wr * (dw[n].real + 1)
    Q[n]= -(dw[n].real + 1.0) / dw[n].imag / 2.0

    for l in range(len(wT)):
        T[n,l] = 1.0 / Qa**2 / ((wT[l]/w - w/wT[l])**2 + (1.0/Q[n])**2)

    dw[n] = wr * dw[n].real / np.pi / 2.0 / 1e6
    print 'Result:  dw = {:.2f} MHz, Q = {:.2f}\n'.format(dw[n],Q[n])

    n = n+1

t2 = time.time()
print 'Elapsed Time:  {} sec'.format(int(t2 - t1))

fig = plt.figure(figsize=(5.25,3))
gs = gridspec.GridSpec(1,2,hspace=0)
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

t2 = np.linspace(t[jmin], t[jmax-1], 500)
for l in range(len(wT)-1,-1,-1):
    T2 = spline(t[jmin:jmax:tt], 10*np.log10(T[:,l]), t2)
    ax0.plot(t2, T2, color=colors[l],
             label=r'{}'.format(int((wT[l]-wr)/1e6/2.0/np.pi)))
ax0.set_xlabel(r'Time [$\mu$s]')
ax0.set_ylabel(r'Transmission [$dB$]')
# ax0.set_xticks(np.arange(0,22,5))
# ax0.set_xlim([-2,22])

gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig('Figures/2dpulseT_inz.eps', dpi=300)
plt.show()
