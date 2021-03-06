import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colorbar as clb
import matplotlib.gridspec as gridspec
import sys
import time
from scipy import integrate
from scipy.interpolate import spline

size = 14
med_size = 15
big_size = 16

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

t1 = time.time()
# data = np.load('Data/2D_350AC_nx80_ny80_80us.npz')
#
# x = data['x']
# y = data['y']
# t = data['t']
# Vd = data['Vd']
# Id = data['Id']
# ne = data['ne']
# Te = data['nt']

eps = 8.8542e-12
me = 9.109e-31
e = 1.6022e-19
kb = 1.3806503e-23
Tg = 300.0
p = 3.0
ninf = p * 101325.0 / 760.0 / kb / Tg

# ttemp = np.linspace(t[-100], t[-1], 1000)
# vtemp = spline(t[-100:], Vd[-100:], ttemp)
# itemp = spline(t[-100:], Id[-100:], ttemp)
# Vrms = np.sqrt( np.sum(vtemp**2) / 1000. )
# Irms = np.sqrt( np.sum(itemp**2) / 1000. )

def mu(x):
    y = np.log(np.maximum(2.34e-1,np.minimum(1.57e2, x)))

    a1 =  58.1133016145
    b1 =  -0.984082217962
    c1 =  -0.164770900119
    d1 =  -0.142058042584
    f1 =  -0.0637079234081
    g1 =   0.0436642742558

    a2 =  70.7846754548
    b2 = -22.7558138237
    c2 =  12.8366493242
    d2 =  -3.57242763244
    f2 =   0.486276623664
    g2 =  -0.0259045697422

    soln = np.zeros(len(x))
    for i in range(len(x)):
        if y[i] < np.log(5.119):
            soln[i] = np.exp(a1 + b1*y[i] + c1*y[i]**2 + d1*y[i]**3
                             + f1*y[i]**4 + g1*y[i]**5)
        else:
            soln[i] = np.exp(a2 + b2*y[i] + c2*y[i]**2 + d2*y[i]**3
                             + f2*y[i]**4 + g2*y[i]**5)

    return soln

# mue = mu([5]) / ninf
# ne_rms = Irms * 1e-2 / (e * mue * Vrms * np.pi * (1.5e-3)**2)
#
# print Vrms, Irms, ne_rms, mue


xl = 7.5e-3
yl = 15e-3
zl = 15e-3

Ej = lambda y, z: 4.0 / (xl * yl * zl) * (np.sin(np.pi * y / yl) * np.sin(np.pi * z / zl))**2

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

data = np.load('Data/4Torr/350V/60x60_58us.npz')

x = data['x']
y = data['y']
t = data['t']
Vd = data['Vd']
Id = data['Id']
ne = data['ne']
Te = data['nt']

nx = 60
ny = 60
nz = 60
Mat = np.zeros([nx, ny, nz], dtype='complex')
Ex = np.zeros([ny, nz])

X = np.linspace(0,xl,nx+2)
Y = np.linspace(0,yl,ny+2)
Z = np.linspace(0,zl,nz+2)

dX = ((X[2:] - X[1:-1]) + (X[1:-1] - X[:-2])) / 2.0
dY = ((Y[2:] - Y[1:-1]) + (Y[1:-1] - Y[:-2])) / 2.0
dZ = ((Z[2:] - Z[1:-1]) + (Z[1:-1] - Z[:-2])) / 2.0

wp = np.zeros([len(x), len(y)], dtype='complex')
wr = np.pi * 2.9989e8 * (1.0 / yl**2 + 1.0 / zl**2)**0.5

nt = 1000
a0 = wr / 2.0 / 2000.
T0 = np.zeros([nt])
T = np.zeros([nt])
wT = np.linspace(13,16,nt)*2.0*np.pi*1e9
print 'Res Freq = {:.2f} GHz'.format(wr / np.pi / 2.0 / 1e9)

m = -1
for i in range(len(x)):
    for j in range(len(y)):
        # Plasm freq Term: wp**2 / (wr**2 + 1j * wr * nu)
        wp[i,j] = e**2 * ne[m,j,i] / me / eps / (wr**2 + 1j * wr * nu(Te[m,j,i]))

# plt.plot(x, ne[-1,0,:])
# plt.yscale('log')
# sys.exit()

for j in range(len(dY)):
    for k in range(len(dZ)):
        Ex[j,k] = Ej(Y[j], Z[k])

print 'Assembling Matrix...'
for i in range(len(dX)):
    xi = (np.abs(x - X[i])).argmin()
    for j in range(len(dY)):
        for k in range(len(dZ)):
            r = np.sqrt((Y[j+1] - yl/2.0)**2 + (Z[k+1] - zl/2.0)**2)
            rj = (np.abs(y - r)).argmin()
            Mat[i,j,k] = wp[xi,rj] * Ex[j,k] * dX[i] * dY[j] * dZ[k]

print 'Integrating Matrix...'
val = int3D(Mat)

Qa = 2000.
dw = ((1 - (1j + 1.0) / Qa + val)**0.5 - 1)
w = wr * (dw.real + 1)
Q = -(dw.real + 1.0) / dw.imag / 2.0

for l in range(nt):
    T[l] = 1.0 / Qa**2 / ((wT[l]/w - w/wT[l])**2 + (1.0/Q)**2)

    T0[l] = 1.0 / Qa**2 / ((wT[l]/wr - wr/wT[l])**2 + (1.0/Qa)**2)

dw = wr * dw / np.pi / 2.0 / 1e6
print 'Result:  dw = {:.2f} MHz, Q = {:.2f}\n'.format(dw,Q)

fig = plt.figure()

gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0,0])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(wT / 2.0 / np.pi / 1e9, 10*np.log10(T0), label='0 V')
ax.plot(wT / 2.0 / np.pi / 1e9, 10*np.log10(T), label='350 V')


### Second Plot ###
data = np.load('Data/4Torr/400V/60x60_48us.npz')

x = data['x']
y = data['y']
t = data['t']
Vd = data['Vd']
Id = data['Id']
ne = data['ne']
Te = data['nt']

nx = 60
ny = 60
nz = 60
Mat = np.zeros([nx, ny, nz], dtype='complex')
Ex = np.zeros([ny, nz])

X = np.linspace(0,xl,nx+2)
Y = np.linspace(0,yl,ny+2)
Z = np.linspace(0,zl,nz+2)

dX = ((X[2:] - X[1:-1]) + (X[1:-1] - X[:-2])) / 2.0
dY = ((Y[2:] - Y[1:-1]) + (Y[1:-1] - Y[:-2])) / 2.0
dZ = ((Z[2:] - Z[1:-1]) + (Z[1:-1] - Z[:-2])) / 2.0

wp = np.zeros([len(x), len(y)], dtype='complex')
wr = np.pi * 2.9989e8 * (1.0 / yl**2 + 1.0 / zl**2)**0.5

nt = 1000
a0 = wr / 2.0 / 2000.
T0 = np.zeros([nt])
T = np.zeros([nt])
wT = np.linspace(13,16,nt)*2.0*np.pi*1e9
print 'Res Freq = {:.2f} GHz'.format(wr / np.pi / 2.0 / 1e9)

m = -1
for i in range(len(x)):
    for j in range(len(y)):
        # Plasm freq Term: wp**2 / (wr**2 + 1j * wr * nu)
        wp[i,j] = e**2 * ne[m,j,i] / me / eps / (wr**2 + 1j * wr * nu(Te[m,j,i]))

for j in range(len(dY)):
    for k in range(len(dZ)):
        Ex[j,k] = Ej(Y[j], Z[k])

print 'Assembling Matrix...'
for i in range(len(dX)):
    xi = (np.abs(x - X[i])).argmin()
    for j in range(len(dY)):
        for k in range(len(dZ)):
            r = np.sqrt((Y[j+1] - yl/2.0)**2 + (Z[k+1] - zl/2.0)**2)
            rj = (np.abs(y - r)).argmin()
            Mat[i,j,k] = wp[xi,rj] * Ex[j,k] * dX[i] * dY[j] * dZ[k]

print 'Integrating Matrix...'
val = int3D(Mat)

Qa = 2000.
dw = ((1 - (1j + 1.0) / Qa + val)**0.5 - 1)
w = wr * (dw.real + 1)
Q = -(dw.real + 1.0) / dw.imag / 2.0

for l in range(nt):
    T[l] = 1.0 / Qa**2 / ((wT[l]/w - w/wT[l])**2 + (1.0/Q)**2)

    T0[l] = 1.0 / Qa**2 / ((wT[l]/wr - wr/wT[l])**2 + (1.0/Qa)**2)

dw = wr * dw / np.pi / 2.0 / 1e6
print 'Result:  dw = {:.2f} MHz, Q = {:.2f}\n'.format(dw,Q)


ax.plot(wT / 2.0 / np.pi / 1e9, 10*np.log10(T), label='400 V')

t2 = time.time()
print 'Elapsed Time:  {} sec'.format(int(t2 - t1))

plt.legend(frameon=False)
plt.xlabel('Frequency [$GHz$]')
plt.ylabel('Transmission [$dB$]')
# plt.ylim([-45,5])
# plt.xlim([9.8,11.7])

gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig('Figures/2DAC_T.eps',dpi=300)
plt.show()
