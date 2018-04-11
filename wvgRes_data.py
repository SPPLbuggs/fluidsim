import numpy as np
import sys
import time
from scipy import integrate
from scipy import interpolate

t1 = time.time()

data = np.load('Data/data.npz')
# data = np.load('Data/2Torr/1600V/60x60_40us.npz')

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

xl = 7.9e-3
yl = 13.6e-3
zl = 15.8e-3

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
nx = 40
ny = 40
nz = 40
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
print 'Res Freq = {:.2f} GHz'.format(wr / np.pi / 2.0 / 1e9)

nt = 40
# tgt = np.logspace(np.log10(3e-2), np.log10(t[-1]), nt)
tgt = np.logspace(np.log10(1e-1), np.log10(20.), nt)
ii = np.unique([np.argmax(t - i > 0) for i in tgt])
nt = len(ii)

dw = np.zeros(nt, dtype='complex')
Qa = 1912.

for j in range(len(dY)):
    for k in range(len(dZ)):
        Ex[j,k] = Ej(Y[j], Z[k])

for n in range(nt):
    m = ii[n]
    print 'Running at t = {:.2f} us and Vd = {:.2f} V'.format(t[m], Vd[m])

    for i in range(len(x)):
        for j in range(len(y)):
            # Plasm freq Term: wp**2 / w**2 / (w**2 + nu**2)
            wp[i,j] = e**2 * ne[m,j,i] / me / eps / (wr**2 + 1j * wr * nu(Te[m,j,i]) * 2.0 * np.pi)

    print 'Assembling Matrix...'
    for i in range(len(dX)):
        # xi = (np.abs(x - X[i])).argmin()
        xi = np.argmax(x - X[i] > 0) - 1

        for j in range(len(dY)):
            for k in range(len(dZ)):
                Rtgt = np.sqrt((Y[j+1] - yl/2.0)**2 + (Z[k+1] - zl/2.0)**2)
                rj = np.argmax(y - Rtgt > 0) - 1
                Wp1 = wp[xi,rj] + (Rtgt - y[rj]) * (wp[xi,rj+1] - wp[xi,rj]) / (y[rj+1] - y[rj])
                Wp2 = wp[xi+1,rj] + (Rtgt - y[rj]) * (wp[xi+1,rj+1] - wp[xi+1,rj]) / (y[rj+1] - y[rj])
                Wp = Wp1 + (X[i] - x[xi]) * (Wp2 - Wp1) / (x[xi+1] - x[xi])
                Mat[i,j,k] = Wp * Ex[j,k] * dX[i] * dY[j] * dZ[k]
                # rj = np.argmin((R - r)**2)
                # rj = (np.abs(y - r)).argmin()
                # Mat[i,j,k] = wp[xi,rj] * Ex[j,k] * dX[i] * dY[j] * dZ[k]
                # Mat[i,j,k] = (splineR(X[i], R) + splineI(X[i], R)* 1j) * Ex[j,k] * dX[i] * dY[j] * dZ[k]

    print 'Integrating Matrix...'
    val = int3D(Mat)

    dw[n] = ((1 - (1j + 1.0) / Qa + val)**0.5 - 1)
    print 'Result:  dw = {:.2f} MHz\n'.format(dw[n] * wr / np.pi / 2.0 / 1e6)

np.savez('Data/wvgResData.npz', dw=dw, t=t[ii])
t2 = time.time()
print 'Elapsed Time:  {} sec'.format(int(t2 - t1))
