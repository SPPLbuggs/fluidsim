import numpy as np
import matplotlib.pyplot as plt

eps0   = 8.85418782e-12
mu0    = 4 * np.pi * 1e-7
c0     = 1.0 / np.sqrt(eps0 * mu0)
e      = 1.60217646e-19
kb     = 1.3806503e-23

x0   = 1e-4
ph0  = e / (eps0 * x0)
t0   = 1e-6

Tg     = 300
p      = 2
ninf   = p * 101325.0 / 760.0 / kb / Tg * x0**3

def getDt(T):
    a = 57.7214952841
    b = -0.353731464348
    c = -0.00505156731795
    d =  0.161173720584
    f = -0.196610345869
    g = -0.0719115643218
    h =  0.11874085868
    i = -0.00669712724784
    j = -0.0236445308504
    k =  0.00917326671764
    l = -0.00135453096483
    m =  7.26379684461e-05

    x = np.log(np.maximum(0.234, np.minimum(157.0, T * ph0)))

    return np.exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5
             + h*x**6 + i*x**7 + j*x**8 + k*x**9
             + l*x**10 + m*x**11) #* x0 / ninf * t0

def getDe(T):
    a = 57.7067018813
    b = -0.0699892875381
    c = -0.117645585949
    d =  0.105390295278
    f = -0.102862612604
    g = -0.0469171521686
    h =  0.0584908312121
    i =  0.000578601715687
    j = -0.0122860884883
    k =  0.00426793748856
    l = -0.000590000082557
    m = 3.00706533201e-05

    x = np.log(np.maximum(0.234, np.minimum(157.0, T * ph0)))

    return np.exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5
             + h*x**6 + i*x**7 + j*x**8 + k*x**9
             + l*x**10 + m*x**11) #* x0 / ninf * t0

T = np.logspace(-1, 2.5, 200)
Dt = getDt(T/ph0)
De = getDe(T/ph0)

plt.plot(T, Dt, label='Dt')
plt.plot(T, De, label='De')
plt.plot(T, De*5.0/3.0, label='5/3 De')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()
