import numpy as np
import matplotlib.pyplot as plt

A = 1.5
r0 = 0.3
w = 0.2

x = np.linspace(0,1,30)
y = A * np.exp(-(x-r0)**2 / w)

A2 = 1.0
r02 = 0.0
w2 = 0.5

for i in range(100):
    for j in range(len(x)):
        temp1 = np.exp(-(x[j] - r02)**2 / w2)
        temp2 = -2.0 * (y[j] - A2 * temp1)
        A2 = A2 - 0.01 * temp2 * temp1
        r02 = r02 - 0.01 * temp2 * temp1 * 2.0 * A2 * (x[j] - r02) / w2
        w2 = w2 - 0.01 * temp2 * temp1 * A2 * (x[j] - r02)**2 / w2**2

print A2, r02, w2

# plt.plot(x,y,'.')
# plt.show()
