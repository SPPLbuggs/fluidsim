import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

y, x = np.meshgrid(np.linspace(-10, 10,100), np.linspace(-10, 10,100))

z = np.sin(x)*np.sin(x)+np.sin(y)*np.sin(y)

v = np.linspace(-10, 10,100)
t = np.sin(v)*np.sin(v)
tt = np.cos(v)*np.cos(v)
###########

fig = plt.figure(figsize=(4.5, 3),facecolor='white')
gs = gridspec.GridSpec(1, 1)

#############################
ax2 = fig.add_subplot(gs[0,0])
quad1 = ax2.pcolormesh(x,y,z,shading='gouraud')
ax2.set_xlabel('time')
ax2.set_ylabel('amplitude')
cb2 = fig.colorbar(quad1,ax=ax2)

def init():
    quad1.set_array([])
    return quad1

def animate(iter):
    t = np.sin(2*v-iter/(2*np.pi))*np.sin(2*v-iter/(2*np.pi))
    tt = np.cos(2*v-iter/(2*np.pi))*np.cos(2*v-iter/(2*np.pi))
    z = np.sin(x-iter/(2*np.pi))*np.sin(x-iter/(2*np.pi))+np.sin(y)*np.sin(y)
    quad1.set_array(z.ravel())
    return quad1

gs.tight_layout(fig)

anim = animation.FuncAnimation(fig,animate,frames=100,interval=50,blit=False,repeat=False)
plt.show()

print 'Finished!!'
