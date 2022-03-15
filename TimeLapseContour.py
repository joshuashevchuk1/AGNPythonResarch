import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys

ff = pc.read_var(trimall=True, ivar=0,magic=["TT"])
f25 = pc.read_var(trimall=True, ivar=25,magic=["TT"])
f50 = pc.read_var(trimall=True, ivar =50,magic=["TT"])
f75 = pc.read_var(trimall=True, ivar =75,magic=["TT"])

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

# Cartesian plot
fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4, figsize=(20, 5))
# X,Y & data2D must all be the same dimesions
ncolors = 256
fig.subplots_adjust(left=0.05, top=0.9, wspace=0.10)

PL1 = ax1.contourf(x2d,y2d,ff.TT,ncolors)
PL2 = ax2.contourf(x2d, y2d, f25.TT-ff.TT, ncolors)
PL3 = ax3.contourf(x2d, y2d, f50.TT-ff.TT, ncolors)
PL4 = ax4.contourf(x2d, y2d, f75.TT-ff.TT, ncolors)

ax1.set_xlabel('Orbit = 0')
ax2.set_xlabel('Orbit = 25')
ax3.set_xlabel('Orbit = 50')
ax4.set_xlabel('Orbit = 75')

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

cax = plt.axes([0.85, 0.1, 0.075, 0.8])
cax.set_aspect(20)
cax.set_ylabel('Temp in code units', fontsize=10)

plt.colorbar(PL4, cax=cax)
#plt.show()
plt.savefig('Contour.png')
