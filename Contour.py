import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys

#n=int(input("Eccentricty value"))
#dirint=int(input("Run int?"))

Orbits = int(input('Orbits?:'))

f0=pc.read_var(trimall=True,ivar=Orbits,magic=['TT'])
ff = pc.read_var(trimall=True, ivar=Orbits,magic=['TT'])
#print np.amax(ff.rho), np.amin(ff.rho)
#print ff.t

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

# Cartesian plot
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(30, 5))
# plt.suptitle('AGNRun'+str(dirint).zfill(0)+'_EnergyE0.'+str(n).zfill(0))
# plt.tight_layout()
# X,Y & data2D must all be the same dimesions
ncolors = 256
fig.subplots_adjust(left=0.05, top=0.9, wspace=0.10)

PL1 = ax1.contourf(x2d, y2d, ff.rho-f0.rho, ncolors)
ax1.set_xlabel('radial density')
PL2 = ax2.contourf(x2d, y2d, ff.TT-f0.TT, ncolors)
ax2.set_xlabel('radial temperature')

cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# cax.set_label('Kelvin')
cax.set_aspect(20)
cax.set_ylabel('Temp in code units', fontsize=10)

plt.colorbar(PL2, cax=cax)
plt.show()
#plt.savefig('Contour.png')
