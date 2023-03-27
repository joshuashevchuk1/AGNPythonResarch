import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys


Orbits = int(input('Orbits?:'))

f0=pc.read_var(trimall=True,ivar=Orbits,magic=['TT'])
ff = pc.read_var(trimall=True, ivar=Orbits,magic=['TT'])

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

print('ff.rho : ',ff.rho)
print('ff.TT : ',ff.TT)

T0 = np.mean(np.transpose(ff.T0),axis=1)
TT = np.mean(np.transpose(ff.TT),axis=1)

plt.plot(rad,T0)
plt.plot(rad,TT)
plt.show()

