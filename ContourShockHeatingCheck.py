import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys


Orbits = int(input('Orbits?:'))

f0=pc.read_var(trimall=True,ivar=0,magic=['TT'])
ff = pc.read_var(trimall=True, ivar=Orbits,magic=['TT'])

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

T0 = np.mean(np.transpose(f0.TT),axis=1)
TT = np.mean(np.transpose(ff.TT),axis=1)

print('T0 : ',T0)
print('TT : ',TT)

plt.plot(rad,T0,label='T0')
plt.plot(rad,TT,label='TT)
plt.show()

