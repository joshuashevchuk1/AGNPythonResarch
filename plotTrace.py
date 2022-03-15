import pencil as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys

ivar = int(input("Orbit:"))
ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"])
ts = pc.read_ts()

f0 = pc.read_var(trimall=True, ivar=0, magic=["TT"])

print np.amax(ff.rho), np.amin(ff.rho)
print ff.t

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

radq2 = ts.xq2
thetaq2 = ts.yq2

xrq2 = radq2*np.cos(thetaq2)
yrq2 = radq2*np.sin(thetaq2)

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)


# Cartesian plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
# X,Y & data2D must all be the same dimesions
ncolors = 256
ax.contourf(x2d, y2d, ff.TT[:]/f0.TT[:],
            ncolors), ax.plot(xrq2, yrq2, color='orange')
# Normalized_Temp=temp_fv[:]/Init_Temp[:]


plt.show()
