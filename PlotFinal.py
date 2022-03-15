import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
import os
from pylab import *
import sys

ivar = int(input("oribit: "))

name = os.path.split(os.getcwd())[1] 

print("name is " + str(name))

ff=pc.read_var(trimall=True,ivar=ivar,magic=["TT"])
ff0=pc.read_var(trimall=True,ivar=0,magic=["TT"])

dfT=ff.TT[:]-ff0.TT[:]

rad = ff.x
theta = ff.y
rad2d, theta2d = np.meshgrid(rad, theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

#plt.plot(ff.x,np.log(np.sum(dfT**2,axis=0)))
plt.contourf(x2d,y2d,dfT,256)
#
print(str(name)+ "-ivar-" + str(ivar))
plt.savefig(name+"-ivar-" + str(ivar)+".png")
