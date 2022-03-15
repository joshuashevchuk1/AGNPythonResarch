import pencil_old as pc
import numpy as np
import matplotlib.pyplot as plt

#
# 12/7/2019
#
# This script plots the toomre Q at every grid cell for
# a selected orbit as a contour plot
# runs in python3 or python. Uses matplotlib for graphing
# this script was created during the update to pencil imports in python
#

pi = np.pi
ivar = int(input('Orbit? :'))


def get_params():
    # return all params needed for this script
    par = pc.read_param()
    h = par.cs0
    cs = par.cs0  # sound speed
    return cs


def make_grid(ivar):
    # takes ivar(orbit from dsnap)
    # make grid elements
    # uncomment below if the disk is not isothermal
    # ff=pc.read_var(trimall=True,ivar=ivar,magic=['TT'])
    ff = pc.read_var(trimall=True, ivar=ivar)
    rad = ff.x  # disk grid points in r
    theta = ff.y  # disk grid points in theta
    ux = ff.ux  # disk grid points in vr
    uy = ff.uy  # disk grid points in vtheta
    rho = ff.rho  # disk surface density
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d*np.cos(theta2d)
    y2d = rad2d*np.sin(theta2d)
    return x2d, y2d, uy, rad, theta, rho


def Calc_ToomreContour(vth, rad, phi, rho, cs):
    # takes grid elements,radial velocity,angular velocity,
    # sound speed, and density
    # return an array for Toomre Q contour to be plotted
    # TC_Array MUST be 2d in order to be plotted for a contour
    # intialize 2d array to be plotted for contour

    # set grav_const to what is set in start.in
    # for now, manually set this

    grav_const = 1.3e-4

    x1, x2 = 64, 192
    TC_Array = [[0 for x in range(x1)] for y in range(x2)]
    # redefine each element of TC_array with values for toomre Q
    j = 0
    i = 0
    di = 1
    dj = 1
    while j <= len(phi)-1:
        while i <= len(rad)-1:
            # append toomre Q value at that point
            TC_Array[j][i] = ((vth[j][i])*(cs**0.5))/(grav_const*pi*rho[j][i])
            print('working at radius:'+str(i))
            print('working at theta:'+str(j))
            i = i+di
        print('===========================')
        print('moving to next theta')
        print('===========================')
        i = 0
        j = j+dj
        print('===========================')
        print('done with TC_Array')
        print('===========================')
    return TC_Array


cs = get_params()
x2d, y2d, vth, rad, phi, rho = make_grid(ivar)
TC_Array = Calc_ToomreContour(vth, rad, phi, rho, cs)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_aspect('equal')
# X,Y & data2D must all be the same dimesions
ncolors = 256
PL1 = ax.contourf(x2d, y2d, TC_Array, ncolors)

plt.title('orbit = '+str(ivar))
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
cax.set_aspect(20)
cax.set_xlabel('Q', fontsize=10)
plt.colorbar(PL1, cax=cax)
plt.savefig('QContour.png')
