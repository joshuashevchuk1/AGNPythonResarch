import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys

# plot denisty, shock and temperature for the local orbit

def plot():
    ivar_array = [0,5,10]
    for i in range(len(ivar_array)):
        plotRuns(ivar)

def plotRuns(ivar):
    plotRunTemperature(ivar)
    plotRunDensity(ivar)
    plotRunShock(ivar)

def plotRunShock(ivar):
    name = os.path.split(os.getcwd())[1]
    ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"], quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"], quiet=True)
    dfs = ff.shock[:] - ff0.shock[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d * np.cos(theta2d)
    y2d = rad2d * np.sin(theta2d)
    plt.contourf(x2d, y2d, dfs, 256)
    plt.savefig(name + "-density-ivar-" + str(ivar) + ".png")
    plt.close()

def plotRunDensity(ivar):
    name = os.path.split(os.getcwd())[1]
    ff=pc.read_var(trimall=True,ivar=ivar,magic=["TT"],quiet=True)
    ff0=pc.read_var(trimall=True,ivar=0,magic=["TT"],quiet=True)
    dfrho=ff.rho[:]-ff0.rho[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d*np.cos(theta2d)
    y2d = rad2d*np.sin(theta2d)
    plt.contourf(x2d,y2d,dfrho,256)
    plt.legend(handles=[Dir_gamma_patches], loc=2)
    plt.savefig(name+"-density-ivar-" + str(ivar)+".png")
    plt.close()

def plotRunTemperature(ivar):
    name = os.path.split(os.getcwd())[1]
    ff=pc.read_var(trimall=True,ivar=ivar,magic=["TT"],quiet=True)
    ff0=pc.read_var(trimall=True,ivar=0,magic=["TT"],quiet=True)
    dfT=ff.TT[:]-ff0.TT[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d*np.cos(theta2d)
    y2d = rad2d*np.sin(theta2d)
    plt.contourf(x2d,y2d,dfT,256)
    plt.savefig(name+"-temperature-ivar-" + str(ivar)+".png")
    plt.close()