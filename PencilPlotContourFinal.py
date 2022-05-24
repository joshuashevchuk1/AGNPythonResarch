import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys


# plot denisty, shock and temperature for the local orbit

def plot():
    root = os.getcwd() # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    print("============")
    print("entering plotRuns")
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            plotRun()
            os.chdir(root)
        except:
            print("============")
            print("ignoring run")
            print("============")
            os.chdir(root)
    print("leaving plotRuns")
    print("============")

def plotRun():
    ivar = [1, 50, 100]
    for i in range(len(ivar)):
        plotRuns(ivar[i], False)
        if i == 100:
            plotRuns(ivar[i], True)


def plotRuns(ivar, color):
    plotRunTemperature(ivar, color)
    plotRunDensity(ivar, color)


def plotRunDensity(ivar, color):
    name = os.path.split(os.getcwd())[1]
    ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"], quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"], quiet=True)
    dfrho = ff.rho[:] - ff0.rho[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d * np.cos(theta2d)
    y2d = rad2d * np.sin(theta2d)
    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(bottom=0.07, top=0.95)
    PL2 = ax1.contourf(x2d, y2d, dfrho, 256)
    ax1.set_aspect('equal')

    if ivar == 100:
        cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        cax.set_aspect(20)
        cax.set_ylabel('Density in code units', fontsize=10)
        plt.colorbar(PL2, cax=cax)

    plt.suptitle('ivar= ' + str(ivar))
    plt.savefig(name + "-density-ivar-" + str(ivar) + ".png")
    plt.close()


def plotRunTemperature(ivar, color):
    name = os.path.split(os.getcwd())[1]
    ff = pc.read_var(trimall=True, ivar=ivar, magic=["TT"], quiet=True)
    ff0 = pc.read_var(trimall=True, ivar=0, magic=["TT"], quiet=True)
    dfT = ff.TT[:] - ff0.TT[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d * np.cos(theta2d)
    y2d = rad2d * np.sin(theta2d)
    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(bottom=0.07, top=0.95)
    PL2 = ax1.contourf(x2d, y2d, dfT, 256)
    ax1.set_aspect('equal')

    if ivar == 100:
        cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        cax.set_aspect(20)
        cax.set_ylabel('temperature in code units', fontsize=10)
        plt.colorbar(PL2, cax=cax)

    plt.suptitle('ivar= ' + str(ivar))
    plt.savefig(name + "-temperature-ivar-" + str(ivar) + ".png")
    plt.close()


plot()
