import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys

def plots():
    root = os.getcwd() # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    maxOrbits = 500
    ivar = 0 
    while (ivar <= maxOrbits):
        plotRuns(dir_run_list,ivar,root) # run correct python file and save plots locally
        moveImages(dir_run_list,ivar,root) # copy over images into PencilAnalysis
        ivar = ivar + 50
    os.chdir(root)    

def plotRuns(dir_run_list,ivar,root):
    print("============")
    print("entering plotRuns")
    for i in range(len(dir_run_list)):
        print("dir_run_list is " + dir_run_list[i])
        os.chdir(dir_run_list[i])
        print("current cwd is " + os.path.split(os.getcwd())[1])
        try:
            plotRun(ivar)
            os.chdir(root)
        except:
            print("============")
            print("ignoring run")
            print("============")
            os.chdir(root)
    print("leaving plotRuns")
    print("============") 

def moveImages(dir_run_list,ivar,root):
     for i in range(len(dir_run_list)):
         os.chdir(dir_run_list[i])
         name = os.path.split(os.getcwd())[1]
         fullName = str(name) + "-ivar-" + str(ivar) + ".png"
         os.chdir(root)
         os.system('cp -rf ' +  dir_run_list[i] + "/" + str(fullName) + " ./" + "Pencil_Analysis")
         os.chdir( root)


def plotRun(ivar):
    name = os.path.split(os.getcwd())[1] 
    ff=pc.read_var(trimall=True,ivar=ivar,magic=["TT"])
    ff0=pc.read_var(trimall=True,ivar=0,magic=["TT"])
    dfT=ff.TT[:]-ff0.TT[:]
    rad = ff.x
    theta = ff.y
    rad2d, theta2d = np.meshgrid(rad, theta)
    x2d = rad2d*np.cos(theta2d)
    y2d = rad2d*np.sin(theta2d)
    plt.plot(ff.x,np.log(np.sum(dfT**2,axis=0)))
    #plt.contourf(x2d,y2d,dfT,256)
    #
    plt.grid(True)
    Dir_gamma_patches = mpatches.Patch(
        color='white', label=r'$\gamma$ :' + str(np.mean(np.log(np.sum(dfT ** 2, axis=0)))))
    plt.legend(handles=[Dir_gamma_patches], loc=2)
    plt.savefig(name+"-log-" + str(ivar)+".png")
    #plt.savefig(name+"-ivar-" + str(ivar)+".png")
    plt.close()

plots()
