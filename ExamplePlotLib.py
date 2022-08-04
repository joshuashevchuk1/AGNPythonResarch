import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json

x=[]
y=[]
z=[]

def plot():
    init()
    plotcontour()

def plotcontour():
    global x,y
    x,y = np.meshgrid(x,y)
    plt.contourf(x,y,z,256)
    plt.show()

def init():
    print('initializing arrays')
    initXarray()
    initYarray()
    initZarray()

def initXarray():
    global x
    x=np.arange(0,5,1)
    print('xarray is ',xarray)

def initYarray():
    global y
    y=np.arange(0,5,1)
    print('yarray is ',yarray)

def initZarray():
    global z
    z=np.arange(0,5,1)
    print('zarray is ',zarray)

if __name__ == "__main__":
   print('starting')
   plot()
