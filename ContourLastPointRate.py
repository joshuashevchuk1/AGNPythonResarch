import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json
import scipy
import scipy.interpolate

x=[] # x=q
y=[] # y=e
z=[] # z=dT/dt
scrapDict={}

def plot():
    global scrapDict
    init()
    plotcontour()
    localDict={"x":x.tolist(),"y":y.tolist(),"z":z}
    scrapDict['test'] = localDict
    saveData()

def plotcontour():
    global x,y,z
    fig = plt.figure()

    x=[0.1,0.3,0.5,0.7]
    y=[1e-4,2e-4]
    x,y=np.meshgrid(x,y)

    z=[[0.127410097,0.4729206,0.16347564,0.00195316089],[0.05088384,0.066684518512,0.02821483,0.0268524]] # numbers got from json files

    print('x is ', x)
    print('y is', y)
    print('z is ', z)

    plt.contourf(x,y,z,256)
    plt.xlabel(r'\epsilon ')
    plt.ylabel('q')
    plt.title('q vs e')
    plt.show()

def saveData():
    global scrapDict
    data = scrapDict
    print('data is ',data)
    with open('data.json', 'w') as f:
        json.dump(data, f)

def init():
    print('initializing arrays')
    initXarray()
    initYarray()
    initZarray()

def initXarray():
    global x
    x=np.arange(0,5,1)
    #x=np.array(x)
    print('xarray is ',x)

def initYarray():
    global y
    y=np.arange(0,2,1)
    #y = np.array(y)
    print('yarray is ',y)

def initZarray():
    global z
    z=[np.arange(0,5,1),np.arange(0,2,1)]
    #z=np.array(z)
    print('zarray is ',z)

if __name__ == "__main__":
   print('starting')
   plot()
