import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json

x=[]
yarray=[]
zarray=[]

def plot():
    init()


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
