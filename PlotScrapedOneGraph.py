import pencil_old as pc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from pylab import *
import sys
import traceback
import json

scrapeDict={}
max_orbits = int(input('max Orbits? : '))


def plotCollectedData(paramDTarray):
    print("entering plotCollectedData")
    print("paramDTArray is ",paramDTarray)
    plt.plot(paramDTarray)
    print("leaving plotCollectedData")

def setScrapedJson():
    global scrapeDict
    with open('ScrapeDeltaTData_At'+str(max_orbits)+'.json') as f:
        scrapeDict = json.load(f)

def plot():
    global scrapeDict

    for key in scrapeDict:
        plt.plot(np.gradient(scrapeDict[key]["DTarray"]), label=
        r"$\varepsilon$ = " +
        str(scrapeDict[key]["ecc_int"]) +
        "; q = " +
        str(scrapeDict[key]["q"]))

    plt.grid(True)
    plt.legend()
    plt.title('data')
    plt.xlabel(r'$t/T_0$')
    plt.ylabel('Temperature')
    plt.tight_layout()
    plt.savefig("ScrapedDeltaData_" + str(max_orbits) + "_plot" + ".png")
    plt.close()

    print("entering plot")

setScrapedJson()
plot()