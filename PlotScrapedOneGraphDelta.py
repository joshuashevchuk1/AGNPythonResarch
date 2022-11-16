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

        N=50

        kernel = np.ones((N,))/N

        temp = np.gradient(scrapeDict[key]["DTarray"]-scrapeDict[key]["DSharray"])
        temp = np.convolve(temp,kernel,mode='valid')
            

        qnum = scrapeDict[key]["q"]

        q = "; q = " + str(scrapeDict[key]["q"])
        ecc = r"$\varepsilon$ = " + str(scrapeDict[key]["ecc_int"])

        if (qnum == 2e-4):
            plt.plot(temp, label=ecc + q,ls=":")
        else: 
            plt.plot(temp, label=ecc + q)

    plt.grid(True)
    plt.legend()
    plt.title('data')
    plt.xlabel(r'$t/T_0$')
    plt.ylabel(r'$\nabla'+'(Temperature)')
    plt.ylim(0,0.05)
    plt.tight_layout()
    plt.savefig("ScrapedDeltaTotalsData_" + str(max_orbits) + "_plot" + ".png")
    plt.close()

    print("entering plot")

setScrapedJson()
plot()
