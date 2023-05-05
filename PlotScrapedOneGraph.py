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

        temp = np.gradient(scrapeDict[key]["DTarray"])
        temp = np.convolve(temp,kernel,mode='valid')
            

        qnum = scrapeDict[key]["q"]

        q = "; q = " + str(scrapeDict[key]["q"])
        ecc = r"$\varepsilon$ = " + str(scrapeDict[key]["ecc_int"])

        if (qnum == 2e-4):
            if (ecc == 0.1):
                plt.plot(temp, label=ecc + q, ls=":",color="blue")
            if (ecc == 0.3):
                plt.plot(temp, label=ecc + q, ls=":", color="red")
            if (ecc == 0.5):
                plt.plot(temp, label=ecc + q, ls=":", color="orange")
            if (ecc == 0.7):
                plt.plot(temp, label=ecc + q, ls=":", color="green")
        else:
            if (ecc == 0.1):
                plt.plot(temp, label=ecc + q,color="blue")
            if (ecc == 0.3):
                plt.plot(temp, label=ecc + q, color="red")
            if (ecc == 0.5):
                plt.plot(temp, label=ecc + q, color="orange")
            if (ecc == 0.7):
                plt.plot(temp, label=ecc + q, color="green")


    plt.grid(True)
    plt.legend()
    plt.title('Log Temperature vs time')
    plt.xlabel(r'$t/T_0$')
    plt.ylabel(r'$\nabla$'+'(Temperature)')
    plt.ylim(0,0.05)
    plt.tight_layout()
    plt.savefig("ScrapedDeltaData_" + str(max_orbits) + "_plot" + ".png")
    plt.close()

    print("entering plot")

setScrapedJson()
plot()
