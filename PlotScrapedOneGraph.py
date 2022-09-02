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
max_orbits = 0

def setMaxOrbits():
    return none

def setScrapedJson():
    global scrapeDict
    with open('ScrapeDeltaTData_At'+str(max_orbits)+'.json', 'r') as f:
        scrapeDict = json.loads(f)

