import os
import numpy as np
import pencil_old as pc
import math
import sys
import traceback
import numpy.polynomial.hermite as herm
import scipy
from scipy.optimize import curve_fit
import logging
import re


class Pencil_Data(object):

    def __init__(self):
        print('================')
        print('initalize vars')
        print('================')

    def grepDATA(self,
                 Directory_Path,
                 Orbit=None,
                 MaxOrbits=10,
                 Calc_Temp=False,
                 Calc_Density=False,
                 step=None,
                 Orbit_standard=None,
                 Calc_Energy=False,
                 Calc_OEnergy=False,
                 Calc_Dynamics=False,
                 Calc_Rates_Energy=False,
                 Calc_ToomreQ=False):

        try:
            vars_dict = {}
            os.chdir('..')
            return vars_dict
        except:
            print('=======================')
            print('No data')
            print(Directory_Path)
            print('=======================')
            traceback.print_exc()
            print('=======================')
            os.chdir('..')
            return False
