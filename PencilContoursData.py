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
                 Calc_Temp=False,
                 Calc_Density=False,
                ):

        self.Directory_Path = Directory_Path
        self.Orbit = Orbit
        self.Calc_Temp = Calc_Temp
        self.Calc_Density = Calc_Density

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
