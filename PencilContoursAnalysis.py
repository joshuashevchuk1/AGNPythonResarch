import numpy
import os
import sys
import PencilContoursData as PD
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import traceback
import logging
import numpy as np
import math
import matplotlib.patches as mpatches
import scipy.interpolate
import CheckArrayTest as chk
import pandas as pd

class Pencil_Analysis(object):

    def __init__(self,
                 data_functions,
                 Orbit=None,
                 dir_run_list=None,
                 var_dir_list=None,
                 temp_dir_list=None,
                 Calc_Temp=False,
                 Calc_Density=False):

        self.Orbit = Orbit
        self.data_functions = data_functions
        self.dir_run_list = dir_run_list
        self.var_dir_list = var_dir_list
        self.temp_dir_list = temp_dir_list
        self.Calc_Temp = Calc_Temp
        self.Calc_Density = Calc_Density

        logging.basicConfig(filename='Pencil_Analysis.log',
                            level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                            datefmt=' %m/%d/%Y %I:%M:%S %p ')
        logging.info('constructor intializied')

    def Make_Vars(self):
                print('================')
                print('Making Vars')
                print('================')

                data_functions = self.data_functions
                Orbit = self.Orbit
                Calc_Temp = self.Calc_Temp
                Calc_Density = self.Calc_Density

                var_dir_list = []
                dir_run_list = next(os.walk('.'))[1]

                i = 0
                di = 1

                while i <= len(dir_run_list) - 1:
                    var_dir_list.append(
                        data_functions.grepDATA(
                            dir_run_list[i],
                            Orbit,
                            Calc_Temp=Calc_Temp,
                            Calc_Density=Calc_Density)
                    )
                    i = i + di
                return var_dir_list