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
                 MaxOrbits=None,
                 Orbit_standard=None,
                 step=None,
                 TempSigma=False,
                 dir_run_list=None,
                 var_dir_list=None,
                 temp_dir_list=None,
                 Calc_DTTm=None,
                 Calc_Temp=False,
                 Calc_Density=False,
                 Calc_Energy=False,
                 Calc_OEnergy=False,
                 Calc_Dynamics=False,
                 Calc_Rates_Energy=False,
                 Calc_ToomreQ=False):

        self.Orbit = Orbit
        self.Orbit_standard = Orbit_standard
        self.MaxOrbits = MaxOrbits
        self.data_functions = data_functions
        self.dir_run_list = dir_run_list
        self.var_dir_list = var_dir_list
        self.temp_dir_list = temp_dir_list
        self.step = step
        self.TempSigma = TempSigma
        self.Calc_DTTm = Calc_DTTm
        self.Calc_Temp = Calc_Temp
        self.Calc_Density = Calc_Density
        self.Calc_Energy = Calc_Energy
        self.Calc_OEnergy = Calc_OEnergy
        self.Calc_Dynamics = Calc_Dynamics
        self.Calc_Rates_Energy = Calc_Rates_Energy
        self.Calc_ToomreQ = Calc_ToomreQ

        logging.basicConfig(filename='Pencil_Analysis.log',
                            level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                            datefmt=' %m/%d/%Y %I:%M:%S %p ')
        logging.info('constructor intializied')

    def Make_Vars(self):
                print('================')
                print('Making Vars')
                print('================')

                data_functions = self.data_functions
                dir_run_list = self.dir_run_list
                var_dir_list = self.var_dir_list
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
                    data_functions.pingGit(dir_run_list[i])
                    data_functions.pingTemp(dir_run_list[i])
                    i = i + di
                return var_dir_list