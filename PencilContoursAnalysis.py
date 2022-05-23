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

        self.Standard_Orbit = None
        self.avgshock_fv = None
        self.avgtemp_fv = None
        self.avgrho_fv = None
        self.Init_Temp = None
        self.shock_fv = None
        self.temp_fv = None
        self.rho_fv = None
        self.yrq2 = None
        self.xrq2 = None
        self.y2d = None
        self.x2d = None
        self.rad_grid = None
        
        
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

    def returnToRoot(self):
        os.chdir('..')

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

    def initLoopEntry(self,data_frame):
        try:
            os.system('mkdir Pencil_Analysis')
            os.chdir('Pencil_Analysis')
        except:
            print('================')
            print('directory already found')
            print('================')
            os.chdir('Pencil_Analysis')
        try:
            print('================')
            print('Making Contour plots for ' +
                  str(data_frame[n]['DirName']))
            print('================')
        except:
            print('================')
            print('Dir error')
            print('================')
            print(n)
            print('================')
            logging.basicConfig(filename='Pencil_Analysis.log',
                                level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                                datefmt=' %m/%d/%Y %I:%M:%S %p ')
            logging.info('========================')
            logging.info(n)
            logging.info('========================')

    def initDirVars(self,data_frame):
        self.rad_grid = data_frame[n]['rad_grid']
        self.x2d = data_frame[n]['x2d']
        self.y2d = data_frame[n]['y2d']
        self.xrq2 = data_frame[n]['xrq2']
        self.yrq2 = data_frame[n]['yrq2']
        self.rho_fv = data_frame[n]['rho_fv']
        self.temp_fv = data_frame[n]['temp_fv']
        self.shock_fv = data_frame[n]['shock_fv']
        self.Init_Temp = data_frame[n]['Init_Temp']
        self.avgrho_fv = data_frame[n]['avgrho_fv']
        self.avgtemp_fv = data_frame[n]['avgtemp_fv']
        self.avgshock_fv = data_frame[n]['avgshock_fv']
        self.Standard_Orbit = data_frame[n]['Standard_Orbit']

    def pingContour(self,data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                print("starting data loop")
                self.initLoopEntry(data_frame)
                self.initDirVars(data_frame)
                self.pingContourDensity(data_frame)
                self.pingContourTemp(data_frame)
                n = n + dn
        except:
            print('================')
            print('ping Contour Loop error')
            print('================')
            traceback.print_exc()
            self.returnToRoot()
            return False

    def pingContourDensity(self, data_frame):
        try:
            self.returnToRoot()
        except:
            print('================')
            print('No Contour density Data found to plot')
            print('================')
            traceback.print_exc()
            self.returnToRoot()
            
    def pingContourTemp(self, data_frame):
        try:
            self.returnToRoot()
        except:
            print('================')
            print('No Contour temp Data found to plot')
            print('================')
            traceback.print_exc()
            self.returnToRoot()