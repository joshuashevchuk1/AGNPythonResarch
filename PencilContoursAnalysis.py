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

        self.Dir_beta_patches = None
        self.Dir_alpha_patches = None
        self.Dir_gamma_patches = None
        self.Dir_rsmooth_patches = None
        self.DirSigma_patches = None
        self.Dirinitial_pressure_patches = None
        self.Diraspect_ratio_patches = None
        self.Dirsound_speed_patches = None
        self.DirEcc_patches = None
        self.DirMass_patches = None
        self.DirName = None
        self.DirEcc = None
        self.DirMass = None
        self.beta = None
        self.alpha = None
        self.gamma = None
        self.rsmooth = None
        self.Sigma = None
        self.initial_pressure = None
        self.aspect_ratio = None
        self.sound_speed = None
        self.eccentricity = None
        self.ncolors = None
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
        print('returning to root')
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

    def initLoopEntry(self, data_frame,n):
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

    def initDirVars(self, data_frame , n):

        self.DirName = data_frame[n]['DirName']
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

        self.ncolors = 256

        self.eccentricity = data_frame[n]['eccentricity']
        self.sound_speed = data_frame[n]['cs']
        self.aspect_ratio = data_frame[n]['aspect_ratio']
        self.initial_pressure = data_frame[n]['initial_pressure']
        self.Sigma = data_frame[n]['Sigma']
        self.rsmooth = data_frame[n]['rsmooth']
        self.gamma = data_frame[n]['gamma']
        self.alpha = data_frame[n]['alpha']
        self.beta = data_frame[n]['beta']
        self.DirMass = data_frame[n]['par1']
        self.DirEcc = (round(eccentricity[0], 1))

        self.DirMass_patches = mpatches.Patch(
            color='white', label='q :' + str(self.DirMass))
        self.DirEcc_patches = mpatches.Patch(
            color='white', label=r'$\varepsilon$ :' + str(self.DirEcc))
        self.Dirsound_speed_patches = mpatches.Patch(
            color='white', label='sound speed :' + str(self.sound_speed))
        self.Diraspect_ratio_patches = mpatches.Patch(
            color='white', label='aspect_ratio :' + str(self.aspect_ratio))
        self.Dirinitial_pressure_patches = mpatches.Patch(
            color='white', label='inital pressure :' + str(self.initial_pressure))
        self.DirSigma_patches = mpatches.Patch(
            color='white', label=r'$\Sigma$ :' + str(self.Sigma))
        self.Dir_rsmooth_patches = mpatches.Patch(
            color='white', label='potential smoothing :' + str(self.rsmooth))
        self.Dir_gamma_patches = mpatches.Patch(
            color='white', label=r'$\gamma$ :' + str(self.gamma))
        self.Dir_alpha_patches = mpatches.Patch(
            color='white', label=r'$\alpha$ :' + str(self.alpha))
        self.Dir_beta_patches = mpatches.Patch(
            color='white', label=r'$\beta$ :' + str(self.beta))

    def pingContour(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                print("starting data loop")
                self.initLoopEntry(data_frame,n)
                self.initDirVars(data_frame,n)
                self.pingContourDensity()
                self.pingContourTemp()
                n = n + dn
        except:
            print('================')
            print('ping Contour Loop error')
            print('================')
            traceback.print_exc()
            self.returnToRoot()
            return False

    def pingContourDensity(self):
        try:
            fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
            fig.subplots_adjust(bottom=0.07, top=0.95)
            PL2 = ax1.contourf(x2d, y2d, rho_fv, ncolors)
            ax1.set_aspect('equal')
            plt.legend(handles=[DirMass_patches,
                                DirEcc_patches,
                                Dirsound_speed_patches,
                                Dirinitial_pressure_patches,
                                DirSigma_patches,
                                Dir_rsmooth_patches,
                                Dir_gamma_patches,
                                Dir_alpha_patches,
                                Dir_beta_patches], loc=2)
            cax = plt.axes([0.85, 0.1, 0.075, 0.8])
            cax.set_aspect(20)
            cax.set_ylabel('Density in code units', fontsize=10)
            plt.colorbar(PL2, cax=cax)
            plt.suptitle('t=' + str(data_frame[n]['ivar']))
            plt.subplots_adjust(bottom=0.05, top=0.95)
            plt.savefig('Standard_Contour_Density_' +
                        self.DirName + '_' + str(Standard_Orbit) + '.png')
            plt.close(fig)
            self.returnToRoot()
        except:
            print('================')
            print('No Contour density Data found to plot')
            print('================')
            traceback.print_exc()
            self.returnToRoot()

    def pingContourTemp(self):
        try:
            fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
            fig.subplots_adjust(bottom=0.07, top=0.95)
            PL2 = ax1.contourf(x2d, y2d, temp_fv, ncolors)
            ax1.set_aspect('equal')
            plt.legend(handles=[DirMass_patches,
                                DirEcc_patches,
                                Dirsound_speed_patches,
                                Dirinitial_pressure_patches,
                                DirSigma_patches,
                                Dir_rsmooth_patches,
                                Dir_gamma_patches,
                                Dir_alpha_patches,
                                Dir_beta_patches], loc=2)
            cax = plt.axes([0.85, 0.1, 0.075, 0.8])
            cax.set_aspect(20)
            cax.set_ylabel('Density in code units', fontsize=10)
            plt.colorbar(PL2, cax=cax)
            plt.suptitle('t=' + str(data_frame[n]['ivar']))
            plt.subplots_adjust(bottom=0.05, top=0.95)
            plt.savefig('Standard_Contour_Density_' +
                        self.DirName + '_' + str(Standard_Orbit) + '.png')
            plt.close(fig)
            self.returnToRoot()
        except:
            print('================')
            print('No Contour temp Data found to plot')
            print('================')
            traceback.print_exc()
            self.returnToRoot()
