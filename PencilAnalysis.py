# ================================================================
#
#   each sub directory in the cwd
#   run selected analysis scripts
#   send all data through github
#
# ================================================================

import numpy
import os
import sys
import PencilData as PD
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


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        return array[idx - 1]
    else:
        return array[idx]


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
        step = self.step
        temp_dir_list = self.temp_dir_list
        TempSigma = self.TempSigma
        Orbit_standard = self.Orbit_standard
        Calc_Temp = self.Calc_Temp
        Calc_Density = self.Calc_Density
        Calc_Energy = self.Calc_Energy
        Calc_Dynamics = self.Calc_Dynamics
        Calc_Rates_Energy = self.Calc_Rates_Energy
        Calc_ToomreQ = self.Calc_ToomreQ

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
                    Calc_Density=Calc_Density,
                    step=step,
                    Orbit_standard=Orbit_standard,
                    Calc_Energy=Calc_Energy,
                    Calc_Dynamics=Calc_Dynamics,
                    Calc_Rates_Energy=Calc_Rates_Energy,
                    Calc_ToomreQ=Calc_ToomreQ)
            )
            data_functions.pingGit(dir_run_list[i])
            data_functions.pingTemp(dir_run_list[i])
            i = i + di

        return var_dir_list

    def pingToomreQ(self, data_frame):

        Calc_ToomreQ = self.Calc_ToomreQ

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Toomre plot for ' +
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

                    Toomre = data_frame[n]['Toomre']
                    x2d = data_frame[n]['x2d']
                    y2d = data_frame[n]['y2d']
                    ncolors = 256

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    PL2 = ax1.contourf(x2d, y2d, Toomre, ncolors)
                    ax1.set_aspect('equal')
                    # plt.contourf(x2d, y2d, Toomre, ncolors)

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('Density in code units', fontsize=10)

                    plt.colorbar(PL2, cax=cax)
                    plt.ylabel('Q eff')
                    plt.xlabel(r'$\tau$')
                    plt.title('mean Toomre Q eff per orbit')
                    plt.savefig('Standard_ToomreQ_' +
                                str(data_frame[n]['DirName']) + '.png')

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Toomre Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping ToomreQ Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingEnergy(self, data_frame):

        Calc_Rates_Energy = self.Calc_Rates_Energy

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Energy plot for ' +
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

                    KE_Sum = data_frame[n]['KE_Sum']
                    UE_Sum = data_frame[n]['UE_Sum']
                    UINT_Sum = data_frame[n]['UINT_Sum']
                    Total = data_frame[n]['Total_Disk_Energy']
                    Total_Avg = data_frame[n]['Total_Disk_Energy_Avg']
                    OE_Sum = data_frame[n]['OE_Sum']

                    KE_fit = data_frame[n]['KE_fit']
                    UE_fit = data_frame[n]['UE_fit']
                    UINT_fit = data_frame[n]['UINT_fit']

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run
                    # ping energy has unique handles do not copy for other plots

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    time = data_frame[n]['time']
                    cut_off = data_frame[n]['indexTimeCutOff']

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))
                    KE_Sum_Label = mpatches.Patch(
                        color='green', label='Kinetic')
                    UE_Sum_Label = mpatches.Patch(
                        color='red', label='Potential')
                    UINT_Sum_Label = mpatches.Patch(
                        color='purple', label='Internal')
                    OE_Sum_Label = mpatches.Patch(
                        color='orange', label='Orbital')
                    Dir_cutoff_patches = mpatches.Patch(
                        color='white', label=r'$t_{c}$ :' + str(time[cut_off - 1]))

                    # ---------------------------------------------------------------------
                    # =========================
                    #
                    # plot all relevant energy data
                    # if calc energy is called, plot those as well
                    #
                    # =========================

                    # =========================
                    # Energy magnitude
                    # =========================

                    plt.figure(figsize=(10, 10))
                    plt.title('Energy Magnitude Comparison')
                    plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
                    plt.plot(UE_Sum, color='red', label='Potential')
                    plt.plot(UINT_Sum, color='purple',
                             label='Internal', ls='-.')
                    plt.plot(np.abs(Total[:]), color='blue',
                             label='Total disk energy')
                    plt.plot(
                        np.abs(OE_Sum[:]), color='orange', label='Orbital', ls=':')
                    plt.yscale('log')
                    plt.xlabel('time (Orbits)', fontweight='bold')
                    plt.ylabel('Joule (Code units)', fontweight='bold')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        KE_Sum_Label,
                                        UE_Sum_Label,
                                        UINT_Sum_Label,
                                        OE_Sum_Label,
                                        Dir_cutoff_patches], loc=2)
                    plt.savefig('Standard_ErgMC_' +
                                str(data_frame[n]['DirName']) + '.png')
                    plt.close()

                    # =========================
                    # Total Energy
                    # =========================

                    plt.figure(figsize=(10, 10))
                    plt.title('Energy Magnitude Comparison')
                    plt.plot(np.abs(Total_Avg[:]), color='blue',
                             label='Total disk energy', ls='--')
                    plt.xlabel('time (Orbits)', fontweight='bold')
                    plt.ylabel('Joule (Code units)', fontweight='bold')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        Dir_cutoff_patches], loc=2)
                    plt.savefig('Standard_ErgTot_' +
                                str(data_frame[n]['DirName']) + '.png')
                    plt.close()

                    # =========================
                    # Potential Energy of the disk
                    # ========================

                    plt.figure(figsize=(10, 10))
                    plt.title('Disk Potential Energy')
                    plt.plot(UE_Sum, color='red',
                             label='Disk Potential Energy', ls='--')
                    if Calc_Rates_Energy == True:
                        # plot rates and check error
                        plt.plot(UE_fit, color='black', label='UE fit', ls=':')
                    plt.xlabel('time (Orbits)', fontweight='bold')
                    plt.ylabel('Joule (Code units)', fontweight='bold')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        UE_Sum_Label,
                                        Dir_cutoff_patches], loc=2)
                    plt.savefig('Standard_ErgUE_' +
                                str(data_frame[n]['DirName']) + '.png')
                    plt.close()

                    # ========================
                    # Kinetic Energy of the disk
                    # ========================

                    plt.figure(figsize=(10, 10))
                    plt.title('Disk Kinetic Energy')
                    plt.plot(KE_Sum, color='green',
                             label='Disk Kinetic Energy', ls='--')
                    if Calc_Rates_Energy == True:
                        # plot rates and check error
                        plt.plot(KE_fit, color='black', label='KE fit', ls=':')
                    plt.xlabel('time (Orbits)', fontweight='bold')
                    plt.ylabel('Joule (Code units)', fontweight='bold')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        KE_Sum_Label,
                                        Dir_cutoff_patches], loc=2)
                    plt.savefig('Standard_ErgKE_' +
                                str(data_frame[n]['DirName']) + '.png')
                    plt.close()

                    # ========================
                    # Internal Energy of the disk
                    # ========================

                    plt.figure(figsize=(10, 10))
                    plt.title('Disk Internal Energy')
                    plt.plot(UINT_Sum, color='purple',
                             label='Disk Internal Energy', ls='--')
                    if Calc_Rates_Energy == True:
                        # plot rates and check error
                        plt.plot(UINT_fit, color='black',
                                 label='UINT fit', ls=':')
                    plt.xlabel('time (Orbits)', fontweight='bold')
                    plt.ylabel('Joule (Code units)', fontweight='bold')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        UINT_Sum_Label,
                                        Dir_cutoff_patches], loc=2)
                    plt.savefig('Standard_ErgUINT_' +
                                str(data_frame[n]['DirName']) + '.png')
                    plt.close()

                    # =========================
                    # ---------------------------------------------------------------------

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Energy Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Energy Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def Make_TempVars(self, data_frame):

        Calc_Temp = self.Calc_Temp

        TempSigma = self.TempSigma
        temp_dir_list = self.temp_dir_list

        if TempSigma == False:
            temp_dir_list = None
        else:
            temp_dir_list = []
            ecc_dir_list = []
            name_dir_list = []

            i = 0
            di = 1
            while i <= len(data_frame) - 1:
                try:
                    if data_frame[i]['tempdata'] == None:
                        i = i + di
                    else:
                        temp_dir_list.append(data_frame[i]['tempdata'])
                        try:
                            eccentricity = data_frame[i]['eccentricity']
                            ecc_dir_list.append(round(eccentricity[0], 3))
                            name_dir_list.append(data_frame[i]['DirName'])
                        except:
                            print('================')
                            traceback.print_exc()
                            print('================')

                        print('================')
                        print('Adding Temp Data')
                        print('================')
                        i = i + di
                except:
                    traceback.print_exc()
                    print('================')
                    print('No Temp Data found to plot')
                    print('================')
                i = i + di

        return temp_dir_list, ecc_dir_list, name_dir_list

    def Get_TempDistribution(self,
                             temp_frame,
                             ecc_frame=None,
                             name_frame=None):
        try:
            print('================')
            print('getting temp distribution')
            print('================')

            # attempt to get temp distribution
            i = 0
            di = 1

            f = 0
            df = 1

            n = 0
            dn = 1

            temp_distribution = []
            temp_ecc = []
            temp_name = []

            while f <= len(temp_frame) - 1:
                try:
                    if math.isnan(np.std(temp_frame[f][:])) == False:
                        temp_distribution.append(np.std(temp_frame[f][:]))
                        if ecc_frame == None or name_frame == None:
                            f = f + df
                        else:
                            temp_ecc.append(ecc_frame[f])
                            temp_name.append(name_frame[f])
                            f = f + df
                    else:
                        f = f + df
                except:
                    print('f failed at ' + str(f))
                    f = f + df
        except:
            print('================')
            print('================')
            print('Temp_distribution error')
            print('================')
            print('================')
            traceback.print_exc()
            print('================')
            print('================')

        return temp_distribution, temp_ecc, temp_name

    def pingGTM_Sigma(self, data_frame):

        GTM_Sigma_AllOrbits = []
        DirMass_GTM = []
        DirEcc_GTM = []

        i = 0
        di = 1

        while i <= len(data_frame) - 1:
            try:
                try:
                    GTM_Sigma_AllOrbits.append(data_frame[i]['GTM_Sigma'])
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                except:
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                    print('================')
                    print('no Data for GTM_Sigma ')
                    print('================')
                    traceback.print_exc()
            except:
                print('================')
                print('error in dir' + str(i))
                print('================')
                traceback.print_exc()
            i = i + di

        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig = plt.figure(figsize=(10, 10))
            # ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(111)

            DirEcc_GTM = np.ravel(DirEcc_GTM)
            DirMass_GTM = np.ravel(DirMass_GTM)
            GTM_Sigma_AllOrbits = np.ravel(GTM_Sigma_AllOrbits)

            ax.scatter(DirEcc_GTM, DirMass_GTM,
                       c=GTM_Sigma_AllOrbits, marker='o')

            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
            ax.set_ylim(0, 1e-4)

            plt.title(
                r'global $\sigma_{T}$ as a function of mass ratio and eccentricity')
            ax.set_ylabel('q', fontsize=14)
            ax.set_xlabel(r'$\varepsilon$', fontsize=14)
            # ax.set_zlabel(r'global $\sigma_{T}$',fontsize=14)
            # plt.show()

            plt.savefig('Global_SigmaTemp(e,q)_Parameter_Space.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No Temp Sigma found to plot')
            print('================')
            os.chdir('..')

        # =========================
        #
        #   get information for slices of contour plots.
        #   we want a accurate represenation of global slices
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(DirEcc_GTM) - 1:
                try:
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
                        print('Making Sigma Temp Slices plot for ' +
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

                    # =========================
                    #
                    # 3D plot slices
                    #
                    # =========================
                    #
                    # array needs to be related to DirEcc_GTM
                    #
                    # for a fixed ecc (0.1 could be the first)
                    # get all data_frame[n] compenents that have this fixed ecc
                    # get there respective q (mass ratio), and GTMSigma_Temp(standard deviation of mean disk temperature)
                    # make new arrays containing q, and sigma T for the constant ecc.
                    # iterate until the maximum eccentricity for the 3D plot has been reached
                    #
                    # =========================

                    Mass_Slice = []
                    GTM_Sigma_Slice = []

                    i = 0
                    di = 1

                    while i <= len(DirEcc_GTM) - 1:
                        if DirEcc_GTM[i] == DirEcc_GTM[n]:
                            Mass_Slice.append(DirMass_GTM[i])
                            GTM_Sigma_Slice.append(GTM_Sigma_AllOrbits[i])
                            i = i + di
                        else:
                            i = i + di

                    plt.scatter(Mass_Slice, GTM_Sigma_Slice, color='purple')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.xlabel('q')
                    plt.ylabel(r'global $\sigma_{T}$')
                    plt.ticklabel_format(
                        style='sci', axis='x', scilimits=(0, 0))
                    plt.ticklabel_format(
                        style='sci', axis='y', scilimits=(0, 0))
                    plt.xlim(0, max(Mass_Slice))
                    plt.ylim(0, max(GTM_Sigma_Slice))

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Global_SigmaT_Slices_' +
                                str(DirEcc_GTM[n]) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)
                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Sigma Temp Slices found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Sigma Temp slices Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingGTM_Sigma_O(self, data_frame):

        GTM_Sigma_AllOrbits = []
        DirMass_GTM = []
        DirEcc_GTM = []

        i = 0
        di = 1

        while i <= len(data_frame) - 1:
            try:
                try:
                    GTM_Sigma_AllOrbits.append(data_frame[i]['GTM_Sigma_O'])
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                except:
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                    print('================')
                    print('no Data for GTM_Sigma in')
                    print('================')
                    traceback.print_exc()
            except:
                print('================')
                print('error in dir' + str(i))
                print('================')
                traceback.print_exc()

            i = i + di

        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig = plt.figure(figsize=(10, 10))
            # ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(111)

            # =========================
            # parameters for contourplot
            # =========================

            DirEcc_GTMi = np.linspace(min(DirEcc_GTM), max(DirEcc_GTM), 500)
            DirMass_GTMi = np.linspace(min(DirMass_GTM), max(DirMass_GTM), 500)
            DirEcc_GTMi, DirMass_GTMi = np.meshgrid(DirEcc_GTMi, DirMass_GTMi)
            GTM_Sigma_AllOrbitsi = scipy.interpolate.griddata(
                (DirEcc_GTM, DirMass_GTM),
                GTM_Sigma_AllOrbits,
                (DirEcc_GTMi,
                 DirMass_GTMi),
                method='linear')

            im = ax.contourf(DirEcc_GTMi, DirMass_GTMi, GTM_Sigma_AllOrbitsi)

            cax = plt.axes([0.85, 0.1, 0.075, 0.8])
            cax.set_aspect(20)
            fig.colorbar(im, cax=cax)

            # =========================

            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
            ax.set_ylim(0, 1e-4)

            plt.suptitle(
                r'global $\sigma_{T}$ as a function of mass ratio and eccentricity')
            ax.set_ylabel('q', fontsize=14)
            ax.set_xlabel(r'$\varepsilon$', fontsize=14)
            # ax.set_zlabel(r'global $\sigma_{T}$',fontsize=14)
            # plt.show()

            plt.savefig('Global_SigmaTemp_O_(e,q)_Parameter_Space.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No Temp Sigma found to plot')
            print('================')
            os.chdir('..')

        # =========================
        #
        #   get information for slices of contour plots.
        #   we want a accurate represenation of global slices
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(DirEcc_GTM) - 1:
                try:
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
                        print('Making Sigma Temp Slices plot for ' +
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

                    # =========================
                    #
                    # 3D plot slices
                    #
                    # =========================
                    #
                    # array needs to be related to DirEcc_GTM
                    #
                    # for a fixed ecc (0.1 could be the first)
                    # get all data_frame[n] compenents that have this fixed ecc
                    # get there respective q (mass ratio), and GTMSigma_Temp(standard deviation of mean disk temperature)
                    # make new arrays containing q, and sigma T for the constant ecc.
                    # iterate until the maximum eccentricity for the 3D plot has been reached
                    #
                    # =========================

                    Mass_Slice = []
                    GTM_Sigma_Slice = []

                    i = 0
                    di = 1

                    while i <= len(DirEcc_GTM) - 1:
                        if DirEcc_GTM[i] == DirEcc_GTM[n]:
                            Mass_Slice.append(DirMass_GTM[i])
                            GTM_Sigma_Slice.append(GTM_Sigma_AllOrbits[i])
                            i = i + di
                        else:
                            i = i + di

                    plt.scatter(Mass_Slice, GTM_Sigma_Slice, color='purple')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.xlabel('q')
                    plt.ylabel(r'global $\sigma_{T}$')
                    plt.ticklabel_format(
                        style='sci', axis='x', scilimits=(0, 0))
                    plt.ticklabel_format(
                        style='sci', axis='y', scilimits=(0, 0))
                    plt.xlim(0, max(Mass_Slice))
                    plt.ylim(0, max(GTM_Sigma_Slice))

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Global_SigmaTi_O_Slices_' +
                                str(DirEcc_GTM[n]) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)
                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Sigma Temp Slices found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Sigma Temp slices Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingGTM_O_DT(self, data_frame):

        GTM_Sigma_AllOrbits = []
        DirMass_GTM = []
        DirEcc_GTM = []

        i = 0
        di = 1

        while i <= len(data_frame) - 1:
            try:
                try:
                    GTM_Sigma_AllOrbits.append(data_frame[i]['DGTemp_Mean'])
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                except:
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                    print('================')
                    print('no Data for GTM_Sigma ')
                    print('================')
                    traceback.print_exc()
            except:
                print('================')
                print('error in dir' + str(i))
                print('================')
                traceback.print_exc()

            i = i + di

        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig = plt.figure(figsize=(10, 10))
            # ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(111)

            # =========================
            # parameters for contourplot
            # =========================

            DirEcc_GTMi = np.linspace(min(DirEcc_GTM), max(DirEcc_GTM), 500)
            DirMass_GTMi = np.linspace(min(DirMass_GTM), max(DirMass_GTM), 500)
            DirEcc_GTMi, DirMass_GTMi = np.meshgrid(DirEcc_GTMi, DirMass_GTMi)
            GTM_Sigma_AllOrbitsi = scipy.interpolate.griddata(
                (DirEcc_GTM, DirMass_GTM),
                GTM_Sigma_AllOrbits,
                (DirEcc_GTMi,
                 DirMass_GTMi),
                method='linear')

            im = ax.contourf(DirEcc_GTMi, DirMass_GTMi, GTM_Sigma_AllOrbitsi)

            cax = plt.axes([0.85, 0.1, 0.075, 0.8])
            cax.set_aspect(20)
            fig.colorbar(im, cax=cax)

            # =========================

            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
            ax.set_ylim(0, 1e-4)

            plt.suptitle(
                r'$\frac{\partial T}{\partial t}$ as a function of mass ratio and eccentricity')
            ax.set_ylabel('q', fontsize=14)
            ax.set_xlabel(r'$\varepsilon$', fontsize=14)
            # ax.set_zlabel(r'global $\sigma_{T}$',fontsize=14)
            # plt.show()

            plt.savefig('Global_dT_O_(e,q)_Parameter_Space.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No DTemp found to plot')
            print('================')
            os.chdir('..')

        # =========================
        #
        #   get information for slices of contour plots.
        #   we want a accurate represenation of global slices
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(DirEcc_GTM) - 1:
                try:
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
                        print('Making DTemp Slices plot for ' +
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

                    # =========================
                    #
                    # 3D plot slices
                    #
                    # =========================
                    #
                    # array needs to be related to DirEcc_GTM
                    #
                    # for a fixed ecc (0.1 could be the first)
                    # get all data_frame[n] compenents that have this fixed ecc
                    # get there respective q (mass ratio), and GTMSigma_Temp(standard deviation of mean disk temperature)
                    # make new arrays containing q, and sigma T for the constant ecc.
                    # iterate until the maximum eccentricity for the 3D plot has been reached
                    #
                    # =========================

                    Mass_Slice = []
                    GTM_Sigma_Slice = []

                    i = 0
                    di = 1

                    while i <= len(DirEcc_GTM) - 1:
                        if DirEcc_GTM[i] == DirEcc_GTM[n]:
                            Mass_Slice.append(DirMass_GTM[i])
                            GTM_Sigma_Slice.append(GTM_Sigma_AllOrbits[i])
                            i = i + di
                        else:
                            i = i + di

                    plt.scatter(Mass_Slice, GTM_Sigma_Slice, color='purple')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.xlabel('q')
                    plt.ylabel(r'dT$')
                    plt.ticklabel_format(
                        style='sci', axis='x', scilimits=(0, 0))
                    plt.ticklabel_format(
                        style='sci', axis='y', scilimits=(0, 0))
                    plt.xlim(0, max(Mass_Slice))
                    plt.ylim(0, max(GTM_Sigma_Slice))

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Global_dT_O_Slices_' +
                                str(DirEcc_GTM[n]) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)
                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No d Temp Slices found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping d Temp slices Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingGTM_Sigma_3D(self, data_frame):

        GTM_Sigma_AllOrbits = []
        DirMass_GTM = []
        DirEcc_GTM = []

        i = 0
        di = 1

        while i <= len(data_frame) - 1:
            try:
                try:
                    GTM_Sigma_AllOrbits.append(data_frame[i]['GTM_Sigma'])
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                except:
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                    print('================')
                    print('no Data for GTM_Sigma in')
                    print('================')
                    traceback.print_exc()
            except:
                print('================')
                print('error in dir' + str(i))
                print('================')
                traceback.print_exc()
            i = i + di

        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            # ax  = fig.add_subplot(111)
            ax.scatter(DirEcc_GTM, DirMass_GTM,
                       GTM_Sigma_AllOrbits, marker='o')

            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
            ax.set_ylim(0, 1e-4)

            plt.title(
                r'global $\sigma_{T}$ as a function of mass ratio and eccentricity')
            ax.set_ylabel('q', fontsize=14)
            ax.set_xlabel(r'$\varepsilon$', fontsize=14)
            ax.set_zlabel(r'global $\sigma_{T}$', fontsize=14)
            # plt.show()

            plt.savefig('Global_SigmaTemp(e,q)_Parameter_Space_3D.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No Temp Sigma found to plot')
            print('================')
            os.chdir('..')

    def pingGTM_Sigma_3D_O(self, data_frame):

        GTM_Sigma_AllOrbits = []
        DirMass_GTM = []
        DirEcc_GTM = []

        i = 0
        di = 1

        while i <= len(data_frame) - 1:
            try:
                try:
                    GTM_Sigma_AllOrbits.append(data_frame[i]['GTM_Sigma_O'])
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                except:
                    DirMass_GTM.append(data_frame[i]['par1'])
                    eccentricity = data_frame[i]['eccentricity']
                    DirEcc_GTM.append(round(eccentricity[0], 1))
                    print('================')
                    print('no Data for GTM_Sigma in')
                    print('================')
                    traceback.print_exc()
            except:
                print('================')
                print('error in dir' + str(i))
                print('================')
                traceback.print_exc()
            i = i + di

        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            # ax  = fig.add_subplot(111)
            ax.scatter(DirEcc_GTM, DirMass_GTM,
                       GTM_Sigma_AllOrbits, marker='o')

            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
            ax.set_ylim(0, 1e-4)

            plt.title(
                r'global $\sigma_{T}$ as a function of mass ratio and eccentricity')
            ax.set_ylabel('q', fontsize=14)
            ax.set_xlabel(r'$\varepsilon$', fontsize=14)
            ax.set_zlabel(r'global $\sigma_{T}$', fontsize=14)
            # plt.show()

            plt.savefig('Global_SigmaTemp_O_(e,q)_Parameter_Space_3D.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No Temp Sigma found to plot')
            print('================')
            os.chdir('..')

    def pingTempDplot(self,
                      dist,
                      ecc,
                      name):
        try:
            try:
                os.system('mkdir Pencil_Analysis')
                os.chdir('Pencil_Analysis')
            except:
                print('================')
                print('directory already found')
                print('================')
                os.chdir('Pencil_Analysis')

            fig, (ax1) = plt.subplots(1, 1, figsize=(20, 10))
            fig.subplots_adjust(bottom=0.07, top=0.95, left=0.07)
            ax1.scatter(name, dist, color='blue')
            ax1.plot(name, dist, ls='--', color='orange', alpha=0.3)
            ax1.set_aspect('auto')
            plt.tight_layout()
            plt.grid(True)
            plt.xlabel('Data Directory')
            plt.ylabel(r'$\sigma$ (global)')
            plt.title('Global Temp ' + r'$\sigma$')
            plt.subplots_adjust(bottom=0.05, top=0.95, left=0.10)
            plt.savefig('Global_Temp_Distribution.png')
            # too many figs may be open.
            # close plots once done
            plt.close(fig)
            os.system('git add *.png')
            os.chdir('..')

        except:
            print('================')
            print('No Temp Distribution found to plot')
            print('================')
            os.chdir('..')

    def pingTTm(self, data_frame):

        Calc_Temp = self.Calc_Temp

        # =========================
        #
        # Get the longitude of Perihelion to check the precession of the orbiter
        # should only be included for global disk runs or anything that may cause point masses to precess
        #
        # =========================

        if self.Calc_Temp == True:
            #
            #
            #
            try:
                n = 0
                dn = 1
                while n <= len(data_frame) - 1:
                    try:
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
                            print('Making TTm vs time plot for ' +
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

                        time = data_frame[n]['t']
                        cut_off = data_frame[n]['indexTimeCutOff']
                        time = time[:cut_off] / (2 * np.pi)

                        print('time is ', time)
                        print('time[:cut_off] is', time[:cut_off])

                        GlobalTemp_Mean = data_frame[n]['GlobalTemp_Mean']
                        TTm = data_frame[n]['TTm']
                        TTm_rate = data_frame[n]['TTm_rate']

                        # =========================
                        # legend handles
                        # =========================

                        # use this information when adding important parameters to the run

                        eccentricity = data_frame[n]['eccentricity']
                        sound_speed = data_frame[n]['cs']
                        aspect_ratio = data_frame[n]['aspect_ratio']
                        initial_pressure = data_frame[n]['initial_pressure']
                        Sigma = data_frame[n]['Sigma']
                        EntropyIndex = data_frame[n]['EntropyIndex']
                        SpecificHeat = data_frame[n]['SpecificHeat']
                        rsmooth = data_frame[n]['rsmooth']
                        gamma = data_frame[n]['gamma']
                        Gamma0 = data_frame[n]['Gamma0']
                        alpha = data_frame[n]['alpha']
                        beta = data_frame[n]['beta']

                        deltaTemp = GlobalTemp_Mean[:] - GlobalTemp_Mean[0]

                        # DEBUG
                        print(GlobalTemp_Mean)
                        print(TTm)

                        plt.plot(
                            time, deltaTemp[:len(time)], color='purple', ls=':')
                        plt.title('Mean disk temperature vs time')
                        plt.xlabel(r'$t/T_0$')
                        plt.ylabel('Temperature')
                        plt.tight_layout()
                        plt.grid(True)
                        plt.xlim([0, time.max()])

                        DirEcc = (round(eccentricity[0], 1))

                        DirMass = data_frame[n]['par1']

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))
                        Dir_TTm_patches = mpatches.Patch(
                            color='white', label=r'$t_{m}$ :' + str(np.round(deltaTemp[:len(time)].max())))
                        Dir_cutoff_patches = mpatches.Patch(
                            color='white', label=r'$t_{c}$ :' + str(time[cut_off - 1]))

                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches,
                                            Dir_TTm_patches,
                                            Dir_cutoff_patches], loc=2)

                        # =========================

                        plt.savefig('Standard_TTm_' +
                                    str(data_frame[n]['DirName']) + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close()
                        os.system('git add *.png')
                        os.chdir('..')
                    except:
                        print('================')
                        print('No TTm  found to plot')
                        print('================')
                        traceback.print_exc()
                        os.chdir('..')
                    n = n + dn
            except:
                print('================')
                print('ping TTm Loop error')
                print('================')
                traceback.print_exc()
                os.chdir('..')
                return False
        else:
            print('================')
            print('No Temperature to plot moving on (TTm)')
            print('================')

    def ping_LongPerihelion(self, data_frame):

        # =========================
        #
        # Get the longitude of Perihelion to check the precession of the orbiter
        # should only be included for global disk runs or anything that may cause point masses to precess
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Perihelion plot for ' +
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

                    plt.plot(data_frame[n]['t'], data_frame[n]
                    ['long_perihelion'], color='red')
                    plt.title('Longitude of perihelion vs time')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\omega(rad)$')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.xlim([0, data_frame[n]['tmax']])

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Standard_Precession_Check_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close()
                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Perhilion Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Torque Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def ping_LongPerihelion_10_Orbits(self, data_frame):

        # =========================
        #
        # Get the longitude of Perihelion to check the precession of the orbiter
        # should only be included for global disk runs or anything that may cause point masses to precess
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Perihelion plot for ' +
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

                    plt.plot(data_frame[n]['t'], data_frame[n]
                    ['long_perihelion'], color='red')
                    plt.title('Longitude of perihelion vs time')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\omega(rad)$')
                    plt.tight_layout()
                    plt.grid(True)
                    plt.xlim([0, 10])

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Standard_Precession_Check_10_Orbits' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close()
                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Perhilion Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Torque Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingTorque(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Torque plot for ' +
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

                    time = data_frame[n]['time']

                    # init for cut_off
                    cut_off = data_frame[n]['indexTimeCutOff']
                    time = time[:cut_off]
                    torqint = data_frame[n]['torqint']
                    torqext = data_frame[n]['torqext']
                    torqtotal = data_frame[n]['torqtotal']
                    #

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.plot(time, torqint[:len(time)], '--', label='Inner')
                    ax1.plot(time, torqext[:len(time)], '--', label='Outer')
                    ax1.plot(time, torqtotal[:len(time)], label='Total')
                    ax1.set_aspect('auto')
                    plt.xlim([0, time.max()])
                    plt.tight_layout()
                    plt.grid(True)

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex  = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))
                    Dir_cutoff_patches = mpatches.Patch(
                        color='white', label=r'$t_{c}$ :' + str(time[cut_off - 1]))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches,
                                        Dir_cutoff_patches], loc=2)

                    # =========================

                    plt.title(r'Torques')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\Gamma$')
                    plt.xlim([0, time[cut_off - 1]])
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Torque_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)
                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Torque Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Torque Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingContour_Dynamic(self, data_frame):

        # get contour plots in slices of orbits

        Calc_Temp = self.Calc_Temp

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Dynamic Contour plots for ' +
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

                    x2d = data_frame[n]['x2d']
                    y2d = data_frame[n]['y2d']
                    xrq2 = data_frame[n]['xrq2']
                    yrq2 = data_frame[n]['yrq2']
                    Dynamic_Density = data_frame[n]['Dynamic_Density']
                    Dynamic_Temperature = data_frame[n]['Dynamic_Temperature']
                    Dynamic_Shock = data_frame[n]['Dynamic_Shock']
                    Standard_Orbit = data_frame[n]['Standard_Orbit']
                    Init_Temp = data_frame[n]['Init_Temp']

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    Dynamic_Loop = 0
                    dLoop = 1
                    ncolors = 256

                    print('================')
                    print('Making Dynamic_Density')
                    print('================')

                    while Dynamic_Loop <= len(Dynamic_Density) - 1:
                        # =========================
                        # =========================
                        # Density Contour
                        # =========================
                        # =========================

                        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        PL2 = ax1.contourf(
                            x2d, y2d, Dynamic_Density[Dynamic_Loop][:], ncolors)
                        ax1.set_aspect('equal')

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Diraspect_ratio_patches = mpatches.Patch(
                            color='white', label='aspect_ratio :' + str(aspect_ratio))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))

                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches], loc=2)

                        # =========================

                        cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                        # cax.set_label('Kelvin')
                        cax.set_aspect(20)
                        cax.set_ylabel('Density in code units', fontsize=10)
                        plt.colorbar(PL2, cax=cax)
                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Dynamic_Density_'
                                    + str(data_frame[n]['DirName']) + '_'
                                    + str(Dynamic_Loop * Standard_Orbit)
                                    + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close(fig)
                        Dynamic_Loop = Dynamic_Loop + dLoop

                    Dynamic_Loop = 0
                    dLoop = 1

                    print('================')
                    print('Making Dynamic_Temperature')
                    print('================')

                    if Calc_Temp == True:

                        while Dynamic_Loop <= len(Dynamic_Temperature) - 1:
                            # =========================
                            # =========================
                            # Temp Contour
                            # =========================
                            # =========================

                            Normalized_Temp = (
                                                      Dynamic_Temperature[Dynamic_Loop][:] - Init_Temp[:]) / Init_Temp[
                                                                                                             :]

                            fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                            fig.subplots_adjust(bottom=0.07, top=0.95)
                            ax1.contourf(
                                x2d, y2d, Dynamic_Temperature[Dynamic_Loop][:] / Init_Temp[:], ncolors)
                            ax1.set_aspect('equal')

                            DirMass_patches = mpatches.Patch(
                                color='white', label='q :' + str(DirMass))
                            DirEcc_patches = mpatches.Patch(
                                color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                            Dirsound_speed_patches = mpatches.Patch(
                                color='white', label='sound speed :' + str(sound_speed))
                            Diraspect_ratio_patches = mpatches.Patch(
                                color='white', label='aspect_ratio :' + str(aspect_ratio))
                            Dirinitial_pressure_patches = mpatches.Patch(
                                color='white', label='inital pressure :' + str(initial_pressure))
                            DirSigma_patches = mpatches.Patch(
                                color='white', label=r'$\Sigma$ :' + str(Sigma))
                            # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                            Dir_rsmooth_patches = mpatches.Patch(
                                color='white', label='potential smoothing :' + str(rsmooth))
                            Dir_gamma_patches = mpatches.Patch(
                                color='white', label=r'$\gamma$ :' + str(gamma))
                            Dir_alpha_patches = mpatches.Patch(
                                color='white', label=r'$\alpha$ :' + str(alpha))
                            Dir_beta_patches = mpatches.Patch(
                                color='white', label=r'$\beta$ :' + str(beta))

                            plt.legend(handles=[DirMass_patches,
                                                DirEcc_patches,
                                                Dirsound_speed_patches,
                                                Dirinitial_pressure_patches,
                                                DirSigma_patches,
                                                Dir_rsmooth_patches,
                                                Dir_gamma_patches,
                                                Dir_alpha_patches,
                                                Dir_beta_patches], loc=2)

                            # =========================

                            # cax=plt.axes([0.85,0.1,0.075,0.8])
                            # cax.set_label('Kelvin')
                            # cax.set_aspect(20)
                            # cax.set_ylabel('Temp in code units',fontsize=10)
                            plt.colorbar(PL2, ax=ax1)
                            plt.suptitle('t=' + str(data_frame[n]['ivar']))
                            plt.subplots_adjust(bottom=0.05, top=0.95)
                            plt.savefig('Dynamic_Temperature_'
                                        + str(data_frame[n]['DirName']) + '_'
                                        + str(Dynamic_Loop * Standard_Orbit)
                                        + '.png')
                            # too many figs may be open.
                            # close plots once done
                            plt.close(fig)

                            Dynamic_Loop = Dynamic_Loop + dLoop

                        Dynamic_Loop = 0
                        dLoop = 1
                    else:
                        print('================')
                        print('No Temperature to plot moving on (Dynamic)')
                        print('================')

                    print('================')
                    print('Making Dynamic_Shock')
                    print('================')

                    while Dynamic_Loop <= len(Dynamic_Shock) - 1:
                        # =========================
                        # =========================
                        # Shock Contour
                        # =========================
                        # =========================

                        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        ax1.contourf(
                            x2d, y2d, Dynamic_Shock[Dynamic_Loop][:], ncolors)
                        ax1.set_aspect('equal')

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Diraspect_ratio_patches = mpatches.Patch(
                            color='white', label='aspect_ratio :' + str(aspect_ratio))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))

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
                        # cax.set_label('Kelvin')
                        cax.set_aspect(20)
                        cax.set_ylabel('shock in code units', fontsize=10)
                        plt.colorbar(PL2, cax=cax)
                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Dynamic_Shock_'
                                    + str(data_frame[n]['DirName']) + '_'
                                    + str(Dynamic_Loop * Standard_Orbit)
                                    + '.png')
                        plt.close(fig)
                        Dynamic_Loop = Dynamic_Loop + dLoop

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Dynamic Contour Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Dynamic Contour Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingContour_Midplane(self, data_frame):

        Calc_Temp = self.Calc_Temp

        # =========================
        #
        #   Handle all plots relevant to contour plots
        #   ( top down disk, midplane , slice, etc)
        #   as well as physical profiles for said contour plot
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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

                    rho_fv = data_frame[n]['rho_fv']
                    shock_fv = data_frame[n]['shock_fv']
                    x_grid = data_frame[n]['x_grid']
                    y_grid = data_frame[n]['y_grid']

                    if Calc_Temp == True:
                        temp_fv = data_frame[n]['temp_fv']
                        Init_Temp = data_frame[n]['Init_Temp']
                        Normalized_Temp = (
                                                  temp_fv[:] - Init_Temp[:]) / Init_Temp[:]
                    else:
                        temp_fv = []
                        Init_Temp = []
                        Normalized_Temp = []
                    ncolors = 256

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    # =========================

                    # =========================
                    # =========================
                    # Density Midplane Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(20, 20))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    PL2 = ax1.contourf(x_grid, y_grid, rho_fv, ncolors)
                    # ax1.set_aspect('equal')
                    ax1.set_aspect(1 / ax1.get_data_ratio())
                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('Density in code units', fontsize=10)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Midplane_Contour_Density_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Temp Contour
                    # =========================
                    # =========================

                    if Calc_Temp == True:
                        #
                        #
                        #

                        Normalized_Temp = (
                                                  temp_fv[:] - Init_Temp[:]) / Init_Temp[:]

                        fig, (ax1) = plt.subplots(1, 1, figsize=(20, 20))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        ax1.contourf(x_grid, y_grid, Normalized_Temp, ncolors)
                        # ax1.set_aspect('equal')
                        ax1.set_aspect(1 / ax1.get_data_ratio())
                        # cax=plt.axes([0.85,0.1,0.075,0.8])
                        # cax.set_label('Kelvin')
                        # cax.set_aspect(20)
                        # cax.set_ylabel('Temp in code units',fontsize=10)
                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches], loc=2)
                        plt.colorbar(PL2, ax=ax1)
                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Midplane_Contour_Temp_' +
                                    str(data_frame[n]['DirName']) + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close(fig)
                    #
                    #
                    #
                    else:
                        print('================')
                        print('Not plotting temperture for midplane contour plots')
                        print('================')

                    # =========================
                    # =========================
                    # Shock Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(20, 20))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.contourf(x_grid, y_grid, shock_fv, ncolors)
                    # ax1.set_aspect('equal')
                    ax1.set_aspect(1 / ax1.get_data_ratio())
                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('shock in code units', fontsize=10)
                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Midplane_Contour_Shock_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Contour Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Contour Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingContour(self, data_frame):

        Calc_Temp = self.Calc_Temp

        # =========================
        #
        #   Handle all plots relevant to contour plots
        #   ( top down disk, midplane , slice, etc)
        #   as well as physical profiles for said contour plot
        #
        # =========================

        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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

                    rad_grid = data_frame[n]['rad_grid']

                    x2d = data_frame[n]['x2d']
                    y2d = data_frame[n]['y2d']
                    xrq2 = data_frame[n]['xrq2']
                    yrq2 = data_frame[n]['yrq2']

                    rho_fv = data_frame[n]['rho_fv']
                    temp_fv = data_frame[n]['temp_fv']
                    shock_fv = data_frame[n]['shock_fv']
                    Init_Temp = data_frame[n]['Init_Temp']

                    avgrho_fv = data_frame[n]['avgrho_fv']
                    avgtemp_fv = data_frame[n]['avgtemp_fv']
                    avgshock_fv = data_frame[n]['avgshock_fv']
                    Standard_Orbit = data_frame[n]['Standard_Orbit']

                    ncolors = 256

                    # =========================
                    # =========================
                    # Density Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    PL2 = ax1.contourf(x2d, y2d, rho_fv, ncolors)
                    ax1.set_aspect('equal')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('Density in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Density_' +
                                str(data_frame[n]['DirName']) + '_' + str(Standard_Orbit) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Temp Contour
                    # =========================
                    # =========================

                    if Calc_Temp == True:
                        #
                        #
                        #

                        Normalized_Temp = (
                                                  temp_fv[:] - Init_Temp[:]) / Init_Temp[:]

                        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        ax1.contourf(x2d, y2d, Normalized_Temp, ncolors)
                        ax1.set_aspect('equal')

                        # =========================
                        # legend handles
                        # =========================

                        # use this information when adding important parameters to the run

                        eccentricity = data_frame[n]['eccentricity']
                        sound_speed = data_frame[n]['cs']
                        aspect_ratio = data_frame[n]['aspect_ratio']
                        initial_pressure = data_frame[n]['initial_pressure']
                        Sigma = data_frame[n]['Sigma']
                        rsmooth = data_frame[n]['rsmooth']
                        gamma = data_frame[n]['gamma']
                        alpha = data_frame[n]['alpha']
                        beta = data_frame[n]['beta']

                        DirMass = data_frame[n]['par1']
                        DirEcc = (round(eccentricity[0], 1))

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))

                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches], loc=2)

                        # =========================

                        # cax=plt.axes([0.85,0.1,0.075,0.8])
                        # cax.set_label('Kelvin')
                        # cax.set_aspect(20)
                        plt.colorbar(PL2, ax=ax1)
                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Standard_Contour_Temp_' +
                                    str(data_frame[n]['DirName']) + '_' + str(Standard_Orbit) + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close(fig)
                    #
                    #
                    #
                    else:
                        print('================')
                        print('Not plotting temperture for Standard contour plots')
                        print('================')

                    # =========================
                    # =========================
                    # Shock Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.contourf(x2d, y2d, shock_fv, ncolors)
                    ax1.set_aspect('equal')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    cax.set_aspect(20)
                    cax.set_ylabel('shock in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Shock_' +
                                str(data_frame[n]['DirName']) + '_' + str(Standard_Orbit) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # -----------------------------------------------------------------------------------
                    # =========================
                    # Handle Contour Traces
                    # =========================
                    # -----------------------------------------------------------------------------------

                    # =========================
                    # =========================
                    # Density Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    PL2 = ax1.contourf(x2d, y2d, rho_fv, ncolors)
                    ax1.plot(xrq2, yrq2, color='orange')
                    ax1.set_aspect('equal')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    cax.set_aspect(20)
                    cax.set_ylabel('Density in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Trace_Density_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    if Calc_Temp == True:
                        #
                        #
                        #

                        # =========================
                        # =========================
                        # Temp Contour
                        # =========================
                        # =========================

                        Normalized_Temp = (
                                                  temp_fv[:] - Init_Temp[:]) / Init_Temp[:]

                        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        ax1.contourf(x2d, y2d, Normalized_Temp, ncolors)
                        ax1.plot(xrq2, yrq2, color='orange')
                        ax1.set_aspect('equal')

                        # =========================
                        # legend handles
                        # =========================

                        # use this information when adding important parameters to the run

                        eccentricity = data_frame[n]['eccentricity']
                        sound_speed = data_frame[n]['cs']
                        aspect_ratio = data_frame[n]['aspect_ratio']
                        initial_pressure = data_frame[n]['initial_pressure']
                        Sigma = data_frame[n]['Sigma']
                        EntropyIndex = data_frame[n]['EntropyIndex']
                        SpecificHeat = data_frame[n]['SpecificHeat']
                        rsmooth = data_frame[n]['rsmooth']
                        gamma = data_frame[n]['gamma']
                        Gamma0 = data_frame[n]['Gamma0']
                        alpha = data_frame[n]['alpha']
                        beta = data_frame[n]['beta']

                        DirMass = data_frame[n]['par1']
                        DirEcc = (round(eccentricity[0], 1))

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Diraspect_ratio_patches = mpatches.Patch(
                            color='white', label='aspect_ratio :' + str(aspect_ratio))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))

                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches], loc=2)

                        # =========================

                        # cax=plt.axes([0.85,0.1,0.075,0.8])
                        # cax.set_label('Kelvin')
                        # cax.set_aspect(20)
                        # cax.set_ylabel('Temp in code units',fontsize=10)
                        plt.colorbar(PL2, ax=ax1)
                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Standard_Contour_Trace_Temp_' +
                                    str(data_frame[n]['DirName']) + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close(fig)
                    #
                    #
                    #
                    else:
                        print('================')
                        print(
                            'Not plotting temperture for Standard contour plot Traces')
                        print('================')

                    # =========================
                    # =========================
                    # Shock Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.contourf(x2d, y2d, shock_fv, ncolors)
                    ax1.plot(xrq2, yrq2, color='orange')
                    ax1.set_aspect('equal')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    EntropyIndex = data_frame[n]['EntropyIndex']
                    SpecificHeat = data_frame[n]['SpecificHeat']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    Gamma0 = data_frame[n]['Gamma0']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Diraspect_ratio_patches = mpatches.Patch(
                        color='white', label='aspect_ratio :' + str(aspect_ratio))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    # DirEntropyIndex             = mpatches.Patch(color='white',label='Entropy Index :'+str(EntropyIndex))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('shock in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax)
                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Trace_Shock_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================

                    # -----------------------------------------------------------------------------------
                    # =========================
                    # Handle radial profiles
                    # =========================
                    # -----------------------------------------------------------------------------------

                    # =========================
                    # =========================
                    # Density Profile
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.plot(rad_grid, avgrho_fv)
                    ax1.set_xlabel('radius')
                    ax1.set_ylabel('average density')
                    ax1.set_aspect('auto')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Profile_Density_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Temp Profile
                    # =========================
                    # =========================

                    if Calc_Temp == True:
                        #
                        #
                        #

                        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                        fig.subplots_adjust(bottom=0.07, top=0.95)
                        ax1.plot(rad_grid, avgtemp_fv)
                        ax1.set_xlabel('radius')
                        ax1.set_ylabel('average temperature')
                        ax1.set_aspect('auto')

                        # =========================
                        # legend handles
                        # =========================

                        # use this information when adding important parameters to the run

                        eccentricity = data_frame[n]['eccentricity']
                        sound_speed = data_frame[n]['cs']
                        aspect_ratio = data_frame[n]['aspect_ratio']
                        initial_pressure = data_frame[n]['initial_pressure']
                        Sigma = data_frame[n]['Sigma']
                        rsmooth = data_frame[n]['rsmooth']
                        gamma = data_frame[n]['gamma']
                        alpha = data_frame[n]['alpha']
                        beta = data_frame[n]['beta']

                        DirMass = data_frame[n]['par1']
                        DirEcc = (round(eccentricity[0], 1))

                        DirMass_patches = mpatches.Patch(
                            color='white', label='q :' + str(DirMass))
                        DirEcc_patches = mpatches.Patch(
                            color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                        Dirsound_speed_patches = mpatches.Patch(
                            color='white', label='sound speed :' + str(sound_speed))
                        Dirinitial_pressure_patches = mpatches.Patch(
                            color='white', label='inital pressure :' + str(initial_pressure))
                        DirSigma_patches = mpatches.Patch(
                            color='white', label=r'$\Sigma$ :' + str(Sigma))
                        Dir_rsmooth_patches = mpatches.Patch(
                            color='white', label='potential smoothing :' + str(rsmooth))
                        Dir_gamma_patches = mpatches.Patch(
                            color='white', label=r'$\gamma$ :' + str(gamma))
                        Dir_alpha_patches = mpatches.Patch(
                            color='white', label=r'$\alpha$ :' + str(alpha))
                        Dir_beta_patches = mpatches.Patch(
                            color='white', label=r'$\beta$ :' + str(beta))

                        plt.legend(handles=[DirMass_patches,
                                            DirEcc_patches,
                                            Dirsound_speed_patches,
                                            Dirinitial_pressure_patches,
                                            DirSigma_patches,
                                            Dir_rsmooth_patches,
                                            Dir_gamma_patches,
                                            Dir_alpha_patches,
                                            Dir_beta_patches], loc=2)

                        # =========================

                        plt.suptitle('t=' + str(data_frame[n]['ivar']))
                        plt.subplots_adjust(bottom=0.05, top=0.95)
                        plt.savefig('Standard_Profile_Temperature_' +
                                    str(data_frame[n]['DirName']) + '.png')
                        # too many figs may be open.
                        # close plots once done
                        plt.close(fig)
                    #
                    #
                    #
                    else:
                        print('================')
                        print(
                            'Not plotting temperture for Standard contour plot radial Profile')
                        print('================')

                    # =========================
                    # =========================
                    # Shock Profile
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.plot(rad_grid, avgshock_fv)
                    ax1.set_xlabel('radius')
                    ax1.set_ylabel('average shock')
                    ax1.set_aspect('auto')

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.suptitle('t=' + str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Profile_Shock_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Contour Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Contour Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingEcc(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Ecc plot for ' +
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

                    time = data_frame[n]['time']
                    eccentricity = data_frame[n]['eccentricity']
                    CN_line = data_frame[n]['CN_line']

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))

                    ax1.plot(time, eccentricity[:len(
                        time)], label='Pencil-code')
                    ax1.plot(time, CN_line[:len(time)],
                             label='CN2008', linestyle=':')
                    plt.title(r'$\varepsilon$ Decay')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\varepsilon$')
                    plt.xlim([0, data_frame[n]['tmax']])
                    plt.tight_layout()
                    plt.grid(True)

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Standard_Ecc_Decay_' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Ecc Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Ecc Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def plotOrbitalEccentricity(self, data_frame, n):

        time = data_frame[n]['t']
        cut_off = data_frame[n]['indexTimeCutOff']
        time = time[:cut_off] / (2 * np.pi)

        # =========================
        # legend handles
        # =========================

        # use this information when adding important parameters to the run

        eccentricity = data_frame[n]['eccentricity']
        print('eccentricity is ',eccentricity)

        # =========================
        plt.tight_layout()
        plt.plot(time, eccentricity[:len(time)], color='orange')
        plt.xlabel('t')
        plt.ylabel('e')
        plt.title(r'$\varepsilon$' + ' '+'vs y')
        plt.grid(True)
        plt.show()
        plt.savefig('Standard_Orbit_Eccentricity' + str(data_frame[n]['DirName']) + '.png')
        # too many figs may be open.
        # close plots once done
        return None

    def plotOrbitalSemiMajor(self, data_frame, n):

        time = data_frame[n]['t']
        semi_major = data_frame[n]['semi_major']
        print('semi_major is ', semi_major)
        cut_off = data_frame[n]['indexTimeCutOff']
        time = time[:cut_off] / (2 * np.pi)
        print('time in semi_major is ',time)

        # =========================
        # legend handles
        # =========================

        # use this information when adding important parameters to the run

        eccentricity = data_frame[n]['eccentricity']

        # =========================
        fig1, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
        ax1.set_aspect('auto')
        plt.tight_layout()
        plt.plot(time, semi_major[:len(time)], color='orange')
        ax1.set_xlabel('t')
        ax1.set_ylabel('a')
        ax1.set_title('a vs t')
        plt.grid(True)
        plt.show()
        plt.savefig('Standard_Orbit_SemiMajor' + str(data_frame[n]['DirName']) + '.png')
        plt.close(fig1)
        # too many figs may be open.
        # close plots once done

    def plotOrbitalTopDown(self, data_frame, n):

        fig1, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
        xrq2 = data_frame[n]['xrq2'][:]
        yrq2 = data_frame[n]['yrq2'][:]

        print('xrq2', xrq2)
        print('yrq2', yrq2)

        # =========================
        plt.tight_layout()
        plt.plot(xrq2, yrq2)
        ax1.set_xlabel('xrq2')
        ax1.set_ylabel('yrq2')
        ax1.set_title('x vs y')
        plt.grid(True)
        ax1.set_aspect('equal')
        plt.savefig('Standard_Orbit_Info' + str(data_frame[n]['DirName']) + '.png')
        plt.close(fig1)
        # too many figs may be open.
        # close plots once done

    def pingOrbital(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Orbital plot for ' +
                              str(data_frame[n]['DirName']))
                        print('================')
                    except:
                        print('================')
                        print('Dir error ', str(data_frame[n]['DirName']))
                        print('================')
                        print(n)
                        print('================')
                        logging.basicConfig(filename='Pencil_Analysis.log',
                                            level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                                            datefmt=' %m/%d/%Y %I:%M:%S %p ')
                        logging.info('========================')
                        logging.info(n)
                        logging.info('========================')

                    xrq2 = data_frame[n]['xrq2'][:]
                    yrq2 = data_frame[n]['yrq2'][:]

                    print('xrq2 outside', xrq2)
                    print('yrq2 outside', yrq2)

                    self.plotOrbitalTopDown(data_frame, n)
                    self.plotOrbitalSemiMajor(data_frame, n)
                    self.plotOrbitalEccentricity(data_frame, n)

                    # use this information when adding important parameters to the run
                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Orbital Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Orbital Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingTempBoxPlot(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame) - 1:
                try:
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
                        print('Making Temp Box plot for ' +
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

                    tempdata = data_frame[n]['tempdata']
                    eccentricity = data_frame[n]['eccentricity']

                    fig1, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    plt.title('Temperature boxplots')
                    ax1.set_title(r'$\epsilon$=' +
                                  str(round(eccentricity[0], 1)))
                    ax1.boxplot(tempdata)
                    ax1.set_ylabel('Temperature')
                    ax1.set_xlabel(r'$\tau$*' + str(self.step))
                    plt.xscale('linear')
                    plt.grid(True)

                    # =========================
                    # legend handles
                    # =========================

                    # use this information when adding important parameters to the run

                    eccentricity = data_frame[n]['eccentricity']
                    sound_speed = data_frame[n]['cs']
                    aspect_ratio = data_frame[n]['aspect_ratio']
                    initial_pressure = data_frame[n]['initial_pressure']
                    Sigma = data_frame[n]['Sigma']
                    rsmooth = data_frame[n]['rsmooth']
                    gamma = data_frame[n]['gamma']
                    alpha = data_frame[n]['alpha']
                    beta = data_frame[n]['beta']

                    DirMass = data_frame[n]['par1']
                    DirEcc = (round(eccentricity[0], 1))

                    DirMass_patches = mpatches.Patch(
                        color='white', label='q :' + str(DirMass))
                    DirEcc_patches = mpatches.Patch(
                        color='white', label=r'$\varepsilon$ :' + str(DirEcc))
                    Dirsound_speed_patches = mpatches.Patch(
                        color='white', label='sound speed :' + str(sound_speed))
                    Dirinitial_pressure_patches = mpatches.Patch(
                        color='white', label='inital pressure :' + str(initial_pressure))
                    DirSigma_patches = mpatches.Patch(
                        color='white', label=r'$\Sigma$ :' + str(Sigma))
                    Dir_rsmooth_patches = mpatches.Patch(
                        color='white', label='potential smoothing :' + str(rsmooth))
                    Dir_gamma_patches = mpatches.Patch(
                        color='white', label=r'$\gamma$ :' + str(gamma))
                    Dir_alpha_patches = mpatches.Patch(
                        color='white', label=r'$\alpha$ :' + str(alpha))
                    Dir_beta_patches = mpatches.Patch(
                        color='white', label=r'$\beta$ :' + str(beta))

                    plt.legend(handles=[DirMass_patches,
                                        DirEcc_patches,
                                        Dirsound_speed_patches,
                                        Dirinitial_pressure_patches,
                                        DirSigma_patches,
                                        Dir_rsmooth_patches,
                                        Dir_gamma_patches,
                                        Dir_alpha_patches,
                                        Dir_beta_patches], loc=2)

                    # =========================

                    plt.savefig('Standard_Temp_Box' +
                                str(data_frame[n]['DirName']) + '.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig1)
                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Temp Data found to plot')
                    print('================')
                    traceback.print_exc()
                    logging.basicConfig(filename='Pencil_Analysis.log',
                                        level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                                        datefmt=' %m/%d/%Y %I:%M:%S %p ')
                    logging.info('========================')
                    logging.info('Temp Box Plot Data error')
                    logging.info('========================')
                    logging.info('========================')
                    logging.info('========================')
                    logging.info('========================')
                    logging.info('========================')
                    os.chdir('..')
                n = n + dn
        except:
            print('================')
            print('ping Temp Box Loop error')
            print('================')
            traceback.print_exc()
            logging.basicConfig(filename='Pencil_Analysis.log',
                                level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                                datefmt=' %m/%d/%Y %I:%M:%S %p ')
            logging.info('========================')
            logging.info('ping temp box plot error')
            logging.info('========================')

            os.chdir('..')
            return False
