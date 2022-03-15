
# se:
#
#   Go through each sub directory in the cwd
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


class Pencil_Post_Analysis(object):

    def __init__(self,
                 data_functions,
                 Orbit=None,
                 step=None,
                 TempSigma=False,
                 dir_run_list=None,
                 var_dir_list=None,
                 temp_dir_list=None,
                 Calc_DTTm=None):

        self.Orbit = Orbit
        self.data_functions = data_functions
        self.dir_run_list = dir_run_list
        self.var_dir_list = var_dir_list
        self.temp_dir_list = temp_dir_list
        self.step = step
        self.TempSigma = TempSigma

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

        var_dir_list = []
        dir_run_list = next(os.walk('.'))[1]

        i = 0
        di = 1

        while i <= len(dir_run_list)-1:
            var_dir_list.append(data_functions.grepDATA(
                dir_run_list[i], Orbit, Calc_Temp=True, Calc_Density=True, step=step))
            data_functions.pingGit(dir_run_list[i])
            data_functions.pingTemp(dir_run_list[i])
            i = i+di

        return var_dir_list

    def Make_TempVars(self, data_frame):

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
            while i <= len(data_frame)-1:
                try:
                    if data_frame[i]['tempdata'] == None:
                        i = i+di
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
                        i = i+di
                except:
                    traceback.print_exc()
                    print('================')
                    print('No Temp Data found to plot')
                    print('================')
                i = i+di

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

            while f <= len(temp_frame)-1:
                try:
                    if math.isnan(np.std(temp_frame[f][:])) == False:
                        temp_distribution.append(np.std(temp_frame[f][:]))
                        if ecc_frame == None or name_frame == None:
                            f = f+df
                        else:
                            temp_ecc.append(ecc_frame[f])
                            temp_name.append(name_frame[f])
                            f = f+df
                    else:
                        f = f+df
                except:
                    print('f failed at ' + str(f))
                    f = f+df
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
            while n <= len(data_frame)-1:
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

                    plt.savefig('Standard_Precession_Check_Post' +
                                str(data_frame[n]['DirName'])+'.png')
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
                n = n+dn
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
            while n <= len(data_frame)-1:
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

                    plt.savefig('Standard_Precession_Check_10_Orbits_Post' +
                                str(data_frame[n]['DirName'])+'.png')
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
                n = n+dn
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
            while n <= len(data_frame)-1:
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

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.plot(data_frame[n]['time'], data_frame[n]
                             ['torqint'], '--', label='Inner')
                    ax1.plot(data_frame[n]['time'], data_frame[n]
                             ['torqext'], '--', label='Outer')
                    ax1.plot(data_frame[n]['time'], data_frame[n]
                             ['torqtotal'], label='Total')
                    ax1.plot(data_frame[n]['time'], data_frame[n]['tanaka_array'],
                             linestyle=':', label='Tanaka et  al. 2002')
                    ax1.set_aspect('auto')
                    plt.xlim([0, data_frame[n]['tmax']])
                    plt.tight_layout()
                    plt.grid(True)

                    plt.title(r'Torques')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\Gamma$')
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Torque_Post' +
                                str(data_frame[n]['DirName'])+'.png')
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
                n = n+dn
        except:
            print('================')
            print('ping Torque Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingContour(self, data_frame):

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
            while n <= len(data_frame)-1:
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
                    rho_fv = data_frame[n]['rho_fv']
                    temp_fv = data_frame[n]['temp_fv']
                    shock_fv = data_frame[n]['shock_fv']
                    Init_Temp = data_frame[n]['Init_Temp']

                    avgrho_fv = data_frame[n]['avgrho_fv']
                    avgtemp_fv = data_frame[n]['avgtemp_fv']
                    avgshock_fv = data_frame[n]['avgshock_fv']

                    ncolors = 256

                    # =========================
                    # =========================
                    # Density Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(15, 15))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    PL2 = ax1.contourf(x2d, y2d, rho_fv, ncolors)
                    ax1.set_aspect('auto')

                    cax = plt.axes([0.90, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('Density in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax, pad=0.2)
                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Density_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Temp Contour
                    # =========================
                    # =========================

                    Normalized_Temp = temp_fv[:]-Init_Temp[:]

                    fig, (ax1) = plt.subplots(1, 1, figsize=(15, 15))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.contourf(x2d, y2d, Normalized_Temp, ncolors)
                    ax1.set_aspect('auto')

                    cax = plt.axes([0.90, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('Temp in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax, pad=0.2)
                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Temp_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Shock Contour
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(15, 15))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.contourf(x2d, y2d, shock_fv, ncolors)
                    ax1.set_aspect('auto')

                    cax = plt.axes([0.90, 0.1, 0.075, 0.8])
                    # cax.set_label('Kelvin')
                    cax.set_aspect(20)
                    cax.set_ylabel('shock in code units', fontsize=10)
                    plt.colorbar(PL2, cax=cax, pad=0.2)
                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Contour_Shock_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

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

                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Profile_Density_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    # =========================
                    # =========================
                    # Temp Profile
                    # =========================
                    # =========================

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
                    fig.subplots_adjust(bottom=0.07, top=0.95)
                    ax1.plot(rad_grid, avgtemp_fv)
                    ax1.set_xlabel('radius')
                    ax1.set_ylabel('average temperature')
                    ax1.set_aspect('auto')

                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Profile_Temperature_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

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

                    plt.suptitle('t='+str(data_frame[n]['ivar']))
                    plt.subplots_adjust(bottom=0.05, top=0.95)
                    plt.savefig('Standard_Profile_Shock_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    os.system('git add *.png')
                    os.chdir('..')
                except:
                    print('================')
                    print('No Temp Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n+dn
        except:
            print('================')
            print('ping Density Contour Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingEcc(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame)-1:
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
                    CN2008_te = data_frame[n]['CN2008_te']

                    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))

                    ax1.plot(time, eccentricity[:len(
                        time)], label='Pencil-code')
                    # ax1.plot(time,CN2008_te[:len(time)],label='CN2008',linestyle=':')
                    plt.title(r'$\varepsilon$ Decay')
                    plt.xlabel(r'$t/T_0$')
                    plt.ylabel(r'$\varepsilon$')
                    # plt.yscale('log')
                    plt.xlim([0, data_frame[n]['tmax']])
                    plt.tight_layout()
                    plt.grid(True)

                    plt.savefig('Standard_Ecc_Decay_Post' +
                                str(data_frame[n]['DirName'])+'.png')
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
                n = n+dn
        except:
            print('================')
            print('ping Ecc Loop error')
            print('================')
            traceback.print_exc()
            os.chdir('..')
            return False

    def pingOrbital(self, data_frame):
        try:
            n = 0
            dn = 1
            while n <= len(data_frame)-1:
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
                    semi_major = data_frame[n]['semi_major']
                    xrq2 = data_frame[n]['xrq2']
                    yrq2 = data_frame[n]['yrq2']

                    fig, ((ax2)) = plt.subplots(1, 1, figsize=(30, 5))
                    fig.tight_layout()

                    plt.title('Orbital path')
                    ax2.set_aspect('equal')
                    # ax3.set_aspect('auto')

                    ax2.set_title('x vs y')
                    #ax3.set_title('a vs t')

                    ax2.plot(xrq2, yrq2)
                    ax2.set_xlabel('x')
                    ax2.set_ylabel('y')
                    # ax3.plot(time,semi_major[:len(time)],color='orange')
                    # ax3.set_xlabel('t')
                    # ax3.set_ylabel('a')

                    ax2.grid(True)
                    # ax3.grid(True)

                    plt.subplots_adjust(wspace=0.18, bottom=0.12)

                    plt.savefig('Standard_Orbit_Info_Post' +
                                str(data_frame[n]['DirName'])+'.png')
                    # too many figs may be open.
                    # close plots once done
                    plt.close(fig)

                    plt.plot(time, semi_major[:len(time)], color='orange')
                    plt.xlabel('t')
                    plt.ylabel('a')
                    plt.grid(True)
                    plt.savefig('Standard_Orbit_SemiMajor' +
                                str(data_frame[n]['DirName'])+'.png')
                    plt.close(fig)

                    os.system('git add *.png')
                    os.chdir('..')

                except:
                    print('================')
                    print('No Orbital Data found to plot')
                    print('================')
                    traceback.print_exc()
                    os.chdir('..')
                n = n+dn
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
            while n <= len(data_frame)-1:
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
                    ax1.set_xlabel(r'$\tau$*'+str(self.step))
                    plt.xscale('linear')
                    plt.grid(True)

                    plt.savefig('Standard_Temp_Box_Post' +
                                str(data_frame[n]['DirName'])+'.png')
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
                    # logging.info(data_frame[n]['DirName'])
                    logging.info('========================')
                    logging.info('========================')
                    # logging.info(str(traceback.print_exc()))
                    logging.info('========================')
                    os.chdir('..')
                n = n+dn
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
