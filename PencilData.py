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


def func(X, A, B):
    # function for power law fitting
    # basic function to be used outside of the class
    # has issues with scipy if used inside the class

    return A * np.exp(-B * X)


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

        print('================')
        print('getting vars')
        print('================')

        self.Directory_Path = Directory_Path
        self.Orbit = Orbit
        self.Calc_Temp = Calc_Temp
        self.Calc_Density = Calc_Density
        self.Orbit_standard = Orbit_standard
        self.Calc_Energy = Calc_Energy
        self.Calc_OEnergy = Calc_OEnergy
        self.Calc_Dynamics = Calc_Dynamics
        self.Calc_Rates_Energy = Calc_Rates_Energy
        self.ToomreQ = Calc_ToomreQ
        self.MaxOrbits = MaxOrbits

        Standard_Orbit = Orbit_standard

        self.step = step

        # ===========================================================================
        #
        #   get all var data from pencil datacubes
        #   for the given directory path
        #   get all fit parameters not outputted in the print.in
        #   calculate fit parameters for analysis
        #   calculate any extra variables
        #   include conditionals for cases where variables are Nan or problematic
        #   regardless of whats being asked the vars will be read from the qvar.namelist
        #   if a var is missing its going to be reported
        #
        #   get time series and read_vars as a python dict for analysis
        #
        # ===========================================================================
        #
        #   DEV LOG
        #
        #   12/27/2018 added time series vars
        #
        #   1/1/2019   generated hermite polynomials for TW in 2D
        #
        #   1/2/2019   started to add generation of fourier-hermite components
        #              added integrator methods
        #
        #   1/3/2019   added generation of first couple n(m,n,q)(enthalpy) fourier-hermite components
        #              further testing needs to be done to determine accuracy of test spiral waves
        #              finished debugging, tempboxplots, torque, eccdecay, and orbital scripts
        #
        #   1/16/2019  updated var looping to include density data
        #              generated Mu_coef arrays
        #
        #   2/15/2019  added TTm from print as a new parameter
        #
        #   2/22/2019  Updated several new parameters to the pencil data script
        #              updated the script for sigma temp calculation
        #              added new correct true_angle and true_anamoly
        #              added dTTm/dt as a parameter
        #              fixed contour plots for sigma temp
        #              Investigation of perihelion is top priority as of now
        #
        #   2/29/2019  Perihelion plots have been fixed
        #
        #   3/15/2019  contours plots changed for dT/dt
        #              also included post analysis scripts for presentations
        #              file known as PencilPostAnalysis.py
        #
        #
        #   4/30/2019  including energy analysis vars
        #              need to include orbit standard vars for measuring long term effects
        #              spiral waves still very tricky to talk about. How energy is being transported in the disk is not well known
        #              soon to include luminous torques/heating torques. New analysis will have to be considered for such a case
        #              since runs are not 3D, very little physics is needed to worry about.
        #              need to include optimized write to file to lower script time.
        #              migration times need to be updated.
        #
        #   5/1/2019   continued working on energy plots
        #
        #   5/29/2019  energy plots on hold as of now. Will work them out on 6/4/2019.
        #
        #   6/3/2019   added correct CN08 fit. This has been a problem recently.
        #	       updated script to include optimzation. Now everything does not have to be calculated at once
        #
        #   6/4/2019   added mean temprature vs time calculation
        #	           investigating optimization schema to fix runtime issues
        #	           optimization schema on hold until more recent problems are dealt with.
        #
        #   6/5/2019   fixed energy plot
        #
        #   6/11/2019  added corrected energy plots
        #	       convolved energy plot data
        #
        #   6/13/2019  applied standard orbit to pencil data
        #	           it shall now be used for dynamic contour plots
        #	           and midplane slices.
        #
        #   6/17/2019  added energy plot fitting
        #              added energy plot rates
        #
        #   6/18/2019  continued to build on fitting procedure
        #
        #   6/20/2019  fixed error in fitting procedure
        #	           fixed bug in energy plots
        #
        #   6/24/2019  fixed an issue with Energy plotting not displaying properly
        #              fixed an issue with Dynamic plotting not working at all
        #
        #   6/27/2019  fixed an issue with ivars with no orbits
        #
        #   7/3/2019   updated googledrives/github to have latest version of python scripts
        #
        #   7/12/2019  updated PencilAnalysis.py to have TTm rates for comparsion to runs with cooling
        #              fixed issue with PencilData crashing when running routine in a Isothermal run
        #              (no conditional was set for case of isothermal)
        #              fixed Dynamics not working in isothermal runs
        #
        #   7/27/2019  added orbital energy for eccentric cases
        #
        #   11/7/2019  Toomre Q parameter is important for discussing disk stability for AGN orbits. Its
        #              gonna to be useful to measure the Toomre Q parameter anyways. You should be able to
        #              use this in the case of self gravity or no self gravity for pencil-code
        #              Added Check to make sure Qvars are not missing. Sometimes the number of orbits is
        #              missing, leading to problems when plotting. This could very well be an issue with stampede
        #              removing some data on the scratch
        #
        #   11/8/2019  Finished ToomreQ parameter issues
        #
        #   12/6/2019  Used pep8 on all files, should be compatible with python3 now (tab error)
        #              fixed some issues involving basic scripts (python3 does not want to save figures
        #              for whatever reason).
        #
        #   5/4/2020   Updated python scripts for latest use of pencil code.
        #              Changed import from pencil to pencil_old
        #              Noticed issue with non-entropy runs for Calculating Toomre Q
        #              Toomre Q requires gravitational constant values. Since this is NOT set in pc.read_params
        #              The gravitational constant has to be manually checked.
        #
        #   5/19/2020  Fixed Toomre Q calculations. It should now display a reasonable Toomre Q calc("should")
        #              Updated kernel size. It should be N=850 by default.
        #
        #   6/5/2020   Simone's birthday is tommorrow. Happy birthday love. Added Calc Orbital Energy logical
        # 	            This is seperate from all other energy calculations and requires a consitency check.
        #	            Fixed any issues with not properly using Calc Orbital Energy Logical
        #	            It is under Calc_OEnergy. So it is changed accordingly.
        #               Fixed synatx issues. Sigmap now properly calculated in case of
        #               rho0 not being 1
        #
        #  6/6/2020    Added relative temperature calculations to all plotting
        #              relative vales should be calculated using the following formula
        #
        #              (value-initial_value)/(initial_value)
        #
        #              worked on ToomreQ colorbar fixes
        #
        # ==========================================================================================================

        print(os.getcwd())
        os.chdir(self.Directory_Path)
        print(os.getcwd())

        try:
            #
            #
            #
            # before anything is done, check to make sure that there are no missing vars.
            #
            #
            #
            #
            try:
                os.chdir('data/proc0')

                var = open('varN.list', 'r')
                vline = var.readlines()

                temp = []
                store = []
                missing = []

                for line in vline:
                    temp.append(line)

                for index in range(len(temp)):
                    temp_string = re.search(
                        '(?<=VAR)[0-9]{0,100}', str(temp[index]))
                    store.append(temp_string.group(0))

                i = 0
                di = 1
                count = 0

                for index in range(len(store)):
                    if i == int(store[index]):
                        i = i + di
                    else:
                        while i <= int(store[index]):
                            if i == int(store[index]):
                                i = i + di
                            else:
                                missing.append('missing VAR' + str(i))
                                count = count + 1
                                i = i + di
                    if i == len(store):
                        missing.append('end with :' + str(count) + ' missing')

                os.chdir('..')
                os.chdir('..')
            #
            #
            #
            #
            except:
                #
                #
                #
                print('=========')
                print('no proc0')
                print('=========')

            # ======================================
            #   pc.read_var()
            # ======================================
            #
            #   get intial read var elements
            #

            # ======================================
            # pc.read_ts()
            # ======================================

            ts = pc.read_ts()
            t = ts.t
            N = 850
            tmax = t.max()

            # DEBUG CALC ENERGY

            # quick set up for intervals of time
            # use Orbit_Len to properly get the length of the interval
            # for Orbits from the time series

            time = t/(2 * np.pi)
            len_t = len(ts.t)

            Max_Orbits = np.round((time[len(time) - 1]))

            n = 1
            dn = 1

            i = 0
            di = 1

            Orbit_Len = []

            while i <= len(time) - 1:
                if np.round(n * 2.0 * np.pi) == np.round(t[i]):
                    Orbit_Len.append(len(t[:i]))
                    n = n + dn
                    i = i + di
                else:
                    i = i + di

            # if Calc_Energy == False:
            #	Orbit_Len = None

            # ======================================
            # pc.read_pararm()
            # ======================================

            par = pc.read_param()
            h = par.cs0
            if (par.iprimary == 1):
                q = par.pmass[1]
            else:
                q = par.pmass[0]

            par1 = par.pmass[1]
            par2 = par.pmass[0]
            gamma = par.gamma
            Gamma0 = (q / h) ** 2
            alpha = par.density_power_law
            beta = par.temperature_power_law
            kernel = np.ones((N,)) / N
            time = np.convolve(t, kernel, mode='valid')
            torqint = np.convolve(ts.torqint_2, kernel, mode='valid')
            torqext = np.convolve(ts.torqext_2, kernel, mode='valid')
            torqtotal = torqext[:] + torqint[:]
            rsmooth = par.r_smooth[1]
            gravC = par.g0
            EntropyIndex = beta - ((gamma - 1) * alpha)
            SpecificHeat = par.cp
            Sigma = par.rho0
            Mstar = par.pmass[0]
            cs = par.cs0  # sound speed

            # ======================================

            radius = ts.xq2
            Omega = ts.yq2
            LinearVelocity = ts.vxq2
            AngularVelocity = ts.vyq2

            # ======================================

            # if temperature is too be calculated

            TTm = []
            TTm_rate = []
            GlobalTemp_Mean = []
            DGTemp_Mean = []
            GTM_Sigma = []
            GTM_Sigma_O = []

            if Calc_Temp == True:
                TTm = ts.TTm
                TTm_rate = np.mean(np.gradient(TTm[:]))
                GlobalTemp_Mean = ts.TTm

                D = 0
                dD = 1

                DGlobalTemp_Mean = []

                while D <= len(GlobalTemp_Mean) - 2:
                    DGlobalTemp_Mean.append(
                        GlobalTemp_Mean[D + 1] - GlobalTemp_Mean[D])
                    D = D + dD

                DGTemp_Mean = sum(DGlobalTemp_Mean) / len(DGlobalTemp_Mean)
                GTM_Sigma = np.std(GlobalTemp_Mean)

                if Orbit == None:
                    GTM_Sigma_O = 0
                else:
                    GTM_Sigma_O = np.std(GlobalTemp_Mean[:int(2 * np.pi * Orbit)])

            else:

                # do not calculate plots with temperature

                print('----------------')
                print('NOT calculating temperature of the disk. Please dont plot it')
                print('----------------')

                TTm = []
                TTm_rate = []
                GlobalTemp_MEan = []
                DGTemp_Mean = []
                GTM_Sigma = []
                GTM_Sigma = []

            # ======================================

            v2 = LinearVelocity ** 2 + AngularVelocity ** 2
            semi_major = 1. / (2 / radius - v2)
            DArclength = radius ** 2 * (AngularVelocity / radius)
            ep1 = (DArclength ** 2) / semi_major
            eccentricity = (1 - ep1) ** 0.5
            Omegap = 1. / semi_major ** 1.5
            Kepler_F = np.sqrt(gravC * Mstar / radius)

            ScaleHeight = cs / Kepler_F

            Hill_Radius = semi_major * ((q / 3.0) ** (1.0 / 3.0))
            xrq2 = radius * np.cos(Omega)
            yrq2 = radius * np.sin(Omega)

            ecc_int = par.eccentricity

            aspect_ratio = ScaleHeight/radius

            Sigmap = Sigma * radius ** -alpha

            eh = eccentricity / ScaleHeight

            mstar = 1
            twave = (mstar * ScaleHeight ** 4) / \
                    (q * Sigmap * semi_major ** 2 * Omegap)

            te = (twave / 0.780) * (1 - 0.14 * eh ** 2 + 0.06 * eh ** 3)
            CN_line = ecc_int * np.exp(-t / te)
            ecc = eccentricity
            ecc_int = par.eccentricity
            ecc_rate = np.gradient(eccentricity)
            ecc_ang = par1 * np.sqrt(semi_major) * np.sqrt(1 - ecc ** 2)

            print('len of ecc is ',str(len(ecc)))
            print('t.max() is ',t.max())

            # for convience in calculating the longitude of the perihelion
            # i introduced ecc as variable and sma as a variable
            # these should be set to none after they have been used for convience

            ecc = eccentricity
            sma = semi_major

            if ecc.all() == 0:
                fixed_ecc = []
                true_angle = []
                # r=1
                i = 0
                di = 1
                while i <= len(ecc) - 1:
                    fixed_ecc.append(ecc[i] + 0.000000000000009)
                    true_angle.append(
                        (1.0 / fixed_ecc[i]) * ((sma[i] * (1.0 - fixed_ecc[i] ** 2.0) / radius[i]) - 1.0))

                    i = i + di
            else:
                true_angle = (1.0 / ecc) * ((sma * (1.0 - ecc ** 2.0) / radius) - 1.0)

            true_anomaly = ((np.arccos(true_angle)) * 180 / np.pi) - 180
            long_perihelion = np.absolute(
                Omega * (180 / np.pi)) + true_anomaly  # arguement of perihelion
            # turn the arguement of perihelion into radians so that we can compare it to other papers
            long_perihelion = long_perihelion * (np.pi / 180)

            indexTimeCutOff = 0;
            indexTimeCutOffLarge = 0;

            for i in range(len(ecc)):
                if ecc_int != 0:
                    if ecc[i] <= 0.01:
                        indexTimeCutOff = i
                        break

            timeCutOff = t[indexTimeCutOff]/(np.pi*2)
            timeCutOffLarge = t[indexTimeCutOffLarge]/(np.pi*2)

            print('number of snapshots',t.max()/(np.pi*2))
            print('indexTimeCutOff',indexTimeCutOff)
            print('t[indexTimeCutOff]',t[indexTimeCutOff])
            print('t[indexTimeCutOff]/(np.pi*2)',t[indexTimeCutOff]/(np.pi*2))
            print('ecc[indexTimeCutOff]',ecc[indexTimeCutOff])
            print('timeCutOff: ', timeCutOff)

            ecc = None
            sma = None

            # ------------------------------------------------------
            # ======================================
            #   pc.read_var()
            # ======================================

            #   get intial read var elements for disk grid elements

            #
            #  NOTE: conditionals for temperature
            #  need to be applied every time
            #  you are asking for disk temperature.
            #

            ff = []
            rad = []
            theta = []
            rad_grid = []
            rad2d = []
            theta2d = []
            y2d = []
            x2d = []
            x_grid = []
            y_grid = []
            Init_Temp = []
            rho = []

            try:
                if Calc_Temp == True:
                    print('entering Calc_temp')
                    ff = pc.read_var(trimall=True, ivar=0, magic=['TT'], quiet=True)
                    rad = ff.x
                    theta = ff.y
                    rad_grid = ff.x
                    rad2d, theta2d = np.meshgrid(rad, theta)
                    x2d = rad2d * np.cos(theta2d)
                    y2d = rad2d * np.sin(theta2d)
                    x_grid = ff.x
                    y_grid = ff.y
                    Init_Temp = ff.TT
                    rho = ff.rho
                else:
                    print('entering no Calc_temp')
                    ff = pc.read_var(trimall=True, ivar=0, quiet=True)
                    rad = ff.x
                    theta = ff.y
                    rad_grid = ff.x
                    rad2d, theta2d = np.meshgrid(rad, theta)
                    x2d = rad2d * np.cos(theta2d)
                    y2d = rad2d * np.sin(theta2d)
                    x_grid = ff.x
                    y_grid = ff.y
                    rho = ff.rho
            except:
                print('================')
                print('No var data to be found')
                print('================')
                traceback.print_exc()

            # get var elements to read for plots
            # let this setup be for looking at only one orbit snap shot for the disk
            # initialize vars so that there is no issue if not calculated

            fv = []
            rho_fv = []
            temp_fv = []
            shock_fv = []
            avgrho_fv = []
            avgtemp_fv = []
            avgshock_fv = []

            try:
                if Calc_Temp == True:
                    ivar = Orbit
                    fv = pc.read_var(trimall=True, ivar=ivar, magic=['TT'], quiet=True)
                    rho_fv = fv.rho
                    temp_fv = fv.TT
                    shock_fv = fv.shock
                    avgrho_fv = np.transpose(rho_fv)
                    avgtemp_fv = np.transpose(temp_fv)
                    avgshock_fv = np.transpose(shock_fv)
                    avgrho_fv = np.mean(avgrho_fv, axis=1)
                    avgtemp_fv = np.mean(avgtemp_fv, axis=1)
                    avgshock_fv = np.mean(avgshock_fv, axis=1)
                    # print(rho_fv)
                else:
                    ivar = Orbit
                    fv = pc.read_var(trimall=True, ivar=ivar, quiet=True)
                    rho_fv = fv.rho
                    shock_fv = fv.shock
                    avgrho_fv = np.transpose(rho_fv)
                    avgshock_fv = np.transpose(shock_fv)
                    avgrho_fv = np.mean(avgrho_fv, axis=1)
                    avgshock_fv = np.mean(avgshock_fv, axis=1)
                    # print(rho_fv)
            except:
                try:
                    if Calc_Temp == True:
                        ivar = 20
                        fv = pc.read_var(trimall=True, ivar=ivar, magic=['TT'], quiet=True)
                        rho_fv = fv.rho
                        temp_fv = fv.TT
                        shock_fv = fv.shock
                        avgrho_fv = np.transpose(rho_fv)
                        avgtemp_fv = np.transpose(temp_fv)
                        avgshock_fv = np.transpose(shock_fv)
                        avgrho_fv = np.mean(avgrho_fv, axis=1)
                        avgtemp_fv = np.mean(avgtemp_fv, axis=1)
                        avgshock_fv = np.mean(avgshock_fv, axis=1)
                    else:
                        ivar = 20
                        fv = pc.read_var(trimall=True, ivar=ivar, quiet=True)
                        rho_fv = fv.rho
                        shock_fv = fv.shock
                        avgrho_fv = np.transpose(rho_fv)
                        avgshock_fv = np.transpose(shock_fv)
                        avgrho_fv = np.mean(avgrho_fv, axis=1)
                        avgshock_fv = np.mean(avgshock_fv, axis=1)
                except:
                    rho_fv = None
                    temp_fv = None
                    shock_fv = None
                    avgrho_fv = None
                    avgtemp_fv = None
                    avgshock_fv = None

            # Generate Temp Data Arrays for maximum amount of Orbit
            # for the given run

            i = 0
            di = step

            f = 0
            df = 1

            placedata = []
            tempdata = []

            if Calc_Temp == False:
                tempdata = None
            else:
                try:
                    print('Entering desnity ivar loop')
                    if step != None:
                        while i <= MaxOrbits - step:
                            ff = pc.read_var(trimall=True, ivar=i, magic=["TT"], quiet=True)
                            while f <= len(ff.TT) - 1:
                                for value in ff.TT[f]:
                                    placedata.append(value)
                                f = f + df
                            f = 0
                            tempdata.append(placedata)
                            placedata = []
                            i = i + di
                    else:
                        print('step is none')
                except:
                    print('================')
                    print('No tempdata to be found')
                    print('================')
                    traceback.print_exc()
                    tempdata = None

            # Generate Density Data Arrays for maximum amount of Orbit
            # for the given run

            i = 0
            di = step

            f = 0
            df = 1

            placedata = []
            densitydata = []

            if Calc_Density == False:
                densitydata = None
            else:
                try:
                    print('Entering temp ivar loop')
                    if step != None:
                        while i <= MaxOrbits - step:
                            ff = pc.read_var(trimall=True, ivar=i, quiet=True)
                            while f <= len(ff.rho) - 1:
                                for value in ff.rho[f]:
                                    placedata.append(value)
                                f = f + df
                            f = 0
                            densitydata.append(placedata)
                            placedata = []
                            i = i + di
                    else:
                        print('step is none')
                except:
                    print('================')
                    print('No densityData to be found')
                    print('================')
                    traceback.print_exc()
                    densitydata = None

            # ======================================
            # |
            # |
            # ======================================
            if Standard_Orbit is None or Calc_Dynamics is False:
                # ======================================

                # keep vars used in Dynamic Loop here
                # these are going to be two dimensional arrays

                Dynamic = None
                Standard_Orbit = None

                Dynamic_Density = None
                Dynamic_Temperature = None
                Dynamic_Shock = None

            else:

                # declare int to be standard orbit

                Dynamic = int(Standard_Orbit)
                Dynamic_Density = []
                Dynamic_Temperature = []
                Dynamic_Shock = []

                # vars for dynamic contour plots
                # vars for midplane slices
                try:
                    while Standard_Orbit <= Max_Orbits - Dynamic:
                        if Calc_Temp == True:
                            fv = pc.read_var(
                                trimall=True, ivar=Standard_Orbit, magic=['TT'], quiet=True)
                            rho_fv = fv.rho
                            temp_fv = fv.TT
                            shock_fv = fv.shock

                            Dynamic_Density.append(rho_fv)
                            Dynamic_Temperature.append(temp_fv)
                            Dynamic_Shock.append(shock_fv)

                        else:
                            fv = pc.read_var(trimall=True, ivar=Standard_Orbit, quiet=True)
                            rho_fv = fv.rho
                            shock_fv = fv.shock

                            Dynamic_Density.append(rho_fv)
                            Dynamic_Shock.append(shock_fv)

                        Standard_Orbit = Standard_Orbit + Dynamic
                except:
                    print('================')
                    print('Dynamic Loop error')
                    print('================')
                    traceback.print_exc()

            # ======================================
            # |
            # |
            # ======================================
            if Calc_Energy == True:
                # ======================================

                # intialize all vars for energy array
                # be careful in the case of isothermal disks

                ecc = eccentricity

                KE_Sum = []
                UE_Sum = []
                UINT_Sum = []
                OE_Sum = []
                Total_Disk_Energy = []
                Total_Disk_Energy_Avg = []
                Total_Disk_Grad = []
                Total_Disk_Grad_Avg = []
                KE_Sum_Avg = []
                UE_Sum_Avg = []
                OE_Sum_Avg = []
                UINT_Sum_Avg = []
                KE_Grad = []
                UE_Grad = []
                OE_Grad = []
                UINT_Grad = []
                KE_Grad_Avg = []
                UE_Grad_Avg = []
                OE_Grad_Avg = []
                UINT_Grad_Avg = []
                KE_fit = []
                UE_fit = []
                UINT_fit = []
                OE_fit = []
                KE_fit_rate = []
                UE_fit_rate = []
                UINT_fit_rate = []
                OE_fit_rate = []
                KE_Sigma = []
                UE_Sigma = []
                UINT_Sigma = []
                OE_Sigma = []
                dKE_fit = []
                dUE_fit = []
                dUINT_fit = []
                dOE_fit = []
                KE_rate = []
                UE_rate = []
                UINT_rate = []
                OE_rate = []
                KE_error = []
                UE_error = []
                UINT_error = []
                OE_error = []
                KE_Grad_Temp = []
                UE_Grad_Temp = []
                UINT_Grad_Temp = []
                OE_Grad_Temp = []
                dKE_Grad = []
                dUE_Grad = []
                dUINT_Grad = []
                dOE_Grad = []
                FE_Sum = []
                TOR_Sum = []
                ECC_Sum = []
                ECC_Rate = []
                ECC_Ang = []
                SEMI_Sum = []
                RAD_Sum = []
                THE_Sum = []
                VEC_Sum = []
                WE_Sum = []

                # ======================================
                # |
                # |
                # ======================================
                #   disk energy and Orbital Energy
                # ======================================
                # |
                # |
                # ======================================
                #
                #
                #

                Int = 0
                dInt = 1

                Orbit = np.round(timeCutOff)

                KE_Sum = []
                UE_Sum = []
                UINT_Sum = []

                while Int <= Orbit:
                    #
                    #
                    #
                    try:
                        if Calc_Temp == True:
                            ff = pc.read_var(
                                trimall=True, ivar=Int, magic=["TT"], quiet=True)
                            TT = ff.TT
                        else:
                            ff = pc.read_var(
                                trimall=True, ivar=Int, quiet=True)
                        rad = ff.x
                        phi = ff.y
                        uu = ff.uu
                        ux = ff.ux
                        uy = ff.uy
                        rho = ff.rho
                        i = 0
                        di = 1
                        j = 0
                        dj = 1
                        KE = []
                        UE = []
                        UINT = []
                        FE = []
                    except:
                        # moving on to next ivar
                        Int = Int + dInt

                    while j <= len(phi) - 1:
                        #
                        #
                        #
                        #
                        while i <= len(rad) - 1:
                            try:
                                local_disk_mass = 0.5 * np.abs(phi[j]) * ((rad[i + 1] ** 2) - (rad[i] ** 2)) * rho[i][j]
                                KE.append(
                                    local_disk_mass * (ux[j][i] ** 2 + (rad[i]**2)*uy[j][j] ** 2)
                                    * 0.5
                                )
                                UE.append(
                                    local_disk_mass
                                    / rad[i]
                                )
                                UINT.append(
                                    TT[i][j]
                                    / gamma
                                )
                                FE.append(
                                    rho[i][j]
                                    / rad[i] ** 2
                                )
                                i = i + di
                            except:
                                # out of bound
                                i = i + di
                        i = 0
                        j = j + dj
                    #
                    #
                    #
                    #
                    FE_Sum.append(np.sum(FE[:]))
                    KE_Sum.append(np.sum(KE[:]))
                    UE_Sum.append(np.sum(UE[:]))
                    UINT_Sum.append(np.sum(UINT[:]))
                    Int = Int + dInt
                i = 0
                di = 1
                Total_Disk_Energy = []

                while i <= len(KE_Sum) - 1:
                    Total_Disk_Energy.append(KE_Sum[i] + UINT_Sum[i] - UE_Sum[i])
                    i = i + di
                #
                #
                #
                #
                # get the Orbital energy if logical is matched
                #
                #
                #
                #
                # ======================================
                if Calc_OEnergy == True:
                    # ======================================
                    #
                    if ecc_int == 0.0:

                        # the case is simple

                        OE = -0.5 * par1 / ts.xq2

                        OE_Sum = []

                        i = 0
                        di = 1

                        while i <= len(Orbit_Len) - 2:
                            OE_Sum.append(
                                np.sum(
                                    OE[
                                    Orbit_Len[i]:Orbit_Len[i + 1]
                                    ]
                                )
                            )
                            i = i + di
                    #
                    #
                    #
                    else:
                        #
                        #	the case is more complex
                        #
                        # -----------------------------------------------
                        #
                        #
                        #
                        #
                        # ----------------------------------------------
                        i = 0
                        di = 1
                        while i <= len(Orbit_Len) - len(torqtotal) - 1:
                            #
                            #
                            #
                            #
                            TOR_Sum.append(
                                torqtotal[
                                    Orbit_Len[i]
                                ]
                            )
                            ECC_Sum.append(
                                ecc[
                                    Orbit_Len[i]
                                ]
                            )
                            ECC_Rate.append(
                                ecc_rate[
                                    Orbit_Len[i]
                                ]
                            )
                            ECC_Ang.append(
                                ecc_ang[
                                    Orbit_Len[i]
                                ]
                            )
                            SEMI_Sum.append(
                                semi_major[
                                    Orbit_Len[i]
                                ]
                            )
                            RAD_Sum.append(
                                radius[
                                    Orbit_Len[i]
                                ]
                            )
                            THE_Sum.append(
                                Omega[
                                    Orbit_Len[i]
                                ]
                            )
                            VEC_Sum.append(
                                np.sqrt(
                                    radius[Orbit_Len[i]] ** 2
                                    +
                                    Omega[Orbit_Len[i]] ** 2
                                )
                            )
                            i = i + di
                        #
                        #
                        #
                        # -----------------------------------------------
                        #
                        #
                        #
                        i = 0
                        i = i + di
                        while i <= len(FE_Sum) - 1:
                            try:
                                WE_Sum.append(
                                    VEC_Sum[i]
                                    *
                                    FE_Sum[i])
                                i = i + di
                            except:
                                #
                                # loop error
                                #
                                print('loop error')
                                i = i + di
                        #
                        #
                        #
                        # -----------------------------------------------
                        #
                        #
                        #
                        OE_Sum = []
                        i = 0
                        di = 1
                        while i <= Orbit - 2:
                            #
                            #
                            #
                            try:
                                #
                                #
                                #
                                OE_Sum.append(
                                    (2 / WE_Sum[i])
                                    *
                                    (
                                            (
                                                    (ECC_Rate[i] * ECC_Sum[i])
                                                    /
                                                    (1 - ECC_Sum[i] ** 2)
                                            )
                                            +
                                            (
                                                    (TOR_Sum[i])
                                                    /
                                                    (ECC_Ang[i])
                                            )

                                    )
                                )
                                #
                                #
                                #
                            except:
                                #
                                # index error
                                #
                                print('loop error')
                                i = i + di
                            #
                            #
                            #
                            i = i + di

                else:
                    print('===================================')
                    print('Orbital Energy not to be calculated')
                    print('===================================')
                # -----------------------------------------------
                # gradients
                # -----------------------------------------------
                KE_Grad = np.gradient(KE_Sum[:])
                UE_Grad = np.gradient(UE_Sum[:])
                UINT_Grad = np.gradient(UINT_Sum[:])
                if Calc_OEnergy == True:
                    OE_Grad = np.gradient(OE_Sum[:])

                KE_Sum_Avg = np.convolve(KE_Sum, kernel, mode='valid')

                KE_rate = np.mean(KE_Grad)
                UE_rate = np.mean(UE_Grad)
                UINT_rate = np.mean(UINT_Grad)
                if Calc_OEnergy == True:
                    OE_rate = np.mean(OE_Grad)

                # time arrays

                KE_time = np.arange(0, len(KE_Sum), 1)
                UE_time = np.arange(0, len(UE_Sum), 1)
                UINT_time = np.arange(0, len(UINT_Sum), 1)
                if Calc_OEnergy == True:
                    OE_time = np.arange(0, len(OE_Sum), 1)

                # attempt to make a fit by using the gradient of every energy array

                i = 0
                di = 1

                while i <= len(KE_time) - 1:
                    KE_fit.append((KE_rate)
                                  * KE_time[i] + KE_Sum[0])
                    i = i + di

                i = 0
                di = 1

                while i <= len(UE_time) - 1:
                    UE_fit.append(UE_rate
                                  * UE_time[i] + UE_Sum[0])
                    i = i + di

                i = 0
                di = 1

                while i <= len(UINT_time) - 1:
                    UINT_fit.append(UINT_rate
                                    * UINT_time[i] + UINT_Sum[0])
                    i = i + di

                i = 0
                di = 1

                if Calc_OEnergy == True:
                    while i <= len(OE_time) - 1:
                        OE_fit.append(
                            OE_Grad[i] * np.sin(OE_Sum[i] + OE_time[i]) + OE_Sum[0])
                        i = i + di

                #
                #
                #

                i = 0
                di = 1

                while i <= len(KE_fit) - 2:
                    dKE_fit.append((KE_fit)[i + 1] - (KE_fit)[i])
                    i = i + di
                i = 0
                di = 1

                while i <= len(UE_fit) - 2:
                    dUE_fit.append((UE_fit)[i + 1] - (UE_fit)[i])
                    i = i + di
                i = 0
                di = 1

                while i <= len(UINT_fit) - 2:
                    dUINT_fit.append((UINT_fit)[i + 1] - (UINT_fit)[i])
                    i = i + di
                i = 0
                di = 1

                if Calc_OEnergy == True:
                    while i <= len(OE_fit) - 2:
                        dOE_fit.append((OE_fit)[i + 1] - (OE_fit)[i])
                        i = i + di

                KE_fit_rate = sum(dKE_fit) / len(dKE_fit)
                UE_fit_rate = sum(dUE_fit) / len(dUE_fit)
                UINT_fit_rate = sum(dUINT_fit) / len(dUINT_fit)
                if Calc_OEnergy == True:
                    OE_fit_rate = sum(dOE_fit) / len(dOE_fit)

                KE_error = np.abs((KE_fit_rate - KE_rate) / KE_rate)
                UE_error = np.abs((UE_fit_rate - UE_rate) / UE_rate)
                UINT_error = np.abs((UINT_fit_rate - UINT_rate) / UINT_rate)
                if Calc_OEnergy == True:
                    OE_error = np.abs((OE_fit_rate - OE_rate) / OE_rate)

                # averages

                KE_Sum_Avg = np.convolve(KE_Sum, kernel, mode='valid')
                UE_Sum_Avg = np.convolve(UE_Sum, kernel, mode='valid')
                UINT_Sum_Avg = np.convolve(UINT_Sum, kernel, mode='valid')
                if Calc_OEnergy == True:
                    OE_Sum_Avg = np.convolve(OE_Sum, kernel, mode='valid')

                KE_Grad_Avg = np.convolve(KE_Grad, kernel, mode='valid')
                UE_Grad_Avg = np.convolve(UE_Grad, kernel, mode='valid')
                UINT_Grad_Avg = np.convolve(UINT_Grad, kernel, mode='valid')
                if Calc_OEnergy == True:
                    OE_Grad_Avg = np.convolve(OE_Grad, kernel, mode='valid')

            else:
                #
                # pass all vars because I dont want to calculate energy
                #

                KE_Sum = []
                UE_Sum = []
                UINT_Sum = []
                OE_Sum = []
                Total_Disk_Energy = []
                Total_Disk_Energy_Avg = []
                Total_Disk_Grad = []
                Total_Disk_Grad_Avg = []
                KE_Sum_Avg = []
                UE_Sum_Avg = []
                OE_Sum_Avg = []
                UINT_Sum_Avg = []
                KE_Grad = []
                UE_Grad = []
                OE_Grad = []
                UINT_Grad = []
                KE_Grad_Avg = []
                UE_Grad_Avg = []
                OE_Grad_Avg = []
                UINT_Grad_Avg = []
                KE_fit = []
                UE_fit = []
                UINT_fit = []
                OE_fit = []
                KE_fit_rate = []
                UE_fit_rate = []
                UINT_fit_rate = []
                OE_fit_rate = []
                KE_Sigma = []
                UE_Sigma = []
                UINT_Sigma = []
                OE_Sigma = []
                dKE_fit = []
                dUE_fit = []
                dUINT_fit = []
                dOE_fit = []
                KE_rate = []
                UE_rate = []
                UINT_rate = []
                OE_rate = []
                KE_error = []
                UE_error = []
                UINT_error = []
                OE_error = []
                KE_Grad_Temp = []
                UE_Grad_Temp = []
                UINT_Grad_Temp = []
                OE_Grad_Temp = []
                dKE_Grad = []
                dUE_Grad = []
                dUINT_Grad = []
                dOE_Grad = []
                FE_Sum = []
                TOR_Sum = []
                ECC_Sum = []
                ECC_Rate = []
                ECC_Ang = []
                SEMI_Sum = []
                RAD_Sum = []
                THE_Sum = []
                VEC_Sum = []
                WE_Sum = []

            # |
            # |
            # ======================================
            # |
            # |
            # ======================================

            Toomre = []
            rad_q = []
            theta_q = []
            rho_q = []

            if Calc_ToomreQ == True:
                #
                # takes ivar (orbit from dsnap)
                # make grid elements
                # uncomment below if the disk is not isothermal
                # ff=pc.read_var(trimall=True,ivar=ivar,magic=['TT'])
                ff_q = pc.read_var(trimall=True, ivar=Orbit, quiet=True)
                rad_q = ff_q.x  # disk grid points in r
                theta_q = ff_q.y  # disk grid points in theta
                rho_q = ff_q.rho  # disk surface density
                #
                pi = np.pi
                grav_const = 1
                #
                nx1, nx2 = 256, 768
                # cparam resolution. This needs to be set here

                Toomre = [[0 for x in range(nx1)] for y in range(nx2)]
                # redefine each element of TC_array with values for toomre Q
                j = 0
                i = 0
                di = 1
                dj = 1
                while j <= len(theta_q) - 1:
                    while i <= len(rad_q) - 1:
                        # append toomre Q value at that point
                        Toomre[j][i] = ((cs)) / (grav_const * pi * rho_q[j][i])
                        print('working at radius:' + str(i))
                        print('working at theta:' + str(j))
                        i = i + di
                    print('===========================')
                    print('moving to next theta')
                    print('===========================')
                    i = 0
                    j = j + dj
                print('===========================')
                print('done with toomre')
                print('===========================')
            # ======================================

            # Last line of all read_var processes

            # ======================================
            # ------------------------------------------------------

            # ======================================
            #
            #   TW2002,TW2004,CN2008
            #
            # ======================================

            tanaka = (-0.85 - alpha - 0.9 * beta) * Gamma0

            tanaka_array = np.repeat(tanaka, len(time))

            #
            # generate hermite polynomials
            # switch Z to value if run is 3D
            # hermite polynomials have a relation with radius and the z-axis
            # this term allows for the generation of 3D waves (bending waves)
            # for H(n) where n > 0
            #

            # radial_coord    = 0
            radial_coord = 1

            n = 1
            dn = 1

            size = 15 + n
            Herm_Array = []

            while n <= size:
                coef = np.repeat(1, n)
                Herm_Array.append(herm.hermval(0, coef))
                n = n + dn

            #   generate enthalpy fourier-hermite coefficents
            #   fourier-hermite coefficents found in TW2002 and TW2004
            #   relevant equations are
            #   TW2002: 1, 2, 3, 5, 6, 13, 14, 15, 17, 20, 22, 35
            #   TW2004:

            initial_pressure = Sigma * cs
            episilonCS = cs / semi_major * Omegap

            # ======================================
            #
            # test out radial coordinate generation
            #
            # ======================================

            m = 0
            dm = 1

            wavenumber_array = []

            while m <= size - 1:
                wavenumber_array.append(m * episilonCS)
                m = m + dm

            q = 0
            dq = 1

            wavenumber_check = []

            while q <= size - 1:
                wavenumber_check.append(2 * q / 3 * wavenumber_array[q])
                q = q + dq

            wavex = []
            wavey = []

            i = 0
            di = 1

            while i <= size - 1:
                wavex.append(x2d[5] + wavenumber_check[i][15])
                wavey.append(y2d[5] + wavenumber_check[i][15])
                i = i + di

            # ======================================

            I = 1
            pi = np.pi

            Num = 3
            lower_limit = 1e10 * -pi
            upper_limit = 1e10 * pi

            Z, w = self.Quadab(Num, lower_limit, upper_limit)

            # limit of I determined by number of coefs needed
            # while I <= 15
            Mu_coef = 0.0
            for k in range(Num):
                Mu_coef += w[k] * \
                           self.fourier_hermite_function((Z[k]), radial_coord, 1, I)

            # Mu_coef must converge

            Mu_coef = Mu_coef * 2 / upper_limit

            # get the second enthalpy fourier-hermite coef

            Num = 3
            lower_limit = 1e10 * -pi
            upper_limit = 1e10 * pi

            Z, w = self.Quadab(Num, lower_limit, upper_limit)

            Mu_coef2 = 0.0
            for k in range(Num):
                Mu_coef2 += w[k] * \
                            self.fourier_hermite_function(
                                (Z[k]), radial_coord, Mu_coef, I)

            Mu_coef2 = Mu_coef2 * 2 / upper_limit

            # get the third enthalpy fourier-hermite coef

            Num = 3
            lower_limit = 1e10 * -pi
            upper_limit = 1e10 * pi

            Z, w = self.Quadab(Num, lower_limit, upper_limit)

            Mu_coef3 = 0.0
            for k in range(Num):
                Mu_coef3 += w[k] * self.fourier_hermite_function(
                    (Z[k]), radial_coord, Mu_coef2, I)

            Mu_coef3 = Mu_coef3 * 2 / upper_limit

            # ======================================
            #
            #   last line should switch to previous directory
            #
            # ======================================

            vars_dict = {
                'DirName': self.Directory_Path, 'Standard_Orbit': Standard_Orbit,
                'ecc_int': ecc_int, 'ecc_rate': ecc_rate, 'ecc_ang': ecc_ang,
                'Dynamic': Dynamic, 'Dynamic_Density': Dynamic_Density, 'Dynamic_Shock': Dynamic_Shock,
                'Dynamic_Temperature': Dynamic_Temperature,
                'x_grid': x_grid, 'y_grid': y_grid,
                'Toomre': Toomre, 'rad_q': rad_q, 'theta_q': theta_q, 'rho_q': rho_q,
                'Total_Disk_Grad': Total_Disk_Grad, 'Total_Disk_Grad_Avg': Total_Disk_Grad_Avg,
                'KE_rate': KE_rate, 'UE_rate': UE_rate, 'UINT_rate': UINT_rate, 'OE_rate': OE_rate,
                'KE_error': KE_error, 'UE_error': UE_error, 'UINT_error': UINT_error, 'OE_error': OE_error,
                'KE_fit': KE_fit, 'UE_fit': UE_fit, 'UINT_fit': UINT_fit, 'OE_fit': OE_fit,
                'KE_fit_rate': KE_fit_rate, 'UE_fit_rate': UE_fit_rate, 'UINT_fit_rate': UINT_fit_rate,
                'OE_fit_rate': OE_fit_rate,
                'OE_Grad': OE_Grad, 'UE_Grad': UE_Grad, 'KE_Grad': KE_Grad, 'UINT_Grad': UINT_Grad,
                'OE_Grad_Avg': OE_Grad_Avg, 'UE_Grad_Avg': UE_Grad, 'KE_Grad_Avg': KE_Grad_Avg,
                'UINT_Grad_Avg': UINT_Grad_Avg,
                'OE_Sum': OE_Sum, 'UE_Sum': UE_Sum, 'KE_Sum': KE_Sum, 'UINT_Sum': UINT_Sum,
                'OE_Sum_Avg': OE_Sum_Avg, 'UE_Sum_Avg': UE_Sum_Avg, 'KE_Sum_Avg': KE_Sum_Avg,
                'UINT_Sum_Avg': UINT_Sum_Avg, 'Total_Disk_Energy': Total_Disk_Energy,
                'Total_Disk_Energy_Avg': Total_Disk_Energy_Avg,
                'CN_line': CN_line, 'TTm': TTm, 'TTm_rate': TTm_rate,
                'rad_grid': rad_grid, 'ivar': ivar, 'rho_fv': rho_fv, 'temp_fv': temp_fv, 'shock_fv': shock_fv,
                'avgrho_fv': avgrho_fv, 'avgtemp_fv': avgtemp_fv, 'avgshock_fv': avgshock_fv,
                'long_perihelion': long_perihelion, 'Init_Temp': Init_Temp, 'true_anomaly': true_anomaly,
                'true_angle': true_angle,
                'DGTemp_Mean': DGTemp_Mean, 'gamma': gamma, 'par1': par1, 'par2': par2,
                'GlobalTemp_Mean': GlobalTemp_Mean, 'GTM_Sigma': GTM_Sigma, 'GTM_Sigma_O': GTM_Sigma_O,
                'wavex': wavex, 'wavey': wavey, 'wavenumber_array': wavenumber_array,
                'wavenumber_check': wavenumber_check,
                'episilonCS': episilonCS, 'Mu_coef2': Mu_coef2, 'Mu_coef3': Mu_coef3,
                'ts': ts, 't': t, 'kernel': kernel, 'time': time,
                'torqint': torqint, 'torqext': torqext,
                'q': q, 'Mstar': Mstar, 'Gamma0': Gamma0, 'alpha': alpha, 'beta': beta, 'rsmooth': rsmooth,
                'gravC': gravC,
                'EntropyIndex': EntropyIndex, 'SpecificHeat': SpecificHeat,
                'Sigma': Sigma,
                'LinearVelocity': LinearVelocity, 'AngularVelocity': AngularVelocity,
                'eccentricity': eccentricity, 'semi_major': semi_major,
                'Hill_Radius': Hill_Radius, 'aspect_ratio': ScaleHeight,
                'Kepler_F': Kepler_F, 'Omegap': Omegap,
                'twave': twave, 'eh': eh, 'tanaka': tanaka, 'tanaka_array': tanaka_array, 'Herm_Array': Herm_Array[:],
                'cs': cs, 'initial_pressure': initial_pressure, 'Mu_coef': Mu_coef,
                'torqtotal': torqtotal, 'tmax': tmax, 'xrq2': xrq2, 'yrq2': yrq2, 'MaxOrbits': MaxOrbits,
                'tempdata': tempdata, 'Calc_Temp': self.Calc_Temp, 'step': step,
                'x2d': x2d, 'y2d': y2d, 'Calc_Density': self.Calc_Density, 'densitydata': densitydata,
                'timeCutOff': timeCutOff,'indexTimeCutOff':indexTimeCutOff, 'indexTimeCutOffLarge':indexTimeCutOffLarge,'timeCutOffLarge':timeCutOffLarge
            }

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

    def pingDATA(self, Directory_Number):

        self.Directory_Number = Directory_Number

        # ======================================
        #
        #   get torque var data from pencil datacubes
        #   requires directory number
        #
        # ======================================
        #
        #   DEV LOG
        #
        #   12/27/2018 added time series vars
        #
        # ======================================

        os.chdir('AGNRun' + str(Directory_Number))

        ts = pc.read_ts()
        t = ts.t / 2 / np.pi
        N = 50

        # ======================================
        # pc.read_pararm()
        # ======================================

        par = pc.read_param()
        h = par.cs0
        if (par.iprimary == 1):
            q = par.pmass[1]
        else:
            q = par.pmass[0]

        Gamma0 = (q / h) ** 2
        alpha = par.density_power_law
        beta = par.temperature_power_law

        # ======================================

        kernel = np.ones((N,)) / N
        time = np.convolve(t, kernel, mode='valid')
        torqint = np.convolve(ts.torqint_2, kernel, mode='valid')
        torqext = np.convolve(ts.torqext_2, kernel, mode='valid')

        # ======================================
        #
        #   last line should switch to previous directory
        #
        # ======================================

        print(os.getcwd())
        os.chdir('..')

        vars_dict = {'ts': ts,
                     't': t,
                     'kernel': kernel,
                     'time': time,
                     'torqint': torqint,
                     'torqext': torqext,
                     'q': q,
                     'Gamma0': Gamma0,
                     'alpha': alpha,
                     'beta': beta}
        return vars_dict

    def pingGit(self, Directory_Path):
        # send git commands
        self.Directory_Path = Directory_Path
        os.chdir(self.Directory_Path)
        try:
            print('=======================')
            print(' adding data to git')
            print('=======================')
            os.system('git add *.in')
            os.system('git add *.py')
            os.system('git add *.png')
            os.system('git add src/cparam.local')
            os.system('git add src/Makefile.local')
            os.chdir('..')

        except:
            print('git failed')
            os.chdir('..')

    def pingEccDecay(self, Directory_Path):

        self.Directory_Path = Directory_Path
        os.chdir(self.Directory_Path)

        try:
            os.system('ipython EccDecayCN2008_editWL.py')
            os.chdir('..')
        except:
            print('=======================')
            print('No EccDecayCN2008_editWL.py')
            print('=======================')
            traceback.print_exc()
            os.chdir('..')

    def pingTorque(self, Directory_Path):

        self.Directory_Path = Directory_Path
        os.chdir(self.Directory_Path)

        try:
            os.system('ipython Torque.py')
            os.chdir('..')
        except:
            print('=======================')
            print('No Torque.py')
            print('=======================')
            traceback.print_exc()
            os.chdir('..')

    def pingTemp(self, Directory_Path):

        self.Directory_Path = Directory_Path
        os.chdir(self.Directory_Path)

        try:
            os.system('ipython PingTemp.py')
            os.chdir('..')
        except:
            print('=======================')
            print('No PingTemp.py')
            print('=======================')
            traceback.print_exc()
            os.chdir('..')

    # ======================================
    #
    #   integrator functions used for fourier-hermite polynomial generation
    #
    # ======================================

    def Quadw(self, Num):

        a = np.linspace(3, 4 * Num - 1, Num) / (4 * Num + 2)
        x = np.cos(np.pi * a + 1 / (8 * Num * Num * np.tan(a)))
        epsilon = 1e-15
        delta = 1.0
        while delta > epsilon:
            p0 = np.ones(Num, float)
            p1 = np.copy(x)
            for k in range(1, Num):
                p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
            dp = (Num + 1) * (p0 - x * p1) / (1 - x * x)
            dx = p1 / dp
            x -= dx
            delta = max(abs(dx))

        w = 2 * (Num + 1) * (Num + 1) / (Num * Num * (1 - x * x) * dp * dp)

        return x, w

    def Quadab(self, Num, a, b):
        Z, w = self.Quadw(Num)
        return 0.5 * (b - a) * Z + 0.5 * (b + a), 0.5 * (b - a) * w

    def fourier_hermite_function(self, Z, radial_coord, Mu_coef, i):
        if radial_coord == 0:
            return 1
        else:
            return np.exp((-Z ** 2.0) / 2.0) * (1.0 / 2.0 * np.pi * i) * Mu_coef

    def make_grid_toomre(self, Orbit):
        # takes ivar (orbit from dsnap)
        # make grid elements
        # uncomment below if the disk is not isothermal
        # ff=pc.read_var(trimall=True,ivar=ivar,magic=['TT'])
        ff_q = pc.read_var(trimall=True, ivar=Orbit)
        rad_q = ff_q.x  # disk grid points in r
        theta_q = ff_q.y  # disk grid points in theta
        ux_q = ff_q.ux  # disk grid points in vr
        uy_q = ff_q.uy  # disk grid points in vtheta
        rho_q = ff_q.rho  # disk surface density
        rad2d_q, theta2d_q = np.meshgrid(rad_q, theta_q)
        x2d_q = rad2d_q * np.cos(theta2d_q)
        y2d_q = rad2d_q * np.sin(theta2d_q)
        return rad_q, theta_q, rho_q
