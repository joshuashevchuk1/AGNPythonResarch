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

global t


class Pencil_Data(object):

    def __init__(self):
        self.sma = None
        self.ecc_int = None
        self.ecc = None
        self.gravC = None
        self.torqtotal = None
        self.torqext = None
        self.torqint = None
        self.kernel = None
        self.beta = None
        self.alpha = None
        self.Gamma0 = None
        self.gamma = None
        self.par2 = None
        self.par1 = None
        self.q = None
        self.h = None
        self.tmax = None
        self.time = None
        self.N = None
        self.t = None

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

        print(os.getcwd())
        os.chdir(self.Directory_Path)
        print(os.getcwd())

        try:
            self.initVars()
            vars_dict = {'t':self.t,'N':self.N,'sma':self.sma,'ecc_int':self.ecc_int,'ecc':self.ecc,
                         'gravC':self.gravC}
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

    def initVars(self):

        # ======================================
        # calculated values
        # ======================================

        ts = pc.read_ts()
        t = ts.t / 2 * math.pi

        N = 850
        tmax = t.max()

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
        rsmooth = par.r_smooth[1]
        radius = ts.xq2
        LinearVelocity = ts.vxq2
        AngularVelocity = ts.vyq2

        v2 = LinearVelocity ** 2 + AngularVelocity ** 2
        semi_major = 1. / (2 / radius - v2)
        DArclength = radius ** 2 * (AngularVelocity / radius)
        ep1 = (DArclength ** 2) / semi_major
        eccentricity = (1 - ep1) ** 0.5

        ecc = eccentricity
        ecc_int = par.eccentricity
        sma = semi_major

        # ======================================

        # ======================================
        # set calculated values
        # ======================================

        self.t = t
        self.N = N
        self.time = time
        self.tmax = tmax
        self.h = h
        self.q = q
        self.par1 = par1
        self.par2 = par2
        self.gamma = gamma
        self.Gamma0 = Gamma0
        self.alpha = alpha
        self.beta = beta
        self.kernel = kernel
        self.gravC = gravC
        self.ecc = ecc
        self.ecc_int = ecc_int
        self.sma = sma
