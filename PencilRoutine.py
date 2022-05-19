# =====================================================================
# ---------------------------------------------------------------------
#
#    analysis in the cwd. Go through Each directory and search
#   for data. Complete Standard Pencil-Code Analysis Routine.
#   Send all data to github for inspection.
#
# ---------------------------------------------------------------------
# =====================================================================
# ---------------------------------------------------------------------

import numpy as np
import os
import PencilData as PD
import PencilAnalysis as PA
import PencilPostAnalysis as PAP
import matplotlib.pyplot as plt
import configparser
import logging
import traceback
import sys
import re

# ---------------------------------------------------------------------
# =====================================================================
# ---------------------------------------------------------------------
# Get Log Config
# ---------------------------------------------------------------------

config = configparser.ConfigParser()
configfile = config.read('Pencil.cfg')

Rawconfig = configparser.RawConfigParser()
Rawconfigfile = Rawconfig.read('Pencil.cfg')

# print(Rawconfig.sections())
PencilLog = Rawconfig.sections()

if 'LOGGING' in PencilLog:
    print('=====================================================')
    print('found logging')
    print('=====================================================')
    try:
        filename = Rawconfig.items('LOGGING')

        for loggingvalue in filename[0]:
            loggingvalue = re.sub('\.logfiles$', '', loggingvalue)

        logging.basicConfig(filename=loggingvalue,
                            level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                            datefmt=' %m/%d/%Y %I:%M:%S %p ')
        logging.info('the program has started')
        print('=====================================================')
        print('making log file')
        print('=====================================================')
    except:
        print('=====================================================')
        print('logging file failed to make')
        print('=====================================================')
        traceback.print_exc()
else:
    print('=====================================================')
    print('no LOGGING section found')
    print('making default logging')
    print('=====================================================')
    try:
        logging.basicConfig(filename='logging.log',
                            level=logging.DEBUG, format=' %(asctime)s,%(message)s ',
                            datefmt=' %m/%d/%Y %I:%M:%S %p ')
        logging.info('the program has started')
        print('=====================================================')
        print('making log file')
        print('=====================================================')
    except:
        print('=====================================================')
        print('logging file failed to make')
        print('=====================================================')
        traceback.print_exc()

#PA                                              =   PA.Pencil_Analysis(PD.Pencil_Data(),Orbit=None,step=10,TempSigma=True)
PA = PA.Pencil_Analysis(PD.Pencil_Data(),
                        Orbit=5,
                        MaxOrbits=5,
                        step=None,
                        TempSigma=False,
                        Calc_Energy=False,
                        Calc_OEnergy=False,
                        Calc_Dynamics=False,
                        Orbit_standard=None,
                        Calc_Rates_Energy=False,
                        Calc_Temp=False,
                        Calc_ToomreQ=False)

#PAP                                             =   PAP.Pencil_Post_Analysis(PD.Pencil_Data(),Orbit=30,step=10,TempSigma=True)

var_dir_list = PA.Make_Vars()

#temp_dir_list, ecc_dir_list, name_dir_list      =   PA.Make_TempVars(var_dir_list)
#temp_distribution , temp_ecc, temp_name         =   PA.Get_TempDistribution(temp_dir_list,ecc_dir_list,name_dir_list)

#PA.pingContour(var_dir_list)
#PA.pingContour_Dynamic(var_dir_list)
#PA.pingContour_Midplane(var_dir_list)
#PA.pingTorque(var_dir_list)
#PA.pingEcc(var_dir_list)
PA.pingOrbital(var_dir_list)
#PA.pingTTm(var_dir_list)
#PA.pingTempBoxPlot(var_dir_list)
#PA.pingGTM_Sigma(var_dir_list)
#PA.pingGTM_Sigma_O(var_dir_list)
#PA.pingGTM_O_DT(var_dir_list)
#PA.pingGTM_Sigma_3D(var_dir_list)
#PA.pingGTM_Sigma_3D_O(var_dir_list)
#PA.ping_LongPerihelion(var_dir_list)
#PA.ping_LongPerihelion_10_Orbits(var_dir_list)
#PA.pingEnergy(var_dir_list)
#PA.pingToomreQ(var_dir_list)

# PAP.pingContour(var_dir_list)
# PAP.pingTorque(var_dir_list)
# PAP.pingEcc(var_dir_list)
# PAP.pingOrbital(var_dir_list)
# PAP.pingTempBoxPlot(var_dir_list)
# PAP.ping_LongPerihelion(var_dir_list)
# PAP.ping_LongPerihelion_10_Orbits(var_dir_list)

i = 0
di = 1

varmaxtimes = []

while i <= len(var_dir_list)-1:
    try:
        varmaxtimes.append(max(var_dir_list[i]['t']))
    except:
        print('no time series')
    i = i+di

print(varmaxtimes)
