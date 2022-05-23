import numpy as np
import os
import PencilContoursData as PD
import PencilContoursAnalysis as PA
import matplotlib.pyplot as plt
import configparser
import logging
import traceback
import sys
import re

config = configparser.ConfigParser()
configfile = config.read('Pencil.cfg')

Rawconfig = configparser.RawConfigParser()
Rawconfigfile = Rawconfig.read('Pencil.cfg')

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

PA = PA.Pencil_Analysis(PD.Pencil_Data(),
                        Orbit=0,
                        Calc_Temp=True,
                        Calc_Density=True)

var_dir_list = PA.Make_Vars()

PA.pingContour(var_dir_list)

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