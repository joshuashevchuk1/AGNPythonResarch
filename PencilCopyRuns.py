#
# test out dir making script
#

import sys
import os

run_start = 927
loop_length = 8

i = 0
di = 1

while i <= loop_length:
    os.system('cp -rf BASE/ AGNRun'+str(run_start+i))
    i = i+di
