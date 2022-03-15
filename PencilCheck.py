#
# test out dir making script
#

import sys
import os

run_start = 900
loop_length = 20

i = 0
di = 1

while i <= loop_length:
    os.system('cp -rf Base/ AGNRun'+str(run_start+i))
    i = i+di
