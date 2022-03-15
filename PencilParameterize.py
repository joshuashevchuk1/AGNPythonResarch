#
# test out dir making script
#

import sys
import os

print('placeholder')


def replace_in_file(filename, key, new_value):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        if line.split('=')[0].strip(' \n') == key:
            lines[i] = key + ' = ' + new_value + '\n'
    f = open(filename, "w")
    f.write("".join(lines))
    f.close()


dir_run_list = next(os.walk('.'))[1]

print(dir_run_list)

os.chdir(str(dir_run_list[1]))

# 3 is a placeholder


length = 20
run_start = 580
run_step = 1
run_length = int(run_start)

while run_start <= length+run_length:
    os.system('pc_newrun AGNRun'+str(run_start-1+run_step))
    print(run_start)
    run_start = run_start+run_step

os.chdir('..')

dir_run_list = next(os.walk('.'))[1]

i = 0
di = 1

ecc_start = 0.1
ecc_step = 0.1
ecc_fixed = int(ecc_start)

while i <= len(dir_run_list)-1:
    os.chdir(str(dir_run_list[i]))
    if ecc_start == 0.5:
        ecc_start = 0.1
    try:
        replace_in_file("start.in", 'eccentricity',
                        str(ecc_start-ecc_fixed+ecc_step))
    except:
        print('nothing to replace')
    os.system('pc_build -H login1.stampede2.tacc.utexas.edu.conf')
    os.chdir('..')
    ecc_start = ecc_start+ecc_step
    i = i+di
