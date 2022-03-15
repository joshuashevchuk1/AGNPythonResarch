#
# test out dir making script
#

import sys
import os
import traceback

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

i = 0
di = 1

ecc_start = 0.1
ecc_step = 0.1
ecc_fixed = int(ecc_start)

while i <= len(dir_run_list)-1:
    os.chdir(str(dir_run_list[i]))
    try:
        print(i)
        replace_in_file('start.in', 'pmass', '1.0,' +
                        str(ecc_start-ecc_fixed+ecc_step)+'e-5')
    except:
        print('nothing to replace')
        traceback.print_exc()
        print(i)
    os.chdir('..')
    ecc_start = ecc_start+ecc_step
    i = i+di
