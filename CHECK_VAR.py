import numpy as np
import pandas as pd
import re
import os
import system

# check to see if var data is missing

try:
    os.chdir('data/proc0')
except:
    print('no proc0')

var = open('varN.list', 'r')
vline = var.readlines()

temp = []
store = []
missing = []

for line in vline:
    temp.append(line)

for index in range(len(temp)):
    temp_string = re.search('(?<=VAR)[0-9]{0,100}', str(temp[index]))
    store.append(temp_string.group(0))

i = 0
di = 1
count = 0

for index in range(len(store)):
    if i == int(store[index]):
        i = i+di
    else:
        while i <= int(store[index]):
            if i == int(store[index]):
                i = i+di
            else:
                missing.append('missing VAR'+str(i))
                count = count + 1
                i = i+di
    if i == len(store):
        missing.append('end with :'+str(count)+' missing')

os.chdir('..')
print(missing)
