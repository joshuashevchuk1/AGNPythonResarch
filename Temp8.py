import pencil_old as pc
import matplotlib.pyplot as plt
import numpy as np
import heapq
import pandas as pd

# 3/3/2021
# use to investigate relative temperature at an orbit

ts = pc.read_ts()
f0 = pc.read_var(trimall=True, ivar = 0 , magic = ["TT"]) 
#f25 = pc.read_var(trimall=True, ivar = 25, magic=["TT"])
f50 = pc.read_var(trimall=True, ivar = 8, magic=["TT"])
#f75 = pc.read_var(trimall=True, ivar = 75, magic=["TT"])

T0 = f0.TT
#T25 = f25.TT
T50 = f50.TT
#T75 = f50.TT

rho50 = f50.rho
rho50_index = np.argsort(rho50)

T0 = np.mean(np.transpose(T0),axis=1)
#T25 = np.mean(np.transpose(T25),axis=1)
T50 = np.mean(np.transpose(T50),axis=1)
#T75 = np.mean(np.transpose(T75),axis=1)

#T25 = T25-T0
T50 = T50-T0
#T75 = T75-T0

rad = f0.x

radq2 = ts.xq2
thetaq2 = ts.yq2

xrq2 = radq2*np.cos(thetaq2)
yrq2 = radq2*np.sin(thetaq2)

radmax = np.max(xrq2)
radmin = np.min(xrq2)

maxindex = 0

T50max=np.max(T50)
T50min=np.min(T50)

index_min = min(range(len(T50)), key=T50.__getitem__)
index_max = max(range(len(T50)), key=T50.__getitem__)

radmax_index_planet=min(range(len(rad)), key=lambda i: abs(rad[i]-radmax))
radmin_index_planet=min(range(len(rad)), key=lambda i: abs(rad[i]-np.abs(radmin)))

T50_planet = T50[radmin_index_planet:radmax_index_planet]

print('len T50')
print(len(T50))
print('len T50_planet')
print(len(T50_planet))
print('T50_planet')
print(T50_planet)
print('T50_planet max and min')
print(np.max(T50_planet))
print(np.min(T50_planet))

index_min_planet = min(range(len(T50_planet)), key=T50_planet.__getitem__)
index_max_planet = max(range(len(T50_planet)), key=T50_planet.__getitem__)

print('radmax_index_planet')
print(radmax_index_planet)
print('radmin_index_planet')
print(radmin_index_planet)
print('index max planet')
print(index_max_planet+radmax_index_planet)
print('index min planet')
print(index_min_planet+radmin_index_planet)

print('rad[radmax_index_planet]')
print(rad[radmax_index_planet])
print('rad[radmin_index_planet]')
print(rad[radmin_index_planet])

position_dict = {'min_planet':rad[index_min_planet+radmin_index_planet],
                 'max_planet':rad[index_max_planet+radmin_index_planet],
                 'radmax':radmax,
                 'radmin':radmin,
                 'max':rad[index_max],
                 'min':rad[index_min]}

position_dataframe = pd.DataFrame(data=position_dict,index=[0])
position_dataframe.to_csv('position_temp.csv')

ygrid = np.linspace(T50min,T50max,511)

Radmax=[]
Radmin=[]
Index_min=[]
Index_max=[]
Index_min_planet=[]
Index_max_planet=[]

for i in range(len(ygrid)):
    Radmax.append(np.abs(radmax))
    Radmin.append(np.abs(radmin))
    Index_min.append(np.abs(rad[index_min]))
    Index_max.append(np.abs(rad[index_max]))
    Index_min_planet.append(np.abs(rad[index_min_planet+radmin_index_planet]))
    Index_max_planet.append(np.abs(rad[index_max_planet+radmin_index_planet]))


print(radmin)
print(rad[index_max])

plt.title('AGNRun1884 (q=2e-5)')
#plt.plot(rad,T0)
#plt.plot(rad,T25,label = 't=25')
plt.plot(rad,T50,label = 't=8')
#plt.plot(rad,T75,label = '75')
plt.plot(Radmax,ygrid,label='Orbit max',ls=':',color='red')
plt.plot(Radmin,ygrid,label='Orbit min',ls=':',color='blue')
#plt.plot(Index_max,ygrid,label='max temperature',ls='-.',color='red')
#plt.plot(Index_min,ygrid,label='min temperature',ls='-.',color='orange')
plt.plot(Index_max_planet,ygrid,label='max temperature planet',ls='--',color='green')
plt.plot(Index_min_planet,ygrid,label='min temperature',ls='--',color='purple')
plt.scatter(rad(rho50_index[1]),T50(rho50_index[1]))
#plt.plot(2,label='radmin')
plt.xlim([0.20999,2.499])
plt.xlabel('radius')
plt.ylabel(r'$\Delta T$')
#plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.grid(True)
#plt.savefig('TempRelative_1884.png')
plt.show()
