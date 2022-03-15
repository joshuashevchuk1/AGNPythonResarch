import pencil_old as pc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

ts = pc.read_ts()
t = ts.t
time = t/2*math.pi
len_t = len(ts.t)

Max_Orbits = np.round((time[len(time)-1]))

n = 1
dn = 1

i = 0
di = 1

Orbit_Len = []

while i <= len(time)-1:
    if np.round(n*2.0*math.pi) == np.round(t[i]):
        Orbit_Len.append(len(t[:i]))
        n = n+dn
        i = i+di
    else:
        i = i+di

print('===================')
print(len(Orbit_Len))
print('===================')

# ======================================
# pc.read_ts()
# ======================================

t = ts.t/2/np.pi
N = 50
tmax = t.max()
ux = ts.ux2m
uy = ts.uy2m
rhomax = ts.rhomax
rhomin = ts.rhomin

# quick set up for intervals of time
# use Orbit_Len to properly get the length of the interval
# for Orbits from the time series

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
Gamma0 = (q/h)**2
alpha = par.density_power_law
beta = par.temperature_power_law
kernel = np.ones((N,))/N
time = np.convolve(t, kernel, mode='valid')
torqint = np.convolve(ts.torqint_2, kernel, mode='valid')
torqext = np.convolve(ts.torqext_2, kernel, mode='valid')
torqtotal = torqext[:]+torqint[:]
rsmooth = par.r_smooth[1]
gravC = par.g0
EntropyIndex = beta-((gamma-1)*alpha)
SpecificHeat = par.cp
Sigma = par.rho0
Mstar = par.pmass[0]
cs = par.cs0  # sound speed
bi = par.xyz0  # disk boundary radi (used in area calculations)
bf = par.xyz1
semi_major_int = par.semimajor_axis
ecc_int = par.eccentricity
mu_mass = (par1+par2)/(par1*par2)  # reduced mass of primary and secondary

# ======================================

radius = ts.xq2
radiusS = ts.xq1
Omega = ts.yq2
LinearVelocity = ts.vxq2
AngularVelocity = ts.vyq2
LinearVelocityS = ts.vxq1
AngularVelocityS = ts.vyq1
TTm = ts.TTm
GlobalTemp_Mean = ts.TTm
AngularMomentum = AngularVelocity*par1

v2 = LinearVelocity**2 + AngularVelocity**2
v2S = LinearVelocityS**2+AngularVelocityS**2
semi_major = 1./(2/radius-v2)
semi_majorS = 1./(2/radiusS-v2S)
DArclength = radius**2*(AngularVelocity/radius)
ep1 = (DArclength**2)/semi_major
eccentricity = (1-ep1)**0.5
Omegap = 1./semi_major**1.5

eccentricity = np.convolve(eccentricity, kernel, mode='valid')

ecc = []
for index in range(len(Orbit_Len)):
    ecc.append(eccentricity[Orbit_Len[index]])

data=pd.read_csv("UINT_Sum.csv")

UINT_Sum_df = pd.DataFrame(data)
UINT=UINT_Sum_df.iloc[:,1].to_list()
t=np.arange(0,len(UINT),1)

UINT_integrate=(UINT[len(UINT)-1])
UINT_percent_increase = UINT_integrate/UINT[0]
data={"UINT Sum":UINT_integrate,"UINT percent increase":UINT_percent_increase,"UINT[0]":UINT[0]}
UINT_integrate_df = pd.DataFrame(data=data,index=[0])
UINT_integrate_df.to_csv('UINT_integrate_df.csv')


#print(len(UINT))
print(UINT_integrate)
print(UINT_percent_increase)
print(UINT[0])

plt.plot(UINT,label='AGNRun1884',ls=":")
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgUINT_1884_paper.png')
