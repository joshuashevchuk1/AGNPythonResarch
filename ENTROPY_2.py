import numpy as np 
import matplotlib.pyplot as plt
import pencil as pc
import math

global S_Array

#quick set up for intervals of time
#use Orbit_Len to properly get the length of the interval 
#for Orbits from the time series

ts=pc.read_ts()
t=ts.t
time=t/2*math.pi
len_t = len(ts.t)
N=50

Max_Orbits=np.round((time[len(time)-1]))

n=1
dn=1

i=0
di=1

Orbit_Len=[]

while i <= len(time)-1:
	if np.round(n*2.0*math.pi) == np.round(t[i]):
                Orbit_Len.append(len(t[:i]))
                n=n+dn
		i=i+di
        else:
               	i=i+di

#
#------------------
#
#   Debug script to test energy iteration quickly using python
#   Not part of standard analysis (Pencil_Analysis.py)
#
#------------------
#
#   
#======================================
#pc.read_pararm()
#======================================

par=pc.read_param()
h=par.cs0
if (par.iprimary == 1):
    q=par.pmass[1]
else:
    q=par.pmass[0]

par1            = par.pmass[1]
par2            = par.pmass[0]
gamma           = par.gamma
Gamma0          = (q/h)**2
alpha           = par.density_power_law
beta            = par.temperature_power_law
kernel          = np.ones((N,))/N
time            = np.convolve(t,kernel,mode='valid')
torqint         = np.convolve(ts.torqint_2,kernel,mode='valid')
torqext         = np.convolve(ts.torqext_2,kernel,mode='valid')
torqtotal       = torqext[:]+torqint[:]
rsmooth         = par.r_smooth[1]
gravC           = par.g0
EntropyIndex    = beta-((gamma-1)*alpha)
SpecificHeat    = par.cp
Sigma           = par.rho0
Mstar           = par.pmass[0]
cs              = par.cs0               #sound speed
bi              = par.xyz0	    #disk boundary radi (used in area calculations)
bf              = par.xyz1
semi_major_int  = par.semimajor_axis
ecc_int	    	= par.eccentricity
mu_mass	    	= (par1+par2)/(par1*par2) # reduced mass of primary and secondary

#======================================

radius          = ts.xq2
Omega           = ts.yq2
LinearVelocity  = ts.vxq2
AngularVelocity = ts.vyq2
TTm		= ts.TTm
GlobalTemp_Mean = ts.TTm
AngularMomentum = AngularVelocity*par1

v2              = LinearVelocity**2 + AngularVelocity**2
semi_major      = 1./(2/radius-v2)
DArclength      = radius**2*(AngularVelocity/radius)
ep1             = (DArclength**2)/semi_major
eccentricity    = (1-ep1)**0.5
Omegap          = 1./semi_major**1.5
Kepler_F        = np.sqrt(gravC*Mstar/radius)
aspect_ratio    = cs/Kepler_F
Hill_Radius     = semi_major*((q/3.0)**(1.0/3.0))
xrq2            = radius*np.cos(Omega)
yrq2            = radius*np.sin(Omega)
#MaxOrbits       = int(round(tmax))
ecc		= eccentricity
ecc_rate	= np.gradient(eccentricity)
ecc_ang		= par1*np.sqrt(semi_major)*np.sqrt(1-ecc**2)

if ecc_int == 0.0:
	print('is zero')
else:
	print('is not zero')
#
#-------------------------------------------
#
# CONTROl ORBITAL INFORMATION HERE!
#
#-------------------------------------------
#


Int=0
dInt=15

#Orbit=MaxOrbits-2
Orbit=800

#
#
#

f0=pc.read_var(trimall=True,ivar=0,magic=["TT"])
s0=f0.ss

S_Array =   []

while Int <= Orbit:

    #ff	=	pc.read_var(trimall=True,ivar=Int,magic=["TT"])
    ff	=	pc.read_var(trimall=True,ivar=Int,magic=["TT"])
    i	=	0
    di	=	1
    j	=	0
    dj	=	1
    ss  =       ff.ss
    rad =       ff.x
    phi =       ff.y
    #S   =       []

    #======
    #
    #   loop through vars and get azimuthal average for entropy. 
    #
    #======

    S_Array.append(np.mean(ss-s0,axis=0))
        
    Int=Int+dInt

i=0
di=1

apoapsis    = []
periapsis   = []

while i <= len(Orbit_Len)-1:
        apoapsis.append(
                semi_major[i]
                *
                (1+ecc[i])
                     )
        periapsis.append(
                semi_major[i]
                *
                (1-ecc[i])
                     )
        i=i+di



time = np.arange(0,len(S_Array)*dInt,1*dInt)

# redine apopasis points and periapsis points to have 
# a limit in number of orbits based on the length of 
# the time array

i=0
di=1

temp_a=[]
temp_p=[]

while i <= len(time)-1:
    temp_a.append(apoapsis[i])
    temp_p.append(periapsis[i])
    i=i+di

apoapsis=temp_a
periapsis=temp_p


#plt.contourf(rad,time,S_Array,256)
plt.plot(apoapsis,time,color='black',ls='--',label='apoapsis')
plt.plot(periapsis,time,color='red',ls=':',label='periapsis')
plt.title('Azimuthal average Entropy vs time and radius ')
plt.ylabel('time')
plt.xlabel('radius')
#plt.colorbar()
plt.legend()
plt.show()
