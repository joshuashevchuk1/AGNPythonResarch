import numpy as np
import matplotlib.pyplot as plt
import pencil as pc
import math

# quick set up for intervals of time
# use Orbit_Len to properly get the length of the interval
# for Orbits from the time series

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
# ------------------
#
#   Debug script to test energy iteration quickly using python
#   Not part of standard analysis (Pencil_Analysis.py)
#
# ------------------


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
Omega = ts.yq2
LinearVelocity = ts.vxq2
AngularVelocity = ts.vyq2
TTm = ts.TTm
GlobalTemp_Mean = ts.TTm
AngularMomentum = AngularVelocity*par1

v2 = LinearVelocity**2 + AngularVelocity**2
semi_major = 1./(2/radius-v2)
DArclength = radius**2*(AngularVelocity/radius)
ep1 = (DArclength**2)/semi_major
eccentricity = (1-ep1)**0.5
Omegap = 1./semi_major**1.5
Kepler_F = np.sqrt(gravC*Mstar/radius)
aspect_ratio = cs/Kepler_F
Hill_Radius = semi_major*((q/3.0)**(1.0/3.0))
xrq2 = radius*np.cos(Omega)
yrq2 = radius*np.sin(Omega)
MaxOrbits = int(round(tmax))

# ======================================
#   disk energy and Orbital Energy
# ======================================

# get the disk and Orbital energy

OE = -0.5*par1/ts.xq2

OE_Sum = []

i = 0
di = 1

while i <= len(Orbit_Len)-2:
    OE_Sum.append(
        np.sum(
            OE[
                Orbit_Len[i]:Orbit_Len[i+1]
            ]
        )
    )
    i = i+di

print("=====================")
print("DEBUG")
print("=====================")
print('len information')
print(len(OE))
print(len(OE_Sum))
print(OE_Sum)
print("=====================")


plt.figure(figsize=(10, 10))
plt.title('Energy Magnitude Comparison')
plt.plot(np.abs(OE_Sum[:]), color='orange', label='Orbital', ls=':')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('Orbital_fix.png')
plt.close()
