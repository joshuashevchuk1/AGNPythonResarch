import numpy as np
import matplotlib.pyplot as plt
import pencil_old as pc
import math
import os
import sys
import traceback

cos = np.cos
sin = np.sin

run_number = 1850

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

print('===================')
print(len(Orbit_Len))
print('===================')

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
ecc = eccentricity
ecc_rate = np.gradient(eccentricity)
ecc_ang = par1*np.sqrt(semi_major)*np.sqrt(1-ecc**2)

if ecc_int == 0.0:
    print('is zero')
else:
    print('is not zero')
Int = 0
dInt = 1

#
# CONTROl ORBITAL INFORMATION HERE!
#

# Orbit=MaxOrbits-2
Orbit = 6

#
#
#

KE_Sum = []
UE_Sum = []
UINT_Sum = []
FE_Sum = []

while Int <= Orbit:

    ff = pc.read_var(trimall=True, ivar=Int, magic=["TT"])
    TT = ff.TT
    rad = ff.x
    phi = ff.y
    uu = ff.uu
    ux = ff.ux
    uy = ff.uy
    rho = ff.rho
    drad = np.gradient(rad)
    dphi = np.gradient(phi)
    radp = []
    for index in range(len(Orbit_Len)):
        radp = radius[Orbit_Len]
    phip = []
    for index in range(len(Orbit_Len)):
        phip = Omega[Orbit_Len]
    dA = []
    for i in range(len(drad)):
        dA.append(rad[i] * drad[i] * np.mean(dphi))
    i = 0
    di = 1
    j = 0
    dj = 1
    KE = []
    UE = []
    UINT = []
    FE = []
    UINT_Specific = []

    while j <= len(phi)-1:
        while i <= len(rad)-1:
            try:
                Rad_Cell = np.sqrt(
                    ((rad[i]*cos(phi[j]))**2) +
                    ((rad[i]*sin(phi[j]))**2)
                )
                RadP_Cell = []
                try:
                    RadP_Cell = np.sqrt(
                        ((radp[Int]*cos(phip[Int]))**2) +
                        ((radp[Int]*cos(phip[Int]))**2)
                    )
                except:
                    RadP_Cell = np.sqrt(
                        ((radp[Int-1]*cos(phip[Int-1]))**2) +
                        ((radp[Int-1]*cos(phip[Int-1]))**2)
                    )
                KE.append(0.5*rho[j][i]*(ux[j][i]**2+uy[j][i]**2)*dA[i])
                UE.append(
                    rho[j][i]*dA[i] *
                    (
                        ((1-q)/Rad_Cell) +
                        (q/(RadP_Cell-Rad_Cell))
                    )
                )
                UINT.append(
                    dA[i] *
                    ((rho[j][i]
                      * TT[j][i])
                     / gamma)
                )
                UINT_Specific.append(
                    dA[i] *
                    (TT[j][i]/gamma))
                i = i+di
            except:
                traceback.print_exc()
                # out of bound

                i = i+di
        i = 0
        j = j+dj

    KE_Sum.append(np.sum(KE[:]))
    UE_Sum.append(np.sum(UE[:]))
    UINT_Sum.append(np.sum(UINT[:]))
    UINT_Specific.append(np.sum(UINT[:]))
    Int = Int+dInt

i = 0
di = 1
Total = []

while i <= len(KE_Sum)-1:
    Total.append(KE_Sum[i]+UINT_Sum[i]-UE_Sum[i])
    i = i+di

print("-----------------------------")
print("=====================")
print("DEBUG")
print("=====================")

# print(len(KE))
# print(len(rad))
# print(len(KE_Sum))
print(KE_Sum)
# print(UE_Sum)
# print(UINT_Sum)
print(Total)
# print('len information')
# print(len(OE))
# print(len(Orbit_Len))
# print(len(OE_Sum))
if ecc_int == 0.0:
    print('is zero')
else:
    print('is not zero')
# print(OE_Sum)
print("=====================")
print("-----------------------------")

plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
#plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
#plt.plot(UE_Sum, color='red', label='Potential')
#plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgTotal'+str(run_number)+'.png')
plt.close()

plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
plt.plot(UE_Sum, color='red', label='Potential')
plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgMC'+str(run_number)+'.png')
plt.close()

plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
#plt.plot(UE_Sum, color='red', label='Potential')
#plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
#plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgKE'+str(run_number)+'.png')
plt.close()

plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
#plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
plt.plot(UE_Sum, color='red', label='Potential')
#plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
#plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgUE'+str(run_number)+'.png')
plt.close()

plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
#plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
#plt.plot(UE_Sum, color='red', label='Potential')
plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
#plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('ErgUINT'+str(run_number)+'.png')
plt.close()
