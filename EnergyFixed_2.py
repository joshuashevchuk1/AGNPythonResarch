import numpy as np
import matplotlib.pyplot as plt
import pencil_old as pc
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
dInt = 5


#
# CONTROl ORBITAL INFORMATION HERE!
#

# Orbit=MaxOrbits-2
Orbit = 50

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
    i = 0
    di = 1
    j = 0
    dj = 1
    KE = []
    UE = []
    UINT = []
    FE = []

    while j <= len(phi)-1:
        while i <= len(rad)-1:
            try:
                KE.append(
                    0.5*np.abs(phi[j])*((rad[i+1]**2)-(rad[i]**2))
                    * (ux[j][i]**2+uy[j][j]**2)
                    * 0.5*rho[i][j]
                )
                UE.append(
                    rho[i][j]
                    / rad[i]
                )
                UINT.append(
                    (rho[i][j]
                     * TT[i][j])
                    / gamma
                )
                FE.append(
                    rho[i][j]
                    / rad[i]**2
                )
                i = i+di
            except:

                # out of bound

                i = i+di
        i = 0
        j = j+dj

    KE_Sum.append(np.sum(KE[:]))
    UE_Sum.append(np.sum(UE[:]))
    UINT_Sum.append(np.sum(UINT[:]))
    FE_Sum.append(np.sum(FE[:]))
    Int = Int+dInt

i = 0
di = 1
Total = []

while i <= len(KE_Sum)-1:
    Total.append(KE_Sum[i]+UINT_Sum[i]-UE_Sum[i])
    i = i+di

# ======================================
#   disk energy and Orbital Energy
# ======================================

# get the disk and Orbital energy

torque_Orbit = []

if ecc_int == 0.0:
        #
        #
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
    #
    #
else:
    #
    #
    #	the energy of the orbiter is much more complicated,
    #	the eccentricity change of the planet must be taken into account
    #
    #

    TOR_Sum = []
    ECC_Sum = []
    ECC_Rate = []
    ECC_Ang = []
    SEMI_Sum = []
    RAD_Sum = []
    THE_Sum = []
    VEC_Sum = []
    WE_Sum = []

    i = 0
    di = 1

    while i <= len(Orbit_Len)-2:

        #
        # get the total torque on the planet for every orbit
        #

        TOR_Sum.append(
            torqtotal[
                Orbit_Len[i]
            ]
        )
        ECC_Sum.append(
            ecc[
                Orbit_Len[i]
            ]
        )
        ECC_Rate.append(
            ecc_rate[
                Orbit_Len[i]
            ]
        )
        ECC_Ang.append(
            ecc_ang[
                Orbit_Len[i]
            ]
        )
        SEMI_Sum.append(
            semi_major[
                Orbit_Len[i]
            ]
        )
        RAD_Sum.append(
            radius[
                Orbit_Len[i]
            ]
        )
        THE_Sum.append(
            Omega[
                Orbit_Len[i]
            ]
        )
        VEC_Sum.append(
            np.sqrt(
                radius[Orbit_Len[i]]**2
                +
                Omega[Orbit_Len[i]]**2
            )
        )
        i = i+di
    #
    #
    #

    i = 0
    i = i+di

    print(len(FE_Sum))
    print(len(VEC_Sum))

    while i <= Orbit-1:
        WE_Sum.append(
            VEC_Sum[i]
            *
            FE_Sum[i])
        i = i+di
    #
    #
    #

    OE_Sum = []

    i = 0
    di = 1

    while i <= Orbit-2:
        OE_Sum.append(
            (2/WE_Sum[i])
            *
            (
                (
                    (ECC_Rate[i]*ECC_Sum[i])
                    /
                    (1-ECC_Sum[i]**2)
                )
                +
                (
                    (TOR_Sum[i])
                    /
                    (ECC_Ang[i])
                )


            )
        )
        i = i+di

print("-----------------------------")
print("=====================")
print("DEBUG")
print("=====================")

# print(len(KE))
# print(len(rad))
# print(len(KE_Sum))
# print(KE_Sum)
# print(UE_Sum)
# print(UINT_Sum)
# print(Total)
#print('len information')
# print(len(OE))
# print(len(Orbit_Len))
# print(len(OE_Sum))
print(len(FE_Sum))
print(len(VEC_Sum))

print(TOR_Sum)
if ecc_int == 0.0:
    print('is zero')
else:
    print('is not zero')
# print(OE_Sum)
print("=====================")
print("-----------------------------")


plt.figure(figsize=(10, 10))
plt.title('Energy Comparision')
plt.plot(KE_Sum, color='green', label='Kinetic', ls='--')
plt.plot(UE_Sum, color='red', label='Potential')
plt.plot(UINT_Sum, color='black', label='Internal', ls='-.')
plt.plot(np.abs(Total[:]), color='blue', label='Total')
plt.plot(np.abs(OE_Sum[:]), color='orange', label='Orbital', ls=':')
plt.yscale('log')
plt.xlabel('time (Orbits)', fontweight='bold')
plt.ylabel('Joule (Code units)', fontweight='bold')
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig('Energy_Fixed_2.png')
plt.close()
