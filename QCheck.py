import numpy as np
import pencil_old as pc
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
# ------------------
#
#   Debug script to test check Toomre Q parameter
#
# ------------------


# ======================================
# pc.read_ts()
# ======================================

t = ts.t/2/np.pi
N = 50
tmax = t.max()

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
dInt = 30

#
# CONTROl ORBITAL INFORMATION HERE!
#

# Orbit=MaxOrbits-2
Orbit = 300

#
#

UE_Sum = []
Q_part = []
disk_angular = []

while Int <= Orbit:

    ff = pc.read_var(trimall=True, ivar=Int, magic=["TT"])
    TT = ff.TT
    ss = ff.ss
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
    UE = []
    ss_R = []
    dang = []

    while j <= len(phi)-1:
        while i <= len(rad)-1:
            try:
                dang.append(
                    uy[j][j]
                )
                UE.append(
                    rho[i][j]
                    / rad[i]
                )
                ss_R.append(
                    np.exp(ss[i][j])
                )
                i = i+di
            except:
                # out of bound
                i = i+di
        i = 0
        j = j+dj

    UE_Sum.append(np.mean(UE[:]))
    Q_part.append(np.mean(ss_R[:]))
    disk_angular.append(np.mean(ss_R[:]))
    Int = Int+dInt

Toomre = []

for part in range(len(Q_part)):
    Toomre.append((
        Q_part[part]
        *
        disk_angular[part]
    )
        /
        UE_Sum[part])

plt.plot(Toomre)
plt.ylabel('Q eff')
plt.xlabel(r'$\tau$')
plt.title('mean Toomre Q eff per orbit')
plt.savefig('Q_Check')
# plt.show()
