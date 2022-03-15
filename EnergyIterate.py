import numpy as np
import matplotlib.pyplot as plt
import pencil as pc

# ======================================
# pc.read_ts()
# ======================================

ts = pc.read_ts()
t = ts.t/2/np.pi
N = 50
tmax = t.max()
ux = ts.ux2m
uy = ts.uy2m
rhomax = ts.rhomax
rhomin = ts.rhomin

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


# ======================================
#   disk energy
# ======================================

# get the disk energy

Int = 0
dInt = 1

Orbit = 3

KE_Sum = []
UE_Sum = []
UINT_Sum = []

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
                i = i+di
            except:

                # out of bound

                i = i+di
        i = 0
        j = j+dj

    KE_Sum.append(np.sum(KE[:]))
    UE_Sum.append(np.sum(UE[:]))
    UINT_Sum.append(np.sum(UINT[:]))
    Int = Int+dInt

print("=====================")
print("DEBUG")
print("=====================")
print(len(KE))
print(len(rad))
print(len(KE_Sum))
print(KE_Sum)
print(UE_Sum)
print(UINT_Sum)
print("=====================")

plt.plot(KE_Sum, color='green')
plt.plot(UE_Sum, color='red')
plt.plot(UINT_Sum, color='orange')
plt.yscale('log')
plt.show()
