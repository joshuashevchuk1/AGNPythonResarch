import pencil as pc
import pylab as plt
import numpy as np

fig, (ax2) = plt.subplots(1, 1, figsize=(10, 10))

ts03 = pc.read.ts()
t03 = ts03.t

par = pc.read.param()

q = 1 - par.pmass[par.iprimary - 1]


def get_ecc(rr, vrad, vphi):
    v2 = vrad ** 2 + vphi ** 2
    a = 1. / (2. / rr - v2)
    phidot = vphi / rr
    h = rr ** 2 * phidot
    e = np.sqrt(1 - h ** 2 / a)
    return a, e


a, e = get_ecc(ts03.xq2, ts03.vxq2, ts03.vyq2)

ax2.plot(t03 / 2 / np.pi, e, label='e')

beta = par.density_power_law

omega = 1 / a ** 1.5
h = par.cs0
twave = 1 / q * 1 / par.rho0 * 1 / a ** 2 * h ** 4 / omega

te = twave / 0.780 * (1 - 0.14 * (e / h) ** 2 + 0.06 * (e / h) ** 3)

ax2.plot(t03 / 2 / np.pi, e[0] * np.exp(-t03 / te), linestyle='--', label='CN08 best fit (Eq 11)')
ax2.set_ylim([0, .3])

pe = 1 + (e / (2.25 * h)) ** 1.2 + (e / (2.84 * h)) ** 6 / (1 - (e / (2.02 * h)) ** 4)
tm = 2 * twave / (2.7 + 1.1 * beta) * 1 / h ** 2 * (pe + pe / np.abs(pe))

ax2.set_title(r'$\varepsilon$' + ' vs t' + ' , ' + 'q=' + str(q),fontsize=30)
ax2.set_xlabel(r'$t$',fontsize=20)
ax2.set_ylabel(r'$\varepsilon$',fontsize=20)
ax2.legend()

indexTimeCutOff = 0
indexTimeCutOffLarge = 0

for i in range(len(e)):
        if e[i] <= 0.01:
            indexTimeCutOff = i
            break

timeCutOff = t03[indexTimeCutOff] / (np.pi * 2)
timeCutOffLarge = t03[indexTimeCutOffLarge] / (np.pi * 2)

ax2.set_xlim([0,t03[indexTimeCutOff]])

plt.grid(True)
plt.savefig("EccDecay_"+max_lim+".png",bbox_inches='tight')
