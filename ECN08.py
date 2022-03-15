import pencil as pc
import numpy as np
import pylab as plt
import math
import os

HEAD = os.getcwd()
print(HEAD)

# ======================================
# pc.read_pararm()
# ======================================

par = pc.read_param()

ecc_int = par.eccentricity

ts = pc.read_ts()

t = ts.t/2/math.pi

r = ts.xq2
phi = ts.yq2
vr = ts.vxq2
vphi = ts.vyq2

v2 = vr**2 + vphi**2
a = 1./(2/r-v2)
h = r**2*(vphi/r)
ep1 = (h**2)/a
e = (1-ep1)**0.5

N = 100

kernel = np.ones((N,))/N

time = np.convolve(t, kernel, mode='valid')
ecc = np.convolve(e, kernel, mode='valid')

i = 0
di = 1

e_fix = []

while i <= len(e)-1:
    e_fix.append(e[i]+0.02)
    i = i+di

Sigmap = 1.
Omegap = 1./a**1.5
mstar = 8.35e2
q = 1e-4
aspect_ratio = 0.05

i = 0
di = 1

eh = []

while i <= len(e_fix)-1:
    eh.append(e_fix[i]/aspect_ratio)
    i = i+di

i = 0
di = 1

twave = []

while i <= len(a)-1:
    twave.append(mstar * aspect_ratio**4 / (q * Sigmap * a[i]**2 * Omegap))
    i = i+di

i = 0
di = 1
te = []

while i <= len(twave)-1:
    te.append(twave[i]/0.780 * (1 - 0.14*eh[i]**2 + 0.06*eh[i]**3))
    i = i+di

line = ecc_int*np.exp(-t/te)

plt.plot(time, ecc, label='data')
plt.plot(time, np.convolve(line, kernel, mode='valid'),
         label='model (CN08)', linestyle=':')

plt.title(
    r'$\varepsilon$ Decay: $q=10^{-4}$, $h=0.05$, $M_\star = 8.5\times 10^{-2}$')
plt.xlabel(r'$t/T_0$')
plt.ylabel(r'$\varepsilon$')
plt.yscale('log')
plt.xlim([0, t.max()])
plt.tight_layout()
plt.grid(True)
plt.legend()

plt.savefig('eccentricity_decay_fix.png')
# plt.show()
