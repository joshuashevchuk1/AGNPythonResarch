import numpy as np
import matplotlib.pyplot as plt

# attempt to recreate cresswell and nelson plots from simulated data

# e0=0.1

e0 = [0.15, 0.2, 0.25, 0.3, 0.5]


def plot_CN08_SEMI(e0):
    disk_mass = 65
    t = np.arange(0, 500, 0.1)
    sigma0 = 2e-3
    alpha = 0.5
    mstar = 1
    q = 2.4e-11
    aspect_ratio = 0.05
    eh = e0/aspect_ratio
    ap0 = 1
    ap = (-0.001)*t+ap0
    Sigmap = (sigma0*(disk_mass))*ap**-alpha
    Omegap = 1./ap**1.5
    Pe = (1+((eh/2.25)**1.2)+((eh/2.84)**6))/(1-(eh/2.02)**4)
    twave = (mstar * aspect_ratio**4) / (q * Sigmap * ap**2 * Omegap)
    tm = 2*(twave/2.7+1.1*alpha)*(aspect_ratio**-2)*(Pe+(Pe/(np.abs(Pe))))
    CN08 = ap*np.exp(-t/tm)
    plt.plot(t, CN08, label='a'+' = '+str(e0), color='blue')


def plot_CN08(e0):
    disk_mass = 62
    t = np.arange(0, 350, 0.1)
    ap0 = 1
    ap = (-0.001)*t+ap0
    Omegap = 1./ap**1.5
    sigma0 = 2e-3
    alpha = 0.5
    Sigmap = (sigma0*(disk_mass))*ap**-alpha
    mstar = 1
    q = 1e-5
    aspect_ratio = 0.05
    eh = e0/aspect_ratio
    twave = (mstar * aspect_ratio**4) / (q * Sigmap * ap**2 * Omegap)
    te = (twave/0.780) * (1 - 0.14*(eh**2) + 0.06*(eh**3))
    CN08 = e0*np.exp(-t/te)
    print('=================')
    print('twave')
    print(twave)
    print('=================')
    plt.plot(t, CN08, label=r'$\epsilon$'+' = '+str(e0), color='red', ls=':')


def plot_powerlaw(e0):
    k = 0.5
    t = np.arange(0, 300, 1)
    powerlaw = 1/(k*t+1/e0)
    plt.plot(t, powerlaw, label=r'$\epsilon$'+' = '+str(e0), color='red')
    return


for index in range(len(e0)):
    # plot_CN08_SEMI(e0[index])
    plot_CN08(e0[index])

plt.legend()
plt.savefig('EC08_CHECK.png')
