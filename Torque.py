import pencil_old as pc
import pylab as plt
import numpy as np
import math
import os

HEAD = os.getcwd()

print(HEAD)

# Orbits=int(input('Orbits?:'))

ts = pc.read_ts()

t = ts.t/2/math.pi

N = 50

kernel = np.ones((N,))/N

time = np.convolve(t, kernel, mode='valid')

torqint = np.convolve(ts.torqint_2, kernel, mode='valid')

torqext = np.convolve(ts.torqext_2, kernel, mode='valid')


par = pc.read_param()

h = par.cs0

if (par.iprimary == 1):

    q = par.pmass[0]

else:

    q = par.pmass[1]

Gamma0 = (q/h)**2

alpha = par.density_power_law

beta = par.temperature_power_law

#tanaka = (-0.85- alpha - 0.9*beta)*Gamma0
#tanaka = (-1.160)*((1e-4/0.05)**2)
tanaka = 0

plt.plot(time, torqint, '--', label='Inner')

plt.plot(time, torqext, '--', label='Outer')

plt.plot(time, torqext+torqint, label='Total')

plt.plot(time, np.repeat(tanaka, len(time)),
         linestyle=':', label='Tanaka et al. 2002')

plt.legend()


plt.title(r'Torques: Mass Ratio ,'+str(8e-6)+', Aspect ratio $h=0.03$')

plt.xlabel(r'$t/T_0$')

plt.ylabel(r'$\Gamma$')

plt.xlim([0, t.max()])


plt.tight_layout()

plt.grid(True)

plt.savefig('TorqueSGDEV.png')

print('len(t)',len(t))
print('len(ts.torqint_2)',len(ts.torqint_2))
print('len(ts.torqext_2)',len(ts.torqext_2))

# plt.show()
