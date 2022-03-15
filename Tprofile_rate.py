import pencil as pc
import numpy as np
import matplotlib.pyplot as plt

ts = pc.read_ts()
t = ts.t
time = t/2*np.pi
len_t = len(ts.t)

Max_Orbits = np.round((time[len(time)-1]))

n = 1
dn = 1

i = 0
di = 1

Orbit_Len = []

while i <= len(time)-1:
    if np.round(n*2.0*np.pi) == np.round(t[i]):
        Orbit_Len.append(len(t[:i]))
        n = n+dn
        i = i+di
    else:
        i = i+di


TTm = ts.TTm

# TTm_rate=np.mean(np.gradient(TTm[:Orbit_Len[300]]))
TTm_rate = (TTm[Orbit_Len[300]]-TTm[0])/(t[300]-t[0])

plt.plot(time, TTm, label='rate: '+str(TTm_rate), color='red')
plt.xlabel('orbits')
plt.xlim([0, 300])
plt.ylim([0.0040, 0.0060])
plt.ylabel('temperature')
plt.legend()
plt.show()
