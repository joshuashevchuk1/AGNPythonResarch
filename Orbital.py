import pencil as pc
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scp
import os

global radq1, thetaq1, t, xrq1, yrq1, xrq2, yrq2, vxq1, vyq1, vxq2, vyq2, r, vq1, vq2, v, G, ap1, ap2, ap3, a, rmax, e, ep1, v, xq2, yq2, xq1, yq1, vx, vy

TIMESERIES = pc.read_ts()

radq1 = TIMESERIES.xq1
thetaq1 = TIMESERIES.yq1
t = TIMESERIES.t

xrq1 = radq1*np.cos(thetaq1)
yrq1 = radq1*np.sin(thetaq1)

radq2 = TIMESERIES.xq2
thetaq2 = TIMESERIES.yq2

xrq2 = radq2*np.cos(thetaq2)
yrq2 = radq2*np.sin(thetaq2)

vxq1 = TIMESERIES.vxq1
vyq1 = TIMESERIES.vyq1

vxq2 = TIMESERIES.vxq2
vyq2 = TIMESERIES.vyq2

xq2 = TIMESERIES.xq2
yq2 = TIMESERIES.yq2

dt = TIMESERIES.dt

t = TIMESERIES.t

dift = np.diff(t)

G = scp.G

v2 = (vxq2**2)+(vyq2**2)
a = 1/(2/xq2-v2)
h = (xq2**2)*(vyq2/xq2)
ep1 = (h**2)/a
e = (1-ep1)**0.5

ts = pc.read_ts()

t = ts.t/2/np.pi

N = 50

kernel = np.ones((N,))/N

time = np.convolve(t, kernel, mode='valid')
ecc = np.convolve(e, kernel, mode='valid')
smajor = np.convolve(a, kernel, mode='valid')

fig, ((ax2, ax3, ax4)) = plt.subplots(1, 3, figsize=(30, 5))
fig.tight_layout()
fig.subplots_adjust(wspace=0.18)

plt.title('Orbit in Cylin Coords')
ax2.set_aspect('equal')
ax3.set_aspect('auto')
ax4.set_aspect('auto')

ax2.set_title('x vs y')
ax3.set_title('a vs t')
ax4.set_title('e vs t')

ax2.set_ylabel('yrq1')
ax2.plot(xrq2, yrq2)
ax2.set_xlabel('xrq2')
ax2.set_ylabel('yrq2')
ax3.plot(time, smajor, color='orange')
ax3.set_xlabel('t')
ax3.set_ylabel('a')
ax4.plot(time, ecc, color='orange')
ax4.set_xlabel('t')
ax4.set_ylabel('e')

ax2.grid(True)
ax3.grid(True)
ax4.grid(True)

plt.savefig('Orbital_Information.png')

# plt.show()
