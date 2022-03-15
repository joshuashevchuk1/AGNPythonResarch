import pencil as pc
import numpy as np
import matplotlib.pyplot as plt
import math

ts = pc.read_ts()

t = ts.t/2/np.pi

r = ts.xq2
phi = ts.yq2
vr = ts.vxq2
vphi = ts.vyq2

v2 = vr**2 + vphi**2
a = 1./(2/r-v2)
h = r**2*(vphi/r)
ep1 = (h**2)/a
e = (1-ep1)**0.5

# print(e)

N = 50

kernel = np.ones((N,))/N

time = np.convolve(t, kernel, mode='valid')
ecc = np.convolve(e, kernel, mode='valid')

semi_major = np.convolve(a, kernel, mode='valid')

radius = np.convolve(r, kernel, mode='valid')

# true_angle=(1.0/e)*(((a*(1.0-e**2.0))/r)-1.0)
true_angle = (1.0/ecc)*(((semi_major*(1.0-ecc**2.0))/radius)-1.0)

f = np.arccos(true_angle)

# f=((math.pi/180)*f)

phi = np.convolve(phi, kernel, mode='valid')

#f   = np.convolve(f,kernel,mode='valid')

angle = np.arctan(radius/phi)

argw = phi-f

# print(argw)
# print(phi)
# print(r)
print(f[:])
# print(phi[:])
print(angle[:])
print(argw)

plt.plot(time, argw)
plt.show()
