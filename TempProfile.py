import pencil as pc
import numpy as np
import matplotlib.pyplot as plt

ts = pc.read_ts()

t = ts.t/2*np.pi

temp = ts.TTm
temp_max = ts.TTmax
temp_min = ts.TTmin

plt.plot(t, temp)
plt.plot(t, temp_max)
plt.plot(t, temp_min)
plt.show()
