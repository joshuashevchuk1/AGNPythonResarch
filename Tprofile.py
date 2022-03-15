import pencil as pc
import numpy as np
import matplotlib.pyplot as plt

ts = pc.read_ts()

TTm = ts.TTm

plt.plot(TTm)
plt.show()
