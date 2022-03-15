import numpy as np
import matplotlib.pyplot as plt
import pencil_old as pc

ff=pc.read_var(trimall=True,ivar=3,magic=["TT"])
f0=pc.read_var(trimall=True,ivar=0,magic=["TT"])

dfT=ff.TT-f0.TT

plt.plot(ff.x,np.log10(np.sum(dfT**2,axis=0)))
plt.show()
