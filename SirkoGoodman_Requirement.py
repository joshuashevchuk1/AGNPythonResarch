import numpy as np
import units as unit

# Q = cs*Omega/(pi*G*sigma)
unit_density = unit.unit_density
sigma = np.arange(4, 5, 0.1)
sigma = 10**(sigma)
Omega = unit.Omega
pi = np.pi
G = unit.G
cs = np.arange(5, 6, 0.1)
cs = (10**cs)
# cs=unit.cs
Omega = 1

Q = (cs*Omega)/(pi*sigma)
Q = np.arange(0.01, 0.5, 0.1)
Q = (10**Q)

print('Expected Q after given sigma range', Q)
