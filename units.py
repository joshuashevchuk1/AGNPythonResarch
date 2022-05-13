import numpy as np
import math

#
#  AGN physical units, Sirko & Goodman model.
#

# Constants

AU = 1.49e13
pc = 3.086e18
G = 6.67408e-8
Msun = 1.99e33
c = 2.998e10
yr = 3.154e7
k = 1.380648e-16
pi = np.pi

#  Input

M = 1e8*Msun    # Central super-massive black hole

#  Gravitational radius

rg = (G/c**2) * M

r0 = 1e5*rg  # Reference radius
aspect_ratio = 0.05     # H/r (in code units this is equal to cs0)

# Calculate units

unit_length = r0
g0 = G*M
Omega = np.sqrt(g0)*r0**(-1.5)  # Keplerian angular frequency
unit_time = 1./Omega
unit_velocity = unit_length/unit_time

print("Units")
print('unit_length=', unit_length, ' cm')
print('unit_time=', unit_time, ' s')
print('orbital period =', 2*math.pi/Omega/yr, ' yrs')
print('unit_velocity=', unit_velocity, ' cm/s')
print('unit_velocity/c=', unit_velocity/c)

#number_density_carbon = 100
#hydrogen_to_carbon_ratio = 1e4
#number_density = hydrogen_to_carbon_ratio * number_density_carbon

atomic_mass_unit = 1.66054e-24
mean_molecular_weight = 2.5
#unit_density =number_density * mean_molecular_weight  * atomic_mass_unit

unit_density = 1e-7  # g/cm**3
surface_density = ((unit_density)/(r0))
# surface_density=3e-5

print('unit_density=', unit_density, ' g/cm3')
print(" ")

H = aspect_ratio * r0   # pressure scale height
cs = Omega * H
gamma = 1

T = cs**2 * (mean_molecular_weight * atomic_mass_unit)/gamma/k

#disk_mass = unit_density*(pi*(25000*AU)**2)
disk_mass = unit_density*(pi*(25000*AU)**2)
Q = ((cs)*Omega)/(pi*G*surface_density)

print("Thermal quantities")
print("Scale height H=", H/AU, " AU")
print("Sound Speed, [m/s]=", cs/1e2)
print("cs/c=", cs/c)
print("Temperature [1e6 K]=", T / 1e6)
print("expected Toomre Q at r0", Q)
print("log Q at r0" ,np.log(Q))
print("disk mass", disk_mass)
print("star mass", M)
print("disk mass/star mass", disk_mass/M)
print("rg in AU", rg/AU)
print(unit_length/pc)
print("r0 in pc",r0/pc)
print("rmin in pc",1e4*rg/pc)
print("delta r is ",r0/pc - 1e4*rg/pc)
print("log unit density is ",np.log(unit_density))
