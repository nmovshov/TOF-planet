#------------------------------------------------------------------------------
# INTERIOR STRUCTURE OF A NON-ROTATING INDEX 1 POLYTROPE
# Example and test of the TOFPlanet class. We construct and converge a
# model of a non-rotating fluid planet with the pressure-density law
#
# $$P = K\rho^2$$
#
# with a polytropic constant $K$. This script demonstrates how to set up a
# TOFPlanet object with a barotrope and converge to a density structure in
# hydrostatic equilibrium. The default starting density profile is that of a
# homogeneous sphere. This scripts also explores the error due to the radius
# grid resolution.
#
# Reminder: the density as a function of radius of a non-rotating, index-1
# polytropic planet in hydrostatic equilibrium is
#
# $$\rho(r) = \rho_c \sin(ar)/(ar)$$
#
# where
#
# $$a = \sqrt{2\pi{G}/K}$$
#
# and
#
# $$\rho_c = 3M/(4\pi{R^3})$$
#
# where $M$ is the planet's mass and $R = \pi/a$, which is independent of mass.
#------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('..')
from TOFPlanet import TOFPlanet

G = 6.67430e-11

# Set up planets with arbitrary mass and radius
M = 317.8*5.9722e+24
R = 71492e3

# Construct a polytrope of index 1
# A non-rotating, index-1 polytrope is completely defined by K. The radius is
# independent of mass. The density structure normalized to the central density
# is also independent of mass, and the absolute value of density is determined
# with the known average density.
n = 1 # polytrope index
K = 2*G/np.pi*R**2 # polytrope constant
def poly1(P): return np.sqrt(P/K)

N = 2**np.arange(12,16)
TPS = []
for n in N:
    tp = TOFPlanet()
    tp.name = f'tp{n}'
    tp.G = G
    tp.mass = M
    tp.GM = G*M
    tp.radius = R
    tp.si = R*np.linspace(1, 1/n, n)
    tp.rhoi = np.ones(n)*M/(4*np.pi/3*R**3)
    tp.period = np.inf
    tp.P0 = 0.0
    tp.set_barotrope(poly1)
    TPS.append(tp)


# Relax to barotrope
for tp in TPS:
    tp.opts['dJtol'] = 1e-10
    tp.relax_to_barotrope()

# Compare computed and analytic density structure
a = np.sqrt(2*np.pi*G/K)
R = np.pi/a
rho_av = 3*M/(4*np.pi*R**3)
rho_c = (np.pi**2/3)*rho_av
r = np.linspace(0,1)*R
rho_exact = rho_c*np.sin(a*r)/(a*r)
rho_exact[0] = rho_c

# Plot
plt.plot(r/R, rho_exact/rho_c, 'k--', label=r'$\sin(ar)/(ar)$')
for tp in TPS:
    plt.plot(tp.si/R, tp.rhoi/rho_c, '-', label=tp.name)

plt.xlabel(r'$r/R$')
plt.ylabel(r'$\rho/\rho_c$')

plt.legend()
plt.show()
