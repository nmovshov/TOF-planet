#------------------------------------------------------------------------------
# GRAVITY COEFFICIENTS OF ROTATING INDEX-1 POLYTROPE
#
# Example and test of the TOFPlanet class. We construct and converge a model of
# a rotating fluid planet with the pressure-density law
#
# $$P = K\rho^2$$
#
# with a polytropic constant $K$. This script demonstrates how to set up a
# TOFPlanet with a barotrope and converge to hydrostatic equilibrium. The
# default starting density profile is that of a homogeneous sphere.
#
# The output is a benchmarking table. The numbers to aim for are apparently
# those of the CLC method (Wisdom, 1996, unpublished) which is lost to time. I
# take the numbers from Wisdom and Hubbard (2016) Table 3. For consistency I
# try to replicate the Wisdom and Hubbard model exactly, including using their
# values for G and K (Guillot, double hearsay) and trying to match their q with
# a rotation period.
#
# Also comparing with CMS results from Wisdom and Hubbard (2016) Table 4, and
# with some q powers I dug out of the ZT red book that I think are supposed to
# match a 3rd-order ToF solution. Also adding some more numbers I got from
# Nadine (Nettelmann 2021 private communication).
#------------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd
sys.path.append('..')
from TOFPlanet import TOFPlanet

# Construct a polytrope of index 1, aiming for exact replication of WH16
G = 6.6738480e-11 # Hubbard to Guillot personal communication
GM = 1.266865361e17 # WH16
M = GM/G
Re = 71492*1e3 # (to match K use K=2*G/pi*R**2 instead)
qrot = 0.089195487 # WH16
wrot = np.sqrt(qrot*GM/Re**3)
Prot = 2*np.pi/wrot
aos = 1.022875431133185 # WH16 Table 3 Re/R
K = 2.003565e5 # Hubbard to Guillot personal communication (no effect on Js)
n = 1
def poly1(P): return np.sqrt(P/K)


# Set up TOFPlanet(s)
N = 2**12 # around N=2^18 ToF7 gets more precise J2 than ToF4
nx = 128

tofour = TOFPlanet(toforder=4)
tofour.opts['verbosity'] = 1
tofour.name = f'{N}-point TOF4'
tofour.G = G # undocumented TOFPlanet property
tofour.GM = GM # undocumented TOFPlanet property
tofour.mass = M
tofour.radius = Re
tofour.period = Prot # trying to match WH16 qrot
tofour.si = Re*np.linspace(1, 1/N, N) # % will be renormalized
tofour.rhoi = np.ones(N)*M/(4*np.pi/3*Re**3) # initial guess
tofour.P0 = 0 # added to surface pressure
tofour.set_barotrope(poly1)

tofsev = TOFPlanet(toforder=7)
tofsev.opts['verbosity'] = 1
tofsev.name = f'{N}-point TOF7'
tofsev.G = G # undocumented TOFPlanet property
tofsev.GM = GM # undocumented TOFPlanet property
tofsev.mass = M
tofsev.radius = Re
tofsev.period = Prot # trying to match WH16 qrot
tofsev.si = Re*np.linspace(1, 1/N, N) # will be renormalized
tofsev.rhoi = np.ones(N)*M/(4*np.pi/3*Re**3) # initial guess
tofsev.P0 = 0 # added to surface pressure
tofsev.set_barotrope(poly1)

# Relax to barotrope (fast for tof4, slow for tof7)
tofour.opts['xlevels'] = nx
tofour.opts['drhotol'] = 1e-6
tofour.opts['dJtol'] = 1e-10
tofour.opts['MaxIterBar'] = 60
tofour.opts['MIXL'] = 0
tofour.relax_to_barotrope()

tofsev.opts['xlevels'] = nx
tofsev.opts['drhotol'] = 1e-6
tofsev.opts['dJtol'] = 1e-10
tofsev.opts['MaxIterBar'] = 60
tofsev.opts['MIXL'] = 0
tofsev.relax_to_barotrope()

# Construct the benchmarking table
# The variables to compare are [Re/R, J2, J4, ..., J14]

# With TOFPlanet
MTOF4 = [tofour.a0/tofour.s0, *tuple(tofour.Js[1:]), np.nan, np.nan, np.nan]
MTOF7 = [tofsev.a0/tofsev.s0, *tuple(tofsev.Js[1:])]

# Wisdom and Hubbard (2016) Table 3
CLC = [1.022875431133185, 1.398851089834637e-2, -5.318281001092471e-4,
                          3.011832290533577e-5, -2.132115710726158e-6,
                          1.740671195871128e-7, -1.568219505602588e-8,
                          1.518099230068580e-9]

# Wisdom and Hubbard (2016) Table 4
CMS512 = [np.nan,         1.398924011471443e-2, -5.318792055591143e-4,
                          3.012230402792236e-5, -2.132458660888379e-6,
                          1.740988882124251e-7, -1.568529225158399e-8,
                          1.518412099478797e-9]

# % Hubbard (2013) Table 5
H13 =    [np.nan, 1.3989253e-2, -5.3187997e-4, 3.0122356e-5, -2.1324628e-6,
                  1.7409925e-7, -1.5685327e-8, 1.5184156e-9]

# % Zharkov & Trubistyn (1978) eq. 34.12
q = tofour._P2q()
ZT3 = [np.nan, (0.173273*q - 0.197027*q**2 + 0.15*q**3),
               (-0.081092*q**2 + 0.15*q**3),
               (0.056329*q**3), np.nan, np.nan, np.nan, np.nan]

# Now the table of comparisons, using absolute diff from WH16-3
cols = ['Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14']
rows = ['CLC', 'TOF4', 'TOF7', 'CMS_WH16', 'CMS_H13', 'ZT78']
A = np.array([CLC, MTOF4, MTOF7, CMS512, H13, ZT3])
E = np.abs(A[1:,:] - A[0,:])
# T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
WH16_diffs = pd.DataFrame(E, columns=cols, index=rows[1:])

# Output
pd.options.display.float_format = '{:.4e}'.format
print(WH16_diffs)
