#---------------------------------------------------------------------------------
#  Fourth-order Theory of Figures gravity calculator
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
from __future__ import division
import sys
import numpy as np
import warnings

def tof4(zvec, dvec, mrot, tol=1e-6, maxiter=100, sskip=0):
    """Return gravity coefficients of mass distribution in hydrostatic equilibrium.

    Parameters
    ----------
    zvec : ndarray, 1d, positive
        Mean radii of level surfaces where density is specified.
    dvec : ndarray, 1d, positive
        Densities on level surfaces. The density should be monotonically
        non-increasing with z, but this is not enforced.
    mrot : float, scalar, nonnegative
        Dimensionless rotation parameter.
    sskip : integer, nonnegative
        skip-n-spline step size
    """

    # Minimal input control (some day)
    zvec = np.array(zvec)
    dvec = np.array(dvec)
    assert zvec.shape == dvec.shape
    assert mrot >= 0
    assert sskip >= 0
    pass

    # Initialize local variables
    p = np.argsort(zvec)
    zvec = zvec[p]
    dvec = dvec[p]
    if zvec[0] == 0:
        zvec[0] = np.spacing(1)
    N = len(zvec)
    ss = 5*[np.zeros(N)]

    # Normalize radii and densities (it's safe to normalize a normal)
    dro = np.hstack((dvec[-1], np.diff(np.flipud(dvec))))
    m = sum(dro*np.flipud(zvec)**3)
    robar = m/zvec[-1]**3
    zvec = zvec/zvec[-1]
    dvec = dvec/robar

    # The loop, following Nettelmann (2017) Appendix B
    Js = np.array([0, 0, 0, 0, 0]) # J0=0 ensures at least one iteration
    for iter in range(maxiter):
        # Equations B.16-B.17
        fs = B1617(ss)

        # Equation B.9
        SS = B9(zvec, dvec, fs)

        # And finally, the system of simultaneous equations B.12-B.15.
        if sskip == 0:
            ss = solve_B1215(ss, SS, mrot)
        else:
            ss = skipnspline_B1215(ss, SS, mrot, zvec, sskip)

        # Now the Js, by eqs. B.1 and B.11
        new_Js, a0 = B111(ss, SS)

        # Check for convergence to terminate
        dJs = np.abs(Js - new_Js)/np.abs(Js + np.spacing(1))
        if np.all(dJs < tol):
            break
        elif iter < maxiter:
            Js = new_Js

    if iter == (maxiter - 1):
        warnings.warn('Figure functions may not be fully converged.')

    # Return
    Js = new_Js
    class out:
        pass
    out.dJs = dJs
    out.iter = iter
    out.a0 = a0
    out.qrot = mrot*a0**3
    out.ss = ss
    out.SS = SS
    return (Js, out)

def mcumtrapz(X, Y):
    # Convert scipy.integrate.cumtrapz to MATLAB-style cumtrapz.
    from scipy.integrate import cumtrapz
    return cumtrapz(Y, X, initial=0)

def B111(ss, SS):
    # Nettelmann 2017 eqs. B.1 and B.11
    N = len(ss[0]) - 1
    s0 = ss[0][N]; s2 = ss[1][N]; s4 = ss[2][N]; s6 = ss[3][N]; s8 = ss[4][N]
    S0 = SS[0][N]; S2 = SS[1][N]; S4 = SS[2][N]; S6 = SS[3][N]; S8 = SS[4][N]
    aos = 1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8
    J0 = -(aos**-0)*S0
    J2 = -(aos**-2)*S2
    J4 = -(aos**-4)*S4
    J6 = -(aos**-6)*S6
    J8 = -(aos**-8)*S8
    Js = np.array([J0, J2, J4, J6, J8])
    return (Js, aos)

def B9(Z, D, fs):
    # Nettelmann 2017 eq. B.9.
    f0 = fs[0]; f2 = fs[1]; f4 = fs[2]; f6 = fs[3]; f8 = fs[4]
    f0p = fs[5]; f2p = fs[6]; f4p = fs[7]; f6p = fs[8]; f8p = fs[9]
    N = len(Z) - 1

    I0 = mcumtrapz(D, Z**(0+3)*f0)
    S0 = D*f0 - Z**-(0+3)*I0

    I2 = mcumtrapz(D, Z**(2+3)*f2)
    S2 = D*f2 - Z**-(2+3)*I2

    I4 = mcumtrapz(D, Z**(4+3)*f4)
    S4 = D*f4 - Z**-(4+3)*I4

    I6 = mcumtrapz(D, Z**(6+3)*f6)
    S6 = D*f6 - Z**-(6+3)*I6

    I8 = mcumtrapz(D, Z**(8+3)*f8)
    S8 = D*f8 - Z**-(8+3)*I8

    I0p = mcumtrapz(D, Z**(2-0)*f0p)
    I0p = I0p[N] - I0p
    S0p = -D*f0p + Z**-(2-0)*(D[N]*f0p[N] - I0p)

    I2p = mcumtrapz(D, Z**(2-2)*f2p)
    I2p = I2p[N] - I2p
    S2p = -D*f2p + Z**-(2-2)*(D[N]*f2p[N] - I2p)

    I4p = mcumtrapz(D, Z**(2-4)*f4p)
    I4p = I4p[N] - I4p
    S4p = -D*f4p + Z**-(2-4)*(D[N]*f4p[N] - I4p)

    I6p = mcumtrapz(D, Z**(2-6)*f6p)
    I6p = I6p[N] - I6p
    S6p = -D*f6p + Z**-(2-6)*(D[N]*f6p[N] - I6p)

    I8p = mcumtrapz(D, Z**(2-8)*f8p)
    I8p = I8p[N] - I8p
    S8p = -D*f8p + Z**-(2-8)*(D[N]*f8p[N] - I8p)

    SS = [S0, S2, S4, S6, S8, S0p, S2p, S4p, S6p, S8p]
    return SS

def B1617(ss):
    """Nettelmann 2017 eqs. B.16 and B.17."""
    s0 = ss[0]; s2 = ss[1]; s4 = ss[2]; s6 = ss[3]; s8 = ss[4]

    f0 = np.ones(len(s0))

    f2 = (3/5)*s2 + (12/35)*s2**2 + (6/175)*s2**3 + (24/35)*s2*s4 + \
         (40/231)*s4**2 + (216/385)*s2**2*s4 - (184/1925)*s2**4

    f4 = (1/3)*s4 + (18/35)*s2**2 + (40/77)*s2*s4 + (36/77)*s2**3 + \
         (90/143)*s2*s6 + (162/1001)*s4**2 + (6943/5005)*s2**2*s4 + \
         (486/5005)*s2**4

    f6 = (3/13)*s6 + (120/143)*s2*s4 + (72/143)*s2**3 + (336/715)*s2*s6 + \
         (80/429)*s4**2 + (216/143)*s2**2*s4 + (432/715)*s2**4

    f8 = (3/17)*s8 + (168/221)*s2*s6 + (2450/7293)*s4**2 + \
         (3780/2431)*s2**2*s4 + (1296/2431)*s2**4

    f0p = (3/2) - (3/10)*s2**2 - (2/35)*s2**3 - (1/6)*s4**2 - \
          (6/35)*s2**2*s4 + (3/50)*s2**4

    f2p = (3/5)*s2 - (3/35)*s2**2 - (6/35)*s2*s4 + (36/175)*s2**3 - \
          (10/231)*s4**2 - (17/275)*s2**4 + (36/385)*s2**2*s4

    f4p = (1/3)*s4 - (9/35)*s2**2 - (20/77)*s2*s4 - (45/143)*s2*s6 - \
          (81/1001)*s4**2 + (1/5)*s2**2*s4

    f6p = (3/13)*s6 - (75/143)*s2*s4 + (270/1001)*s2**3 - (50/429)*s4**2 + \
          (810/1001)*s2**2*s4 - (54/143)*s2**4 - (42/143)*s2*s6

    f8p = (3/17)*s8 - (588/1105)*s2*s6 - (1715/7293)*s4**2 + \
          (2352/2431)*s2**2*s4 - (4536/12155)*s2**4

    fs = [f0, f2, f4, f6, f8, f0p, f2p, f4p, f6p, f8p]
    return fs

def solve_B1215(ss0, SS, mrot):
    # Solve the system B.12-B.15 for unknowns s2,s4,s6,s8.
    from scipy.optimize import fsolve

    N = len(ss0[0])
    Y = np.zeros((N,5))
    for k in range(N):
        es0 = [x[k] for x in ss0[1:]]
        ES = [x[k] for x in SS]
        es = fsolve(B1215, es0, args=(ES,mrot))
        XX = np.concatenate((np.array((0.0, )), es))
        XX[0] = -1/5*XX[1]**2 - 2/105*XX[1]**3 - 1/9*XX[2]**2 - 2/35*XX[1]**2*XX[2]
        Y[k,:] = XX
        pass

    s0 = Y[:,0]; s2 = Y[:,1]; s4 = Y[:,2]; s6 = Y[:,3]; s8 = Y[:,4]
    ss = [s0, s2, s4, s6, s8]
    return ss

def skipnspline_B1215(ss0, SS, mrot, zvec, sskip):
    # Solve the system B.12-B.15 for unknowns s2,s4,s6,s8.
    from scipy.optimize import fsolve
    from scipy.interpolate import interp1d

    N = len(ss0[0])
    ind = np.arange(0, N, sskip)
    Y = np.zeros((len(ind),5))
    for k in range(len(ind)):
        es0 = [x[ind[k]] for x in ss0[1:]]
        ES = [x[ind[k]] for x in SS]
        es = fsolve(B1215, es0, args=(ES,mrot))
        XX = np.concatenate((np.array((0.0, )), es))
        XX[0] = -1/5*XX[1]**2 - 2/105*XX[1]**3 - 1/9*XX[2]**2 - 2/35*XX[1]**2*XX[2]
        Y[k,:] = XX
        pass

    X = zvec[ind]
    s0 = interp1d(X, Y[:,0], 'cubic', fill_value='extrapolate')(zvec)
    s2 = interp1d(X, Y[:,1], 'cubic', fill_value='extrapolate')(zvec)
    s4 = interp1d(X, Y[:,2], 'cubic', fill_value='extrapolate')(zvec)
    s6 = interp1d(X, Y[:,3], 'cubic', fill_value='extrapolate')(zvec)
    s8 = interp1d(X, Y[:,4], 'cubic', fill_value='extrapolate')(zvec)
    ss = [s0, s2, s4, s6, s8]
    return ss

def B1215(s, S, m):
    # Compute the RHS of B.12-B.15.

    s2 = s[0]; s4 = s[1]; s6 = s[2]; s8 = s[3] # s0 not needed
    S0 = S[0]; S2 = S[1]; S4 = S[2]; S6 = S[3]; S8 = S[4]
    S2p = S[6]; S4p = S[7]; S6p = S[8]; S8p = S[9] # S0p not needed

    # B.12
    A2 = 0
    A2 = A2 + S0*(-1*s2 + 2/7*s2**2 + 4/7*s2*s4 - 29/35*s2**3 + 100/693*s4**2 +
                  454/1155*s2**4 - 36/77*s2**2*s4)
    A2 = A2 + S2*(1 - 6/7*s2 - 6/7*s4 + 111/35*s2**2 - 1242/385*s2**3 + 144/77*s2*s4)
    A2 = A2 + S4*(-10/7*s2 - 500/693*s4 + 180/77*s2**2)
    A2 = A2 + S2p*(1 + 4/7*s2 + 1/35*s2**2 + 4/7*s4 - 16/105*s2**3 + 24/77*s2*s4)
    A2 = A2 + S4p*(8/7*s2 + 72/77*s2**2 + 400/693*s4)
    A2 = A2 + m/3*(-1 + 10/7*s2 + 9/35*s2**2 - 4/7*s4 + 20/77*s2*s4 - 26/105*s2**3)

    # B.13
    A4 = 0
    A4 = A4 + S0*(-1*s4 + 18/35*s2**2 - 108/385*s2**3 + 40/77*s2*s4 +
                  90/143*s2*s6 + 162/1001*s4**2 + 16902/25025*s2**4 -
                  7369/5005*s2**2*s4)
    A4 = A4 + S2*(-54/35*s2 - 60/77*s4 + 648/385*s2**2 - 135/143*s6 +
                  21468/5005*s2*s4 - 122688/25025*s2**3)
    A4 = A4 + S4*(1 - 100/77*s2 - 810/1001*s4 + 6368/1001*s2**2)
    A4 = A4 + S6*(-315/143*s2)
    A4 = A4 + S2p*(36/35*s2 + 108/385*s2**2 + 40/77*s4 + 3578/5005*s2*s4 -
                   36/175*s2**3 + 90/143*s6)
    A4 = A4 + S4p*(1 + 80/77*s2 + 1346/1001*s2**2 + 648/1001*s4)
    A4 = A4 + S6p*(270/143*s2)
    A4 = A4 + m/3*(-36/35*s2 + 114/77*s4 + 18/77*s2**2 - 978/5005*s2*s4 +
                   36/175*s2**3 - 90/143*s6)

    # B.14
    A6 = 0
    A6 = A6 + S0*(-s6 + 10/11*s2*s4 - 18/77*s2**3 + 28/55*s2*s6 + 72/385*s2**4 +
                  20/99*s4**2 - 54/77*s2**2*s4)
    A6 = A6 + S2*(-15/11*s4 + 108/77*s2**2 - 42/55*s6 - 144/77*s2**3 + 216/77*s2*s4)
    A6 = A6 + S4*(-25/11*s2 - 100/99*s4 + 270/77*s2**2)
    A6 = A6 + S6*(1 - 98/55*s2)
    A6 = A6 + S2p*(10/11*s4 + 18/77*s2**2 + 36/77*s2*s4 + 28/55*s6)
    A6 = A6 + S4p*(20/11*s2 + 108/77*s2**2 + 80/99*s4)
    A6 = A6 + S6p*(1 + 84/55*s2)
    A6 = A6 + m/3*(-10/11*s4 - 18/77*s2**2 + 34/77*s2*s4 + 82/55*s6)

    # B.15
    A8 = 0
    A8 = A8 + S0*(-1*s8 + 56/65*s2*s6 + 72/715*s2**4 + 490/1287*s4**2 - 84/143*s2**2*s4)
    A8 = A8 + S2*(-84/65*s6 - 144/143*s2**3 + 336/143*s2*s4)
    A8 = A8 + S4*(-2450/1287*s4 + 420/143*s2**2)
    A8 = A8 + S6*(-196/65*s2)
    A8 = A8 + S8*(1)
    A8 = A8 + S2p*(56/65*s6 + 56/143*s2*s4)
    A8 = A8 + S4p*(1960/1287*s4 + 168/143*s2**2)
    A8 = A8 + S6p*(168/65*s2)
    A8 = A8 + S8p*(1)
    A8 = A8 + m/3*(-56/65*s6 - 56/143*s2*s4)

    A = [A2, A4, A6, A8]
    return A

def _test():
    N = 2048
    zvec = np.linspace(1, 1.0/N, N)
    dvec = np.linspace(1/N,2,N)
    mrot = 0.1
    sskip = 32
    Js, out = tof4(zvec, dvec, mrot, 1e-5, sskip=sskip)
    print("After {} iterations:".format(out.iter+1))
    print("J0 = {}".format(Js[0]))
    print("J2 = {}".format(Js[1]))
    print("J4 = {}".format(Js[2]))
    print("J6 = {}".format(Js[3]))
    print("J8 = {}".format(Js[4]))
    print("q = {}".format(out.qrot))
    print("")

if __name__ == '__main__':
    _test()
    sys.exit(0)
