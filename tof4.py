#------------------------------------------------------------------------------
#  Fourth-order Theory of Figures gravity calculator
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#------------------------------------------------------------------------------
import sys
import numpy as np
import warnings
from scipy.interpolate import interp1d

def tof4(zvec, dvec, mrot, **kwargs):
    """Return gravity coefficients of self-gravitating rotating fluid.

    Required parameters
    -------------------
    zvec : ndarray, 1d, positive
        Mean radii of level surfaces where density is defined.
    dvec : ndarray, 1d, positive
        Density on level surfaces; density should be monotonically
        non-increasing with z (this is *assumed*, not enforced).
    mrot : float, scalar, nonnegative
        Dimensionless rotation parameter, w^2*s0^3/GM.
    Optional parameters
    -------------------
    xlevels : scalar or vector, nonnegative, integer (xlevels=-1)
        Levels whose shape will be explicitly calculated. The shape functions
        will be explicitly calculated for these levels, and spline-interpolated
        in between. This can result in significant speedup with minimal loss of
        precision. A scalar value is interpreted as a number of xlevels to be
        uniformaly distributed among the density levels. A vector value is
        interpreted as indices of levels to be used as xlevels. Skip-n-spline
        is recommended for very high resolution density profiles, N>~10^4.
        Disable skip-n-spline by passing xlevels=-1 rather than xlevels=N, to
        avoid the spline ovearhead.
    tol : scalar, positive, (tol=1e-6)
        Convergence tolerance for relative changes in Js in successive
        iterations (but keep in mind truncation error of ToF).
    maxiter : scalar, positive, integer, (maxiter=100)
        Maximum number of algorithm iterations.
    calc_moi : bool, scalar (default False)
        Flag to calculate and return the normalized moment of inertia in output
        struct; off by default since it takes an extra half second or so.
    ss_guesses : list of length svec.size arrays (default 5*[np.zeros(N)])
        Initial guess for shape functions. This is not all that helpful in
        speeding up convergence. It's occasionally helpful to preserve state
        between successive calls.

    Outputs
    -------
    Js : 1-by-5 vector, real
        Even harmonic gravity coefficients J0 to J8 (J0 is included as a sanity
        check and test of convergence).
    out : struct
        A structure holding other quantities calculated in the process,
        including the shape functions that define the full solution.

    Algorithm
    ---------
    Theory of figures equations and coefficients in Nettelmann 2017 Appendix B.
    """

    # Minimal input control
    zvec = np.array(zvec)
    dvec = np.array(dvec)
    p = np.argsort(zvec) # ToF equations easier to read w inside-out radii
    zvec = zvec[p]
    dvec = dvec[p]
    if zvec[0] == 0:
        zvec[0] = np.spacing(1)
    assert zvec.shape == dvec.shape
    assert mrot >= 0

    # Spread default and passed **kwargs
    opts = {**default_opts(), **kwargs}

    # Normalize radii and densities (it's safe to normalize a normal)
    dro = np.hstack((dvec[-1], np.diff(np.flipud(dvec))))
    m = sum(dro*np.flipud(zvec)**3)
    robar = m/zvec[-1]**3
    zvec = zvec/zvec[-1]
    dvec = dvec/robar

    # Define and initialize local variables
    N = len(zvec)
    if opts['ss_guesses'] is None:
        ss = 5*[np.zeros(N)]
    else:
        ss = opts['ss_guesses']
    if np.isscalar(opts['xlevels']):
        sskip = int(max(np.fix(N/opts['xlevels']), 1))
        xind = range(0, N, sskip)
    else:
        raise(Exception("Vector xlevels feature not yet implemented"))

    # The loop, following Nettelmann (2017) Appendix B
    Js = np.array([0, 0, 0, 0, 0]) # J0=0 ensures at least one iteration
    for it in range(opts['maxiter']):
        # Equations B.16-B.17
        fs = B1617(ss)

        # Equation B.9
        SS = B9(zvec, dvec, fs)

        # And finally, the system of simultaneous equations B.12-B.15.
        ss = skipnspline_B1215(ss, SS, mrot, zvec, xind)

        # Now the Js, by eqs. B.1 and B.11
        new_Js, a0 = B111(ss, SS)

        # Check for J2 convergence to terminate
        dJs = np.abs(Js - new_Js)
        if (it > 0) and (dJs[1] < opts['tol']):
            break
        Js = new_Js

    if (it == (opts['maxiter'] - 1)) and (opts['verbosity'] > 0):
        warnings.warn('Figure functions may not be fully converged.')

    # Return
    Js = new_Js
    class out:
        pass
    out.dJs = dJs
    out.it = it
    out.a0 = a0
    out.qrot = mrot*a0**3
    out.ss = ss
    out.SS = SS
    out.A0 = B4(ss,SS,mrot)
    if opts['calc_moi']:
        out.NMoI = NMoI(zvec, dvec, ss, a0)
    else:
        out.NMoI = None
    return (Js, out)


### Helper functions
def default_opts():
    """Return dict of default values for tof4's **kwargs."""
    opts = {}
    opts['xlevels'] = -1
    opts['tol'] = 1e-10
    opts['maxiter'] = 100
    opts['calc_moi'] = False
    opts['ss_guesses'] = None
    opts['verbosity'] = 1
    return opts

def NMoI(zi, rhoi, ss, a0):
    # Fast moment of inertia calculation (normalized by a0)
    rhoi = np.hstack((rhoi[0], (rhoi[:-1] + rhoi[1:])/2))
    deltas = np.flipud(np.hstack((rhoi[-1], np.diff(np.flipud(rhoi)))))
    num = 0
    den = 0
    mus = np.linspace(0, 1, 101)
    h = mus[1]/2 # for trapz integral
    p2term = 1 - 0.5*(3*mus**2 - 1)
    for k in range(len(zi)):
        if ss is None: # hack for spherically symmetric rho(z)
            ximax = zi[k]*np.ones_like(mus)
        else:
            ximax = (zi[k]/a0)*level_surface(k, mus, ss)
        fun1 = deltas[k]*(ximax**5)*p2term
        fun2 = deltas[k]*(ximax**3)
        num = num + h*(fun1[0] + 2*sum(fun1[1:-1]) + fun1[-1]) # trapezoid rule
        den = den + h*(fun2[0] + 2*sum(fun2[1:-1]) + fun2[-1]) # trapezoid rule

    return (2/5)*(num/den)

def level_surface(k, mus, ss):
    # Normalized r(cos(theta)) on kth-level surface
    s0 = ss[0][k]; s2 = ss[1][k]; s4 = ss[2][k]; s6 = ss[3][k]; s8 = ss[4][k]
    shp = s0*Pn(0,mus) + s2*Pn(2,mus) + s4*Pn(4,mus) + s6*Pn(6,mus) + s8*Pn(8,mus)
    return 1 + shp

def Pn(n, x):
    # Fast implementation of ordinary Legendre polynomials of low even degree.
    if n == 0:
        y = np.ones_like(x)
    elif n == 2:
        y = 0.5*(3*x**2 - 1)
    elif n == 4:
        y = (1/8)*(35*x**4 - 30*x**2 + 3)
    elif n == 6:
        y = (1/16)*(231*x**6 - 315*x**4 + 105*x**2 - 5)
    elif n == 8:
        y = (1/128)*(6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)
    else:
        raise(Exception("Unimplemented order"))
    return y

def mcumtrapz(X, Y):
    # Convert scipy.integrate.cumtrapz to MATLAB-style cumtrapz.
    from scipy.integrate import cumtrapz
    return cumtrapz(Y, X, initial=0)

def B111(ss, SS):
    # Return Js from SS, with Req/Rm a necessary bonus.
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

def skipnspline_B1215(ss0, SS, mrot, zvec, xind):
    # Update the system B.12-B.15 for new s2,s4,s6,s8.

    # Skip
    Zs = np.array([S[xind] for S in SS]).T
    zs = np.array([s[xind] for s in ss0]).T
    newzs = B1215(zs, Zs, mrot)
    newz0 = ((-1/5)*newzs[:,0]**2 - (2/105)*newzs[:,0]**3 -
                (1/9)*newzs[:,1]**2 - 2/35*newzs[:,0]**2*newzs[:,1])
    newz0 = np.reshape(newz0, (newz0.size,1))
    Y = np.hstack((newz0,newzs))

    # And spline
    if len(xind) < len(zvec):
        X = zvec[xind]
        s0 = interp1d(X, Y[:,0], 'cubic', fill_value='extrapolate')(zvec)
        s2 = interp1d(X, Y[:,1], 'cubic', fill_value='extrapolate')(zvec)
        s4 = interp1d(X, Y[:,2], 'cubic', fill_value='extrapolate')(zvec)
        s6 = interp1d(X, Y[:,3], 'cubic', fill_value='extrapolate')(zvec)
        s8 = interp1d(X, Y[:,4], 'cubic', fill_value='extrapolate')(zvec)
    else:
        s0 = Y[:,0]
        s2 = Y[:,1]
        s4 = Y[:,2]
        s6 = Y[:,3]
        s8 = Y[:,4]
    ss = [s0, s2, s4, s6, s8]
    return ss

def B4(ss, SS, m):
    # Compute the RHS of B.4 in Nettelmann (2017).
    s2 = ss[1]; s4 = ss[2]
    S0 = SS[0]; S2 = SS[1]; S4 = SS[2]
    S0p = SS[5]; S2p = SS[6]; S4p = SS[7]

    A0 = np.zeros_like(s2)
    A0 = A0 + S0*(1 + (2/5)*s2**2 - (4/105)*s2**3 + (2/9)*s4**2 + (43/175)*s2**4 - (4/35)*s2**2*s4)
    A0 = A0 + S2*(-(3/5)*s2 + (12/35)*s2**2 - (234/175)*s2**3 + (24/35)*s2*s4)
    A0 = A0 + S4*((6/7)*s2**2 - (5/9)*s4)
    A0 = A0 + S0p*(1)
    A0 = A0 + S2p*((2/5)*s2 + (2/35)*s2**2 + (4/35)*s2*s4 - (2/25)*s2**3)
    A0 = A0 + S4p*((4/9)*s4 + (12/35)*s2**2)
    A0 = A0 + (m/3)*(1 - (2/5)*s2 - (9/35)*s2**2 - (4/35)*s2*s4 + (22/525)*s2**3)

    return A0

def B1215(s, S, m):
    # Compute the RHS of B.12-B.15 and "solve" for sn.

    s2 = s[:,1]; s4 = s[:,2]; s6 = s[:,3]; # s0 and s8 not needed
    S0 = S[:,0]; S2 = S[:,1]; S4 = S[:,2]; S6 = S[:,3]; S8 = S[:,4]
    S2p = S[:,6]; S4p = S[:,7]; S6p = S[:,8]; S8p = S[:,9] # S0p not needed

    # B.12 (not including -s2S0)
    A2 = 0
    A2 = A2 + S0*(2/7*s2**2 + 4/7*s2*s4 - 29/35*s2**3 + 100/693*s4**2 +
                  454/1155*s2**4 - 36/77*s2**2*s4)
    A2 = A2 + S2*(1 - 6/7*s2 - 6/7*s4 + 111/35*s2**2 - 1242/385*s2**3 + 144/77*s2*s4)
    A2 = A2 + S4*(-10/7*s2 - 500/693*s4 + 180/77*s2**2)
    A2 = A2 + S2p*(1 + 4/7*s2 + 1/35*s2**2 + 4/7*s4 - 16/105*s2**3 + 24/77*s2*s4)
    A2 = A2 + S4p*(8/7*s2 + 72/77*s2**2 + 400/693*s4)
    A2 = A2 + m/3*(-1 + 10/7*s2 + 9/35*s2**2 - 4/7*s4 + 20/77*s2*s4 - 26/105*s2**3)

    # B.13 (not including -s4S0)
    A4 = 0
    A4 = A4 + S0*(18/35*s2**2 - 108/385*s2**3 + 40/77*s2*s4 +
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

    # B.14 (not including -s6S0)
    A6 = 0
    A6 = A6 + S0*(10/11*s2*s4 - 18/77*s2**3 + 28/55*s2*s6 + 72/385*s2**4 +
                  20/99*s4**2 - 54/77*s2**2*s4)
    A6 = A6 + S2*(-15/11*s4 + 108/77*s2**2 - 42/55*s6 - 144/77*s2**3 + 216/77*s2*s4)
    A6 = A6 + S4*(-25/11*s2 - 100/99*s4 + 270/77*s2**2)
    A6 = A6 + S6*(1 - 98/55*s2)
    A6 = A6 + S2p*(10/11*s4 + 18/77*s2**2 + 36/77*s2*s4 + 28/55*s6)
    A6 = A6 + S4p*(20/11*s2 + 108/77*s2**2 + 80/99*s4)
    A6 = A6 + S6p*(1 + 84/55*s2)
    A6 = A6 + m/3*(-10/11*s4 - 18/77*s2**2 + 34/77*s2*s4 + 82/55*s6)

    # B.15 (not including -s8S0)
    A8 = 0
    A8 = A8 + S0*(56/65*s2*s6 + 72/715*s2**4 + 490/1287*s4**2 - 84/143*s2**2*s4)
    A8 = A8 + S2*(-84/65*s6 - 144/143*s2**3 + 336/143*s2*s4)
    A8 = A8 + S4*(-2450/1287*s4 + 420/143*s2**2)
    A8 = A8 + S6*(-196/65*s2)
    A8 = A8 + S8*(1)
    A8 = A8 + S2p*(56/65*s6 + 56/143*s2*s4)
    A8 = A8 + S4p*(1960/1287*s4 + 168/143*s2**2)
    A8 = A8 + S6p*(168/65*s2)
    A8 = A8 + S8p*(1)
    A8 = A8 + m/3*(-56/65*s6 - 56/143*s2*s4)

    new_s = np.array((A2/S0, A4/S0, A6/S0, A8/S0)).T
    return new_s

def Upu(ss, SS, m):
    # Following Nettelmann 2017 eqs. B3 and B.4, assuming equipotential.
    s2 = ss[1]
    s4 = ss[2]
    S0 = SS[0]; S2 = SS[1]; S4 = SS[2]
    S0p = SS[5]; S2p = SS[6]; S4p = SS[7]
    A0 = np.zeros(s2.shape)
    A0 = A0 + S0*(1 + 2/5*s2**2 - 4/105*s2**3 + 2/9*s4**2 + 43/175*s2**4 - 4/35*s2**2*s4)
    A0 = A0 + S2*(-3/5*s2 + 12/35*s2**2 - 234/175*s2**3 + 24/35*s2*s4)
    A0 = A0 + S4*(-5/9*s4 + 6/7*s2**2)
    A0 = A0 + S0p*(1)
    A0 = A0 + S2p*(2/5*s2 + 2/35*s2**2 + 4/35*s2*s4 - 2/25*s2**3)
    A0 = A0 + S4p*(4/9*s4 + 12/35*s2**2)
    A0 = A0 + m/3*(1 - 2/5*s2 - 9/35*s2**2 - 4/35*s2*s4 + 22/525*s2**3)

    return -A0

def _test():
    from timeit import default_timer as timer
    N = 1024
    zvec = np.linspace(1, 1/N, N)
    dvec = -3000*zvec**2 + 3000
    mrot = 0.08
    nx = 128
    tic = timer()
    Js, out = tof4(zvec, dvec, mrot, xlevels=nx)
    toc = timer()
    print(f"With N = {N} and nx = {nx}.")
    print(f"Elapsed time {toc - tic:0.2} seconds.")
    print(f"After {out.it+1} iterations:")
    print(f"J0 = {Js[0]}")
    print(f"J2 = {Js[1]}")
    print(f"J4 = {Js[2]}")
    print(f"J6 = {Js[3]}")
    print(f"J8 = {Js[4]}")
    print(f"a0 = {out.a0}")
    print(f"q = {out.qrot}")
    print(f"I = {out.NMoI}")
    print("")

if __name__ == '__main__':
    _test()
