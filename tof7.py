#------------------------------------------------------------------------------
#  Seventh-order Theory of Figures gravity calculator
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#------------------------------------------------------------------------------
import sys, os
import json
import numpy as np
import warnings
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

c7file = os.path.join(os.path.dirname(__file__),'tof7_coeffs.json')
with open(c7file,'r') as f:
    C7 = json.load(f)

def tof7(zvec, dvec, mrot, **kwargs):
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
    tol : scalar, positive, (tol=1e-8)
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
    Js : 1-by-8 vector, real
        Even harmonic gravity coefficients J0 to J14 (J0 is included as a
        sanity check and test of convergence).
    out : struct
        A structure holding other quantities calculated in the process,
        including the shape functions that define the full solution.

    Algorithm
    ---------
    Theory of figures equations and coefficients from Nettelmann et al. (2021).
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
        ss = 8*[np.zeros(N)]
    else:
        ss = opts['ss_guesses']
    if np.isscalar(opts['xlevels']):
        sskip = int(max(np.fix(N/opts['xlevels']), 1))
        xind = range(0, N, sskip)
    else:
        raise(Exception("Vector xlevels feature not yet implemented"))

    # The loop, following Nettelmann (2017) Appendix B
    Js = np.array([0,0,0,0,0,0,0,0]) # J0=0 ensures at least one iteration
    for it in range(opts['maxiter']):
        # Equations B.16-B.17
        fs = skipnspline_B1617(ss, zvec, xind)

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
    return 1 + shp # TODO: improve precision by using higher order ss

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
    return cumtrapz(Y, X, initial=0)

def B111(ss, SS):
    # Return Js from SS, with Req/Rm a necessary bonus.
    N = len(ss[0]) - 1
    s0 = ss[0][N]; s2 = ss[1][N]; s4 = ss[2][N]; s6 = ss[3][N]; s8 = ss[4][N]
    s10 = ss[5][N]; s12 = ss[6][N]; s14 = ss[7][N]
    S0 = SS[0][N]; S2 = SS[1][N]; S4 = SS[2][N]; S6 = SS[3][N]; S8 = SS[4][N]
    S10 = SS[5][N]; S12 = SS[6][N]; S14 = SS[7][N]
    aos = (1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8 -
                    (63/256)*s10 + (231/1024)*s12 - (429/2048)*s14)
    J0  = -(aos**-0)*S0
    J2  = -(aos**-2)*S2
    J4  = -(aos**-4)*S4
    J6  = -(aos**-6)*S6
    J8  = -(aos**-8)*S8
    J10 = -(aos**-10)*S10
    J12 = -(aos**-12)*S12
    J14 = -(aos**-14)*S14
    Js  = np.array([J0, J2, J4, J6, J8, J10, J12, J14])
    return (Js, aos)

def B9(Z, D, fs):
    # The ToF7 equivalent of Nettelmann 2017 eq. B.9.
    f0 = fs[0]; f2 = fs[1]; f4 = fs[2]; f6 = fs[3]; f8 = fs[4]
    f10 = fs[5]; f12 = fs[6]; f14 = fs[7]
    f0p = fs[8]; f2p = fs[9]; f4p = fs[10]; f6p = fs[11]; f8p = fs[12]
    f10p = fs[13]; f12p = fs[14]; f14p = fs[15]
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

    I10 = mcumtrapz(D, Z**(10+3)*f10)
    S10 = D*f10 - Z**-(10+3)*I10

    I12 = mcumtrapz(D, Z**(12+3)*f12)
    S12 = D*f12 - Z**-(12+3)*I12

    I14 = mcumtrapz(D, Z**(14+3)*f14)
    S14 = D*f14 - Z**-(14+3)*I14

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

    I10p = mcumtrapz(D, Z**(2-10)*f10p)
    I10p = I10p[N] - I10p
    S10p = -D*f10p + Z**-(2-10)*(D[N]*f10p[N] - I10p)

    I12p = mcumtrapz(D, Z**(2-12)*f12p)
    I12p = I12p[N] - I12p
    S12p = -D*f12p + Z**-(2-12)*(D[N]*f12p[N] - I12p)

    I14p = mcumtrapz(D, Z**(2-14)*f14p)
    I14p = I14p[N] - I14p
    S14p = -D*f14p + Z**-(2-14)*(D[N]*f14p[N] - I14p)

    SS = [S0,S2,S4,S6,S8,S10,S12,S14,S0p,S2p,S4p,S6p,S8p,S10p,S12p,S14p]
    return SS

def skipnspline_B1617(ss, zvec, xind):
    """Update the ToF7 equivalent of B.16-B.17 for new f_2n."""

    # Skip
    zs = np.array([s[xind] for s in ss[1:]]).T
    fs = B1617(zs)

    # And spline
    if len(xind) < len(zvec):
        X = zvec[xind]
        for k in range(len(fs)):
            Y = fs[k]
            fs[k] = interp1d(X, Y, 'cubic',fill_value='extrapolate')(zvec)

    return fs

def B1617(sarray):
    """The ToF7 equivalent of Nettelmann 2017 eqs. B.16 and B.17."""

    N = sarray.shape[0]
    ones = np.ones((N,1))

    f0 = np.zeros((N,1))
    block = np.array(C7['f']['f0'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f0 = f0 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f0 = f0 + c*ones
    f0 = f0.squeeze()

    f2 = np.zeros((N,1))
    block = np.array(C7['f']['f2'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f2 = f2 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f2 = f2 + c*ones
    f2 = f2.squeeze()

    f4 = np.zeros((N,1))
    block = np.array(C7['f']['f4'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f4 = f4 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f4 = f4 + c*ones
    f4 = f4.squeeze()

    f6 = np.zeros((N,1))
    block = np.array(C7['f']['f6'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f6 = f6 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f6 = f6 + c*ones
    f6 = f6.squeeze()

    f8 = np.zeros((N,1))
    block = np.array(C7['f']['f8'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f8 = f8 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f8 = f8 + c*ones
    f8 = f8.squeeze()

    f10 = np.zeros((N,1))
    block = np.array(C7['f']['f10'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f10 = f10 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f10 = f10 + c*ones
    f10 = f10.squeeze()

    f12 = np.zeros((N,1))
    block = np.array(C7['f']['f12'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f12 = f12 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f12 = f12 + c*ones
    f12 = f12.squeeze()

    f14 = np.zeros((N,1))
    block = np.array(C7['f']['f14'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f14 = f14 + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f14 = f14 + c*ones
    f14 = f14.squeeze()

    f0p = np.zeros((N,1))
    block = np.array(C7['f']['f0p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f0p = f0p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f0p = f0p + c*ones
    f0p = f0p.squeeze()

    f2p = np.zeros((N,1))
    block = np.array(C7['f']['f2p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f2p= f2p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f2p = f2p + c*ones
    f2p = f2p.squeeze()

    f4p = np.zeros((N,1))
    block = np.array(C7['f']['f4p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f4p = f4p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f4p = f4p + c*ones
    f4p = f4p.squeeze()

    f6p = np.zeros((N,1))
    block = np.array(C7['f']['f6p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f6p = f6p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f6p = f6p + c*ones
    f6p = f6p.squeeze()

    f8p = np.zeros((N,1))
    block = np.array(C7['f']['f8p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f8p = f8p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f8p = f8p + c*ones
    f8p = f8p.squeeze()

    f10p = np.zeros((N,1))
    block = np.array(C7['f']['f10p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f10p = f10p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f10p = f10p + c*ones
    f10p = f10p.squeeze()

    f12p = np.zeros((N,1))
    block = np.array(C7['f']['f12p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f12p = f12p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f12p = f12p + c*ones
    f12p = f12p.squeeze()

    f14p = np.zeros((N,1))
    block = np.array(C7['f']['f14p'])
    for k in range(block.shape[0]):
        c = block[k,-1]
        pows = block[k,:-1]
        if np.any(pows):
            f14p = f14p + c*np.prod(sarray**pows, axis=1).reshape(N,1)
        else:
            f14p = f14p + c*ones
    f14p = f14p.squeeze()

    fs = [f0,f2,f4,f6,f8,f10,f12,f14,f0p,f2p,f4p,f6p,f8p,f10p,f12p,f14p]
    return fs

def skipnspline_B1215(ss0, SS, mrot, zvec, xind):
    # Update the ToF7 equivalent of the system B.12-B.15 for new s_2n.

    # Skip
    Zs = np.array([S[xind] for S in SS]).T
    zs = np.array([s[xind] for s in ss0]).T
    newzs = B1215(zs, Zs, mrot)
    newz0 = A24(newzs)
    newz0 = np.reshape(newz0, (newz0.size,1))
    Y = np.hstack((newz0,newzs))

    # And spline
    if len(xind) < len(zvec):
        X = zvec[xind]
        s0 =  interp1d(X, Y[:,0], 'cubic', fill_value='extrapolate')(zvec)
        s2 =  interp1d(X, Y[:,1], 'cubic', fill_value='extrapolate')(zvec)
        s4 =  interp1d(X, Y[:,2], 'cubic', fill_value='extrapolate')(zvec)
        s6 =  interp1d(X, Y[:,3], 'cubic', fill_value='extrapolate')(zvec)
        s8 =  interp1d(X, Y[:,4], 'cubic', fill_value='extrapolate')(zvec)
        s10 = interp1d(X, Y[:,5], 'cubic', fill_value='extrapolate')(zvec)
        s12 = interp1d(X, Y[:,6], 'cubic', fill_value='extrapolate')(zvec)
        s14 = interp1d(X, Y[:,7], 'cubic', fill_value='extrapolate')(zvec)
    else:
        s0  = Y[:,0]
        s2  = Y[:,1]
        s4  = Y[:,2]
        s6  = Y[:,3]
        s8  = Y[:,4]
        s10 = Y[:,5]
        s12 = Y[:,6]
        s14 = Y[:,7]
    ss = [s0, s2, s4, s6, s8, s10, s12, s14]
    return ss

def B1215(s, S, m):
    # Compute ToF7 equivalent of RHS of B.12-B.15 and "solve" for sn.

    s = s[:,1:] # s0 is not part of this system
    N = s.shape[0]
    mterm = (m/3)*np.ones((N,1))
    terms = np.hstack((S,mterm))
    fields = list(C7['A']['A2'].keys()) # same for all As
    nt = len(fields)
    assert terms.shape[1] == nt

    # A2 (not including -s2S0)
    A2 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A2'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A2 = A2 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A2 = A2 + c*terms[:,j]

    # A4 (not including -s4S0)
    A4 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A4'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A4 = A4 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A4 = A4 + c*terms[:,j]

    # A6 (not including -s6S0)
    A6 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A6'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A6 = A6 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A6 = A6 + c*terms[:,j]

    # A8 (not including -s8S0)
    A8 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A8'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A8 = A8 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A8 = A8 + c*terms[:,j]

    # A10 (not including -s10S0)
    A10 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A10'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A10 = A10 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A10 = A10 + c*terms[:,j]

    # A12 (not including -s12S0)
    A12 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A12'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A12 = A12 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A12 = A12 + c*terms[:,j]

    # A14 (not including -s14S0)
    A14 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A14'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A14 = A14 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A14 = A14 + c*terms[:,j]

    S0 = S[:,0]
    new_s = np.array((A2/S0,A4/S0,A6/S0,A8/S0,A10/S0,A12/S0,A14/S0)).T
    return new_s

def B4(ss, SS, m):
    # Compute ToF7 equivalent of B.4.
    S = np.array([z for z in SS]).T
    s = np.array([z for z in ss]).T
    s = s[:,1:] # s0 is not part of this system
    N = s.shape[0]
    mterm = (m/3)*np.ones((N,1))
    terms = np.hstack((S,mterm))
    fields = list(C7['A']['A0'].keys())
    nt = len(fields)
    assert terms.shape[1] == nt

    A0 = np.zeros(N)
    for j in range(nt):
        block = np.array(C7['A']['A0'][fields[j]]).reshape(-1,8)
        for row in block:
            c = row[-1]
            pows = row[:-1]
            if np.any(pows):
                A0 = A0 + c*terms[:,j]*np.prod(s**pows,axis=1)
            else:
                A0 = A0 + c*terms[:,j]

    return A0

def A24(sn):
    """Nettelmann et al 2021 eq. A.24."""
    s2 = sn[:,0]; s4 = sn[:,1]; s6 = sn[:,2];

    s0 = 0*s2
    s0 = s0 - (1/5)*s2**2
    s0 = s0 - (2/105)*s2**3
    s0 = s0 - (1/9)*s4**2 - (2/35)*s2**2*s4
    s0 = s0 - (2/525)*s2**5 - (20/693)*s2*s4**2
    s0 = s0 + (127/55125)*s2**6 - (2/175)*s2**4*s4 - (6/1001)*s4**3 - \
              (10/143)*s2*s4*s6 - (1/13)*s6**2
    s0 = s0 - (2/2625)*s2**7 - (82/10395)*s2**3*s4**2 - (8/525)*s2**5*s4 - \
              (14/715)*s2*s6**2 - (20/1257)*s4**2*s6

    return s0

def _test():
    from timeit import default_timer as timer
    N = 1024
    zvec = np.linspace(1, 1/N, N)
    dvec = -3000*zvec**2 + 3000
    mrot = 0.08
    nx = 128
    tic = timer()
    Js, out = tof7(zvec, dvec, mrot, xlevels=nx)
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
