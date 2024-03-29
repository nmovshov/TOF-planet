#------------------------------------------------------------------------------
# Interior model of rotating fluid planet.
#------------------------------------------------------------------------------
import numpy as np
import tof4
import tof7
import warnings
from timeit import default_timer as timer

class TPError(Exception):
    pass

class TOFPlanet:
    """Interior model of rotating fluid planet.

    This class implements a model of a rotating fluid planet using Theory of
    Figures to calculate the hydrostatic equilibrium shape and resulting
    gravity field. A TOFPlanet object is defined by a density profile rho(s),
    supplied by the user and stored in the column vectors si and rhoi,
    indexed from the surface in. To complete the definition the user must also
    specify a mass, equatorial radius, and rotation period. With these a
    gravity field and equilibrium shape can be determined, with a call to
    relax_to_HE().

    Note that the oblate shape calculated with relax_to_HE() preserves the mass
    of the planet but not the equatorial radius. If fixradius is True (default:
    True) the si will be re-normalized to match the reference equatorial
    radius, modifying the implied mass. If fixmass is True (default: True) the
    density is then re-normalized to match the reference mass, modifying the
    reference 1-bar density. It is not possible to define mass, radius, and
    density simultaneously.

    A call to relax_to_HE() with fixradius True and/or fixmass False modifies
    the implied rotation period, given a fixed rotation parameter m. To
    preserve the reference mass, radius, and rotation period, it is necessary
    to run the equilibrium shape calculation iteratively, recalculating the
    rotation parameter between iterations. This is done with a call to
    relax_to_rotation(). Normally, a call to relax_to_rotation() with default
    flags is the recommended way to obtain a self-consistent model planet
    suitable for direct comparison with observation and/or third-party models.

    Alternatively the user may supply a barotrope, rho(P,r) function and call
    relax_to_barotrope() to iteratively find a density profile consistent with
    the calculated equilibrium pressure. To this end a boundary pressure must
    also be given and an initial density profile guess is still required (can
    be a simple one). Again, we can't simultaneously impose an exact mass,
    radius, and barotrope. By default the reference mass and radius will be
    honored, by re-normalizing the converged density profile, and this will
    modify the _effective_ barotrope.
    """
    def __init__(self, obs=None, **kwargs):
        self.G = 6.67430e-11 # m^3 kg^-1 s^-2 (2018 NIST reference)
        self.opts = _default_opts(kwargs) # holds user configurable options
        if obs:
            self.set_observables(obs)
        self.M = None
        self.NMoI = None
        self.Pi = None
        self.ss = None
        self.SS = None

    def set_observables(self,obs):
        """Copy physical properties from an observables struct."""
        self.name = obs.pname
        self.mass = obs.M
        self.radius = obs.a0
        self.P0 = obs.P0
        self.period = obs.P
        self.GM = self.G*self.mass
        self.rhobar = self.mass/(4*np.pi/3*obs.s0**3)

    def set_barotrope(self, fun):
        """Store a function with signature density = f(pressure, radius)."""
        self.baro = fun

    def relax_to_rotation(self, fixmass=True):
        """Call relax_to_He repeatedly for simultaneous shape and rotation."""

        self.Js = np.hstack((-1, np.zeros(self.opts['toforder'])))
        self.opts['MaxIterHE'] = 3 # optimal from token benchmark
        self.opts['verbosity'] = 0 # silence HE warnings

        self.wrot = 2*np.pi/self.period
        self.mrot = self.wrot**2*self.si[0]**3/self.GM
        it = 1
        while (it < self.opts['MaxIterRot']) and (self.mrot > 0):
            it = it + 1
            old_Js = self.Js
            old_m = self.mrot
            self.relax_to_HE(fixradius=True,fixmass=fixmass,fixrot=True)
            dJs = np.abs((self.Js[1] - old_Js[1]))
            drot = np.abs(old_m - self.mrot)
            if (drot < self.opts['drottol']) and (dJs < self.opts['dJtol']):
                break
        if it == self.opts['MaxIterRot']:
            warnings.warn('Rotation period may not be fully converged.')
        return it

    def relax_to_barotrope(self, fixmass=True):
        """Call relax_to_He repeatedly to converge shape/rotation/density."""

        self.Js = np.hstack((-1, np.zeros(self.opts['toforder'])))
        self.opts['MaxIterHE'] = 2
        self.opts['verbosity'] = 0 # silence HE warnings

        self.wrot = 2*np.pi/self.period
        self.mrot = self.wrot**2*self.si[0]**3/self.GM

        it = 1
        while (it < self.opts['MaxIterBar']):
            it = it + 1
            old_Js = self.Js
            old_m = self.mrot
            old_ro = self.rhoi

            self.relax_to_HE(fixradius=True,fixmass=fixmass,fixrot=True)
            self.update_densities(renorm=True)

            dJs = np.abs((self.Js[1] - old_Js[1]))
            drot = np.abs(old_m - self.mrot)
            dro = np.var(self.rhoi[1:]/old_ro[1:])
            if (drot < self.opts['drottol'] and
                dJs < self.opts['dJtol'] and
                dro < self.opts['drhotol']):
                break
        if it == self.opts['MaxIterBar']:
            warnings.warn('Barotrope may not be fully converged.')
        return it

    def relax_to_HE(self, fixradius=True, fixmass=False, fixrot=False,
                    moi=False, pressure=False):
        """Call tof<n> once to obtain equilibrium shape and gravity."""

        if (self.opts['verbosity'] > 1):
            print('  Relaxing to hydrostatic equilibrium...')
        tic = timer()
        if self.opts['toforder'] == 4:
            tofun = tof4.tof4
        elif self.opts['toforder'] == 7:
            tofun = tof7.tof7
        else:
            raise ValueError('Unimplemented tof order')

        self.Js, out = tofun(self.si, self.rhoi, self.mrot,
            tol=self.opts['dJtol'], maxiter=self.opts['MaxIterHE'],
            xlevels=self.opts['xlevels'], verbosity=self.opts['verbosity'],
            ss_guesses=self.ss,
            calc_moi=moi)
        toc = timer() - tic

        if (self.opts['verbosity'] > 1):
            print('  Relaxing to hydrostatic equilibrium...done.')
            print(f' Elapsed time {toc:g} sec.')

        self.NMoI = out.NMoI
        self.A0 = np.flip(out.A0)
        self.ss = out.ss # consider flipping used other than ss_guesses
        self.SS = out.SS # consider flipping if used anywhere
        self.aos = out.a0

        if fixradius:
            self._fixradius()
        if fixmass:
            self._fixmass()
        if fixrot:
            self._fixrot()
        if pressure:
            self._pressurize()

    def update_densities(self,renorm=True):
        # Call baro to set rhoi from Pi
        self._pressurize()
        self.rhoi = self.baro(self.Pi)
        if np.any(np.diff(self.rhoi) < 0):
            raise TPError("barotrope created density inversion!")
        if renorm:
            self._fixradius()
            self._fixmass()
            self._fixrot()

    def level_surfaces(self,mus):
        # Normalized r(cos(theta)) from shape functions
        ss = self.ss
        s0 = ss[0]; s2 = ss[1]; s4 = ss[2]; s6 = ss[3]; s8 = ss[4]
        shp = (np.outer(s0,_Pn(0,mus)) + np.outer(s2,_Pn(2,mus)) +
               np.outer(s4,_Pn(4,mus)) + np.outer(s6,_Pn(6,mus)) +
               np.outer(s8,_Pn(8,mus)))
        return np.flipud(1 + shp)

    ### Private methods slash pseudo properties
    def _fixradius(self):
        self.si = self.si*self.radius/(self.si[0]*self.aos)
        self.s0 = self.si[0]
        self.a0 = self.s0*self.aos
        self.M = _mass_int(self.si, self.rhoi)

    def _fixmass(self):
        self.rhoi = self.rhoi*self.mass/self.M
        self.M = _mass_int(self.si, self.rhoi)

    def _fixrot(self):
        self.mrot = self.wrot**2*self.s0**3/(self.G*self.M)

    def _pressurize(self):
        r = self.si
        rho = self.rhoi
        U = -self.G*self.mass/self.s0**3*self.si**2*self.A0
        gradU = np.zeros_like(r)
        gradU[0] = (U[0] - U[1])/(r[0] - r[1])
        gradU[1:-1] = (U[0:-2] - U[2:])/(r[0:-2] - r[2:])
        gradU[-1] = (U[-2] - U[-1])/(r[-2] - r[-1])
        intgrnd = rho*gradU
        P = np.zeros_like(r)
        P[0] = self.P0
        for k in range(P.size-1): # note integrate downward
            P[k+1] = P[k] + 0.5*(r[k] - r[k+1])*(intgrnd[k] + intgrnd[k+1])
        self.Pi = P

    def _m2P(self):
        return 2*np.pi/(np.sqrt(self.mrot*(self.GM)/self.s0**3))

    def _P2m(self):
        return (2*np.pi/self.period)**2*self.s0**3/(self.GM)

    def _P2q(self):
        return (2*np.pi/self.period)**2*self.a0**3/(self.GM)

    def ai(self):
        return self.s0*self.level_surfaces(0) # level surfaces equatorial radii

    def bi(self):
        return self.s0*self.level_surfaces(1) # level surfaces polar radii

    ### Visualizers
    def plot_rho_of_r(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = np.append(self.si/self.si[0], 0)
        y = np.append(self.rhoi, self.rhoi[-1])/1000
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'Level surface radius, $s/s_0$', fontsize=12)
        plt.ylabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.show(block=False)

    def plot_P_of_r(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.si/self.si[0]
        y = self.Pi/1e9
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'Level surface radius, $s/s_0$', fontsize=12)
        plt.ylabel(r'$P$ [GPa]', fontsize=12)
        plt.show(block=False)

    def plot_P_of_rho(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.rhoi/1000
        y = self.Pi/1e9
        plt.plot(x, y, lw=2, label=self.name)
        plt.xlabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.ylabel(r'$P$ [GPa]', fontsize=12)
        plt.show(block=False)

    def plot_rho_of_P(self):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,6))
        x = self.Pi/1e9
        y = self.rhoi/1000
        plt.loglog(x, y, lw=2, label=self.name)
        plt.xlabel(r'$P$ [GPa]', fontsize=12)
        plt.ylabel(r'$\rho$ [1000 kg/m$^3$]', fontsize=12)
        plt.xlim(left=1e-3)
        plt.show(block=False)

    ## Reporters/exporters
    def to_ascii(self, fname=None):
        fname = fname if fname else (self.name+'.txt')
        with open(fname, 'wt') as fid:
            # The header
            fid.write('# Equilibrium shape and gravity solved with ')
            fid.write(f"{self.opts['toforder']}th-order Theory of Figures\n")
            fid.write('#\n')
            fid.write(f'# Model name: {self.name}\n')
            fid.write('#\n')
            fid.write(f'# Scalar quantities:\n')
            fid.write(f'# N levels = {len(self.si)}\n')
            fid.write(f'# xlevels = {self.opts["xlevels"]}\n')
            fid.write(f'# Mass M = {self.M} kg\n')
            fid.write(f'# Mean radius       s0 = {self.s0:0.6e} m\n')
            fid.write(f'# Equatorial radius a0 = {self.a0:0.6e} m\n')
            fid.write(f'# Rotation period P = {self.period:0.6g} s\n')
            fid.write(f'# Normalized MOI = {self.NMoI:0.6f}\n')
            fid.write(f'#\n')
            fid.write(f'# Calculated gravity zonal harmonics (x 10^6):\n')
            fid.write(f'# J0  = {self.Js[0]*1e6:12.6f}\n')
            fid.write(f'# J2  = {self.Js[1]*1e6:12.6f}\n')
            fid.write(f'# J4  = {self.Js[2]*1e6:12.6f}\n')
            fid.write(f'# J6  = {self.Js[3]*1e6:12.6f}\n')
            fid.write(f'# J8  = {self.Js[4]*1e6:12.6f}\n')
            fid.write(f'#\n')
            fid.write('# Column data description (MKS):\n')
            fid.write('# i     - level surface index (increasing with depth)\n')
            fid.write('# s_i   - mean radius of level surface i\n')
            fid.write('# rho_i - density on level surfaces i\n')
            fid.write('# P_i   - pressure on level surface i\n')
            fid.write('#\n')

            # The data
            fid.write('# Column data:\n')
            fid.write(f'# {"i":>4}')
            fid.write(f'{"s_i":>8}')
            fid.write(f'{"rho_i":>14s}')
            fid.write(f'{"P_i":>10}')
            fid.write('\n')
            for k in range(len(self.si)):
                fid.write(f'  {k:-4d}  ')
                fid.write(f'{self.si[k]:10.4e}  ')
                fid.write(f'{self.rhoi[k]:10.4e}  ')
                fid.write(f'{self.Pi[k]:10.4e}  ')
                fid.write(f'\n')
        return

# Class-related functions
def _Pn(n, x):
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

def _mass_int(svec, dvec):
    """Trapz-integrate mass from rho(r) data."""
    from scipy.integrate import trapz
    return -4*np.pi*trapz(dvec*svec**2, x=svec)

def _default_opts(kwargs):
    """Return options dict used by TOFPlanet class methods."""
    opts = {'toforder':4,
            'dJtol':1e-10,
            'drottol':1e-6,
            'drhotol':1e-6,
            'MaxIterHE':60,
            'MaxIterBar':60,
            'MaxIterRot':20,
            'xlevels':-1,
            'verbosity':1
            }
    for kw, v in kwargs.items():
        if kw in opts:
            opts[kw] = v
    return opts

class _default_planet:
    """Use this to prefill critical TP fields with reasonable values."""
    pname = 'planet'
    M  = 1898.187e24
    a0 = 71492e3
    s0 = 69911e3
    P0 = 1e5
    P = 0.41354*24*3600

def _a_jupi(N):
    tp = TOFPlanet(obs=_default_planet)
    a = -15*tp.mass/8/np.pi/tp.radius**3
    zvec = np.linspace(1, 1/N, N)
    dvec = a*zvec**2 - a
    tp.si = zvec*tp.radius
    tp.rhoi = dvec
    return tp

def _test_rot(N,nx,torder):
    tp = TOFPlanet(obs=_default_planet,xlevels=nx)
    zvec = np.linspace(1, 1/N, N)
    dvec = -3000*zvec**2 + 3000
    tp.si = zvec*tp.s0
    tp.rhoi = dvec
    tp.opts['toforder'] = torder
    tic = timer()
    it = tp.relax_to_rotation()
    toc = timer()
    print()
    print(f"{N=}, {nx=}")
    print(f"{it} iterations of tof{torder} in {toc-tic:.3g} sec.")
    print("J0 = {}".format(tp.Js[0]))
    print("J2 = {}".format(tp.Js[1]))
    print("J4 = {}".format(tp.Js[2]))
    print("J6 = {}".format(tp.Js[3]))
    print("J8 = {}".format(tp.Js[4]))
    print("I = {}".format(tp.NMoI))
    print("")

def _test_baro(N,nx,torder):
    tp = TOFPlanet(obs=_default_planet,xlevels=nx)
    zvec = np.linspace(1, 1/N, N)
    dvec = -3000*zvec**2 + 3000
    tp.si = zvec*tp.radius
    tp.rhoi = dvec
    tp.opts['toforder'] = torder
    def poly1(P): return np.sqrt(P/2e5)
    tp.set_barotrope(poly1)

    tic = timer()
    it = tp.relax_to_barotrope()
    toc = timer()
    print()
    print(f"{N=}, {nx=}")
    print(f"{it} iterations of tof{torder} in {toc-tic:.3g} sec.")
    print("J0 = {}".format(tp.Js[0]))
    print("J2 = {}".format(tp.Js[1]))
    print("J4 = {}".format(tp.Js[2]))
    print("J6 = {}".format(tp.Js[3]))
    print("J8 = {}".format(tp.Js[4]))
    print("I = {}".format(tp.NMoI))
    print("")

if __name__ == '__main__':
    _test_baro(4096,-1,4)
