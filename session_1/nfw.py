# Last-modified: 16 Jul 2014 15:11:18
import numpy as np
from scipy import special
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d

class NFW(object):

    """ Define an NFW halo with mass, concentration, background density and
    overdensity factor. All distance scales are physical rather than comoving.

    Parameters
    ----------
    mass: float, scalar
        Halo virial mass.
    c: float, scalar
        Halo concentration.
    rho_mean: float, scalar
        Background density (or a.k.a mean density,
        :math:`\Omega_m \\times \\rho_{crit} `)
    DELTA_HALO: scalar, optional
        Overdensity factor, default to be 200, which defines halos
        as having an average density of 200*rho_mean.

    """

    def __init__(self, mass, c, rho_mean, DELTA_HALO=200.0):
        self.mass = mass
        self.c = c
        self.rho_vir = rho_mean*DELTA_HALO
        self.r_vir = self._r_vir()
        self.r_s = self.r_vir/self.c
        self.rho_0 = self.mass/self._mass()

    def __str__(self):
        return(("M_h: %g C: %g R_vir: %g R_s: %g Mass within 2*R_vir: %g") % (
            self.mass, self.c, self.r_vir, self.r_s,
            self.total_mass(2.*self.r_vir)))

    def profile(self, r):
        """ Returns the normalized NFW profile.
        """
        return(self._profile(r, rho_0=self.rho_0))

    def total_mass(self, r_max):
        """ Total mass inside of some radius :math:`r_{max}`.

        Parameters
        ----------
        r_max: float, scalar
            Maximum radius to which NFW profile is to be integrated to
            get the enclosed mass.

        Returns
        -------
        mass: float, scalar
            Enclosed mass.

        """
        return(4.0*np.pi*self.rho_0*self.r_s**3*(
            np.log((self.r_s + r_max)/self.r_s) - r_max/(self.r_s + r_max)))

    def DelSig(self, r):
        """ Following the analytic formalism in Wright&Brainerd00.
        """
        _r = np.atleast_1d(r)
        x = _r/self.r_s
        _g = self._glt_or_ggt(x)
        return(self.r_s * self.rho_0 * _g)

    def _mass(self, rho_0=1.0):
        """ Internal function, returns the virial mass normalized by rho_0.
        """
        return(4.0*np.pi*rho_0*self.r_s**3*(np.log(1.+self.c)-self.c/(
            1.+self.c)))

    def _r_vir(self):
        """ Returns the virial radius.
        """
        # given mean density, we can solve for the virial radius from input
        # mass.
        return(np.power((3.0/(4.0*np.pi*self.rho_vir)) * self.mass, 1./3.))

    def _profile(self, r, rho_0=1.0):
        """ Internal function, returns the NFW profile normalized by rho_0.
        """
        _reff = r/self.r_s
        return(rho_0/(_reff*(1.0+_reff)**2))

    def _glt_or_ggt(self, x):
        """ Eqn. 15 and 16 in Wright & Brainerd 2000

        x is the scaled radius r/r_s, and is expected to be an ndarray

        """
        x2 = x*x
        pub = (4./x2) * np.log(x/2.0)
        _g = np.empty_like(x)
        at = np.empty_like(x)
        ilt = x < 1.0
        igt = x > 1.0
        ieq = x == 1.0
        ineq = ilt | igt
        if np.any(ilt):
            at[ilt] = np.arctanh(
                np.sqrt((1.0-x[ilt])/(1.0+x[ilt]))) / np.sqrt(1.-x2[ilt])
        if np.any(igt):
            at[igt] = np.arctan(
                np.sqrt((x[igt]-1.0)/(x[igt]+1.0))) / np.sqrt(x2[igt] - 1.)
        if np.any(ieq):
            _g[ieq] = pub[ieq]+10./3.
        _g[ineq] = 4.0*at[ineq] * \
            (2./x2[ineq] + 1./(x2[ineq]-1.0)) + pub[ineq] - 2.0/(x2[ineq]-1.0)
        return(_g)
