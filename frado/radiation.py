"""Reference quadrature, structurally equivalent to MATLAB integral2.

SciPy dblquad integrates phi first, then R; integral2 may use a different
adaptive algorithm, so strict numerical identity is not guaranteed.
"""
import numpy as np
from scipy.integrate import dblquad
from .constants import TOTAL_FIT, ABS_FIT

class Radiation:
    def __init__(self, cfg, scales, shielding):
        self.cfg, self.scales, self.shielding = cfg, scales, shielding
        self.backend = cfg.radiation_backend
        if self.backend == 'numba':
            try:
                from .radiation_numba import evaluate_integrals, evaluate_absorption
            except ImportError as exc:
                raise ImportError(
                    "radiation_backend='numba' requires Numba; "
                    "install with 'pip install numba', or select 'scipy'."
                ) from exc
            self._compiled_integrals = evaluate_integrals
            self._compiled_absorption = evaluate_absorption

    def _opacity(self, R, fit):
        k = self.scales.opacity_argument_offset + .25*np.log10(1-np.sqrt(self.cfg.inner_radiating_radius/R))-.75*np.log10(R)
        return 10**np.sum(fit[:,0]*np.sin(fit[:,1]*k+fit[:,2]))

    def absorbed_power(self, state):
        """Calculate dust heating at a state (without computing radiation force).

        The reference SciPy route deliberately calls the unchanged evaluate().
        The optimized route integrates the single absorption kernel directly.
        """
        if self.backend == 'scipy':
            return self.evaluate(state)[1]
        x, _, y, _, z, _ = state
        rho = np.hypot(x, y)
        bounds = self.shielding.bounds(rho, z)
        if bounds[0] >= bounds[1] or bounds[2] >= bounds[3]:
            return 0.0
        q_abs = self._compiled_absorption(self.cfg, self.scales, bounds, rho, z)
        return self.scales.absorption_coefficient*abs(z)*q_abs

    def evaluate(self, state):
        x, _, y, _, z, _ = state
        rho = np.hypot(x,y)
        r = np.sqrt(rho*rho+z*z)
        rin, rout, phiin, phiout = self.shielding.bounds(rho,z)
        if rin >= rout or phiin >= phiout:
            return np.zeros(3), 0.0
        if self.backend == 'numba':
            radial, vertical, q_abs = self._compiled_integrals(
                self.cfg, self.scales, (rin,rout,phiin,phiout), rho, z)
            absorbed = self.scales.absorption_coefficient*abs(z)*q_abs
            c = self.scales.radiative_coefficient*abs(z)
            ax = c*x/rho*radial if rho > 0 else 0.0
            ay = c*y/rho*radial if rho > 0 else 0.0
            az = c*z*vertical
            return np.array([ax,ay,az]), absorbed
        rmin = self.cfg.inner_radiating_radius
        def kernel(phi,R,kind):
            dist = r*r+R*R-2*R*rho*np.cos(phi)
            if dist <= 0: raise ValueError('Singularity in radiation integral')
            base = R*(1-np.sqrt(rmin/R))/R**3
            if kind == 'abs':
                return self._opacity(R, ABS_FIT)*base/dist**1.5
            base *= self._opacity(R, TOTAL_FIT)/dist**2
            if kind == 'rho': return base*(rho-R*np.cos(phi))
            return base
        def integrate(kind):
            value, _ = dblquad(lambda phi,R: kernel(phi,R,kind), rin,rout,
                lambda R: phiin, lambda R: phiout,
                epsabs=self.cfg.quad_epsabs, epsrel=self.cfg.quad_epsrel)
            return value
        radial = integrate('rho')
        vertical = integrate('z')
        absorbed = self.scales.absorption_coefficient*abs(z)*integrate('abs')
        c = self.scales.radiative_coefficient*abs(z)
        ax = c*x/rho*radial if rho > 0 else 0.0
        ay = c*y/rho*radial if rho > 0 else 0.0
        az = c*z*vertical
        return np.array([ax,ay,az]), absorbed
