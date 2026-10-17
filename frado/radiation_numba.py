"""Opt-in QUADPACK backend with a compiled, low-level quadrature integrand.

The integration *method* (nested adaptive Gauss-Kronrod/QUADPACK) and
physical expressions are unchanged from the SciPy/Python reference backend.
Numba compiles only the scalar callback, avoiding Python work per sample.
The buffer and callback are local to each integration for future thread safety.
"""
import ctypes
import math
import numpy as np
from scipy import LowLevelCallable
from scipy.integrate import nquad
from numba import cfunc, carray, types
from .constants import TOTAL_FIT, ABS_FIT

# QUADPACK nquad integrates x[0]=phi first, x[1]=R second.
# Parameters = [rho, r_squared, rmin, opacity_offset, kind, flattened fit(21)].
_SIGNATURE = types.double(types.intc, types.CPointer(types.double), types.voidptr)


@cfunc(_SIGNATURE)
def _radiation_integrand(n, xx, user_data):
    data = carray(user_data, (26,), dtype=types.float64)
    phi, R = xx[0], xx[1]
    rho, r_squared, rmin, offset, kind = data[0], data[1], data[2], data[3], data[4]
    cphi = math.cos(phi)
    frac = 1.0-math.sqrt(rmin/R)
    k = offset + 0.25*math.log10(frac)-0.75*math.log10(R)
    fit_sum = 0.0
    for i in range(7):
        j = 5+3*i
        fit_sum += data[j]*math.sin(data[j+1]*k+data[j+2])
    opacity = 10.0 ** fit_sum
    d2 = r_squared+R*R-2.0*R*rho*cphi
    if d2 <= 0.0:
        return math.nan
    base = opacity*R*frac/(R**3)
    if kind == 2.0:  # heating / absorption
        return base / (d2**1.5)
    base /= d2*d2
    if kind == 0.0:  # cylindrical component
        return base * (rho-R*cphi)
    return base       # vertical component


def evaluate_integrals(cfg, scales, bounds, rho, z, kinds=(0, 1, 2)):
    """Return selected dimensionless quadratures WITHOUT physical prefactors."""
    rin, rout, phiin, phiout = bounds
    r_squared = rho*rho+z*z
    params = np.empty(26, dtype=np.float64)
    params[:4] = (rho, r_squared, cfg.inner_radiating_radius,
                  scales.opacity_argument_offset)
    pointer = params.ctypes.data_as(ctypes.c_void_p)
    callback = LowLevelCallable(_radiation_integrand.ctypes, pointer)
    options = {'epsabs': cfg.quad_epsabs, 'epsrel': cfg.quad_epsrel}
    values = []
    for kind in kinds:
        params[4] = kind
        params[5:] = (ABS_FIT if kind == 2 else TOTAL_FIT).ravel()
        value, error = nquad(callback, [[phiin, phiout], [rin, rout]], opts=options)
        if not math.isfinite(value):
            raise ValueError('Non-finite radiation integral with compiled backend')
        values.append(value)
    return tuple(values)


def evaluate_absorption(cfg, scales, bounds, rho, z):
    """Heating-only integral for prelaunch screening; no acceleration work."""
    return evaluate_integrals(cfg, scales, bounds, rho, z, kinds=(2,))[0]
