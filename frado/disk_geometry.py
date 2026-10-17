"""Reproduce the original disk_shape.m bisection-midpoint interpolation (not linear)."""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import PchipInterpolator
from .config import Config

# (log10 M, 100 * Eddington ratio): (a, b) for H=a*R**b
FITS = {
 (6,1):(.001294,1.154),(6,10):(.001826,1.15),(6,100):(.0132,1.),
 (7,1):(.001977,1.095),(7,10):(.0008029,1.207),(7,100):(.002524,1.128),
 (8,1):(.002615,1.04),(8,10):(.001016,1.169),(8,100):(.002472,1.119),
 (9,1):(.003223,1.),(9,10):(.001942,1.091),(9,100):(.00141,1.154),
 (10,1):(.00241,1.),(10,10):(.0003,1.215),(10,100):(.00196,1.11),
}


def _midpoint_interpolate(lo, hi, target, lo_value, hi_value, tol=1e-3):
    # MATLAB computes one midpoint even at exactly matching grid values.
    mid_value = (lo_value + hi_value)/2
    while (hi-lo)/2 > tol:
        mid = (hi+lo)/2
        mid_value = (hi_value+lo_value)/2
        if mid < target:
            lo, lo_value = mid, mid_value
        else:
            hi, hi_value = mid, mid_value
    return mid_value

@dataclass
class DiskGeometry:
    radius: np.ndarray
    height: np.ndarray
    interpolator: PchipInterpolator

    def at(self, r): return self.interpolator(r)


def build_disk(cfg: Config) -> DiskGeometry:
    if not (6 <= cfg.bh_mass_log <= 10 and -2 <= cfg.mdot_log <= 0):
        raise ValueError('Disk fit is defined only for 6<=logM<=10, -2<=logmdot<=0')
    radius = np.linspace(100, 100000, 20000)
    ml, mh = int(np.floor(cfg.bh_mass_log)), int(np.ceil(cfg.bh_mass_log))
    dl, dh = int(np.floor(cfg.mdot_log)), int(np.ceil(cfg.mdot_log))
    def curve(m, d):
        key = (m, int(round(100*10.0**d)))
        a, b = FITS[key]
        return a * radius**b
    at_dh = _midpoint_interpolate(ml, mh, cfg.bh_mass_log, curve(ml,dh), curve(mh,dh))
    at_dl = _midpoint_interpolate(ml, mh, cfg.bh_mass_log, curve(ml,dl), curve(mh,dl))
    h = _midpoint_interpolate(dl, dh, cfg.mdot_log, at_dl, at_dh)
    return DiskGeometry(radius, h, PchipInterpolator(radius,h,extrapolate=False))
