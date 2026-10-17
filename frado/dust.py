"""Original wavelength convention and signed trapezoidal integration preserved."""
from pathlib import Path
import numpy as np
from .config import Config
from .constants import Scales

DEFAULT_OPACITY = Path(__file__).with_name('k_average_tot.dat')

def emission_power(cfg: Config, scales: Scales, opacity_path=DEFAULT_OPACITY):
    table = np.loadtxt(opacity_path)
    if table.ndim != 2 or table.shape[0] != 374 or table.shape[1] < 4:
        raise ValueError('Expected 374 opacity rows and at least four columns')
    wavelength = table[:,1]*1e-6
    absorption = table[:,3]*0.1
    with np.errstate(over="ignore"):
        integrand = absorption/wavelength**5 / np.expm1(scales.planck_coefficient/wavelength/cfg.sublimation_temperature)
    # In the MATLAB data wavelength decreases: positive contribution uses lambda[i]-lambda[i+1].
    return scales.emission_coefficient * np.sum((wavelength[:-1]-wavelength[1:])*(integrand[:-1]+integrand[1:])/2)
