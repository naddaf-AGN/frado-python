from dataclasses import dataclass
import numpy as np
from .config import Config

G = 6.674e-11
MP = 1.6726219e-27
C = 2.99792458e8
SIGMA_T = 6.652e-29
MSUN = 1.988e30
SIGMA_SB = 5.6703e-8
KB = 1.380649e-23
H = 6.62607015e-34
YEAR = 31536000.0

# The two sets of seven sinusoidal opacity fits are transcribed verbatim.
TOTAL_FIT = np.array([[20.21,.3984,.1393],[8.95,.572,2.871],[.2394,2.6,-3.005],
 [.08779,4.203,-1.033],[.04031,5.845,.5002],[.02038,8.644,-1.891],[.01048,9.489,-2.852]])
ABS_FIT = np.array([[17.39,.4623,-.08806],[6.475,.7709,2.154],[.1488,3.057,-4.933],
 [.04243,3.999,.2016],[.04352,5.519,1.389],[.1273,8.843,-3.535],[.1274,8.9,-.6908]])

@dataclass(frozen=True)
class Scales:
    gravitational_radius_m: float
    velocity_kms_factor: float
    kepler_velocity: float
    escape_velocity: float
    orbital_period_coefficient: float
    temperature_coefficient: float
    opacity_argument_offset: float
    gravitational_coefficient: float
    radiative_coefficient: float
    absorption_coefficient: float
    planck_coefficient: float
    emission_coefficient: float
    eddington_accretion_rate: float


def scales(cfg: Config) -> Scales:
    mass = cfg.mass
    rg = G * mass * MSUN / C**2
    ct = 3 * MP * C**5 / (2 * SIGMA_T * SIGMA_SB * G * MSUN * cfg.efficiency)
    psi = .005 * cfg.metallicity
    return Scales(
        gravitational_radius_m=rg,
        velocity_kms_factor=rg * (C / 1e3) / C,
        kepler_velocity=C**3 / (MSUN * mass * G),
        escape_velocity=2**.5 * C**3 / (MSUN * mass * G),
        orbital_period_coefficient=2*np.pi*G*mass*MSUN/C**3,
        temperature_coefficient=ct,
        opacity_argument_offset=.25*(np.log10(ct)+cfg.mdot_log-cfg.bh_mass_log),
        gravitational_coefficient=C**6/(G**2*MSUN**2*mass**2),
        radiative_coefficient=4*H*C*ct*cfg.mdot*(psi/(1+psi))/mass/rg,
        absorption_coefficient=4*H*C**2*ct*cfg.mdot/mass,
        planck_coefficient=H*C/KB,
        emission_coefficient=8*np.pi*H*C**2,
        eddington_accretion_rate=4*np.pi*G*mass*MSUN*MP/(cfg.efficiency*C*SIGMA_T),
    )
