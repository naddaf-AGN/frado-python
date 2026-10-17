"""Public, immutable FRADO configuration; no changes to physics are hidden here."""
from dataclasses import dataclass
from typing import Literal
import numpy as np

TerminationMode = Literal['disk_impact', 'first_excursion']
NumericalCheckMode = Literal['off', 'screen', 'recheck']


@dataclass(frozen=True)
class Config:
    mdot_log: float = -1.7
    bh_mass_log: float = 8.5
    metallicity: float = 5.0
    initial_radius: float = 2340.0
    final_radius: float = 4000.0
    radius_step: float = 10.0
    sublimation_temperature: float = 1500.0
    efficiency: float = 0.1
    patch_size: float = 3.0
    maximum_orbits: float = 5.0
    output_samples: int = 500_000
    relative_tolerance: float = 1e-7
    absolute_tolerance: float = 1e-6
    inner_radiating_radius: float = 6.0
    outer_blr_radius: float = 1e5
    quad_epsabs: float = 1e-8
    quad_epsrel: float = 1e-12
    # QUADPACK in both modes. 'numba' compiles the scalar integrand only.
    # SciPy is the default radiation backend.
    radiation_backend: Literal['scipy', 'numba'] = 'scipy' 
    # Cap the integrator's internal step to resolve short disk crossings.
    max_step_orbit_fraction: float = 0.02

    # Physical stopping policy.
    termination_mode: TerminationMode = 'disk_impact'
    # A close approach occurs when 0 < z-H(rho) <= near_miss_fraction * H(rho)
    # at a local *minimum* of the vertical clearance following an excursion.
    near_miss_fraction: float = 0.05

    # Optional conservative screening of grazing trajectories and a targeted
    # second integration using finer ODE tolerances / smaller maximum step.
    numerical_check: NumericalCheckMode = 'screen'
    grazing_screen_fraction: float = 0.001
    grazing_speed_threshold: float = 0.3
    recheck_tolerance_factor: float = 0.1
    recheck_step_factor: float = 0.5
    recheck_time_tolerance_orbits: float = 0.005

    def __post_init__(self):
        if self.radiation_backend not in ('scipy', 'numba'):
            raise ValueError('radiation_backend must be scipy or numba')
        if self.numerical_check not in ('off', 'screen', 'recheck'):
            raise ValueError('numerical_check must be off, screen, or recheck')
        if self.termination_mode not in ('disk_impact', 'first_excursion'):
            raise ValueError('termination_mode must be disk_impact or first_excursion')
        if not (0.0 < self.max_step_orbit_fraction <= 1):
            raise ValueError('max_step_orbit_fraction must lie in (0,1]')
        if not (0.0 < self.near_miss_fraction < 1):
            raise ValueError('near_miss_fraction must lie strictly between 0 and 1')
        if not (0 < self.grazing_screen_fraction < 1):
            raise ValueError('grazing_screen_fraction must lie in (0,1)')
        if not (0 < self.recheck_tolerance_factor < 1 and 0 < self.recheck_step_factor < 1):
            raise ValueError('Recheck refinement factors must lie in (0,1)')
        positive = ('grazing_speed_threshold', 'recheck_time_tolerance_orbits', 'metallicity', 'sublimation_temperature', 'efficiency', 'patch_size',
                    'maximum_orbits', 'radius_step', 'relative_tolerance', 'absolute_tolerance',
                    'quad_epsabs', 'quad_epsrel')
        for name in positive:
            val = getattr(self, name)
            if not np.isfinite(val) or val <= 0:
                raise ValueError(f'{name} must be positive and finite')
        if not (np.isfinite(self.initial_radius) and np.isfinite(self.final_radius)):
            raise ValueError('Launch radii must be finite')
        if self.initial_radius < 100 or self.final_radius > 100_000 or self.final_radius < self.initial_radius:
            raise ValueError('Launch radii must satisfy 100 <= initial <= final <= 100000 Rg')
        if self.inner_radiating_radius <= 0 or self.outer_blr_radius <= 0:
            raise ValueError('Radiating inner radius and BLR outer radius must be positive')
        if not isinstance(self.output_samples, (int, np.integer)) or self.output_samples < 2:
            raise ValueError('output_samples must be an integer >= 2')
        if not (6 <= self.bh_mass_log <= 10 and -2 <= self.mdot_log <= 0):
            raise ValueError('Disk fits require 6<=bh_mass_log<=10 and -2<=mdot_log<=0')

    @property
    def mdot(self):
        return 10.0 ** self.mdot_log

    @property
    def mass(self):
        return 10.0 ** self.bh_mass_log

    @property
    def launch_radii(self):
        steps = (self.final_radius - self.initial_radius) / self.radius_step
        if steps < 0 or not np.isclose(steps, round(steps)):
            raise ValueError('Radius interval must contain an integer number of positive steps')
        return np.linspace(self.initial_radius, self.final_radius, int(round(steps)) + 1)
