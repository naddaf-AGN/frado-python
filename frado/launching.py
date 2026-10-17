"""Physical dust-survival check before cloud launching.

A point at the disk surface is *eligible for dusty launching* only when its
absorbed power is strictly less than the dust thermal-emission threshold.
This is a necessary condition for dusty launching, not a guarantee that the
cloud acquires a positive vertical acceleration or escapes the disk.
"""
from dataclasses import dataclass
import numpy as np


class NoDustAtLaunchError(ValueError):
    """Dust is already at/above the sublimation threshold at the launch point."""

    def __init__(self, radius, absorbed_power, emission_power):
        self.radius = float(radius)
        self.absorbed_power = float(absorbed_power)
        self.emission_power = float(emission_power)
        super().__init__(
            f'No dusty launch at R={radius:g} Rg: absorbed power '
            f'{absorbed_power:.8g} >= thermal emission {emission_power:.8g}'
        )


@dataclass(frozen=True)
class LaunchAssessment:
    radius_rg: float
    height_rg: float
    absorbed_power: float
    emission_power: float
    heating_ratio: float
    can_launch: bool

    @property
    def status(self):
        return 'LAUNCH_ALLOWED' if self.can_launch else 'NO_DUST_AT_LAUNCH'


def assess_launch(radius_rg, geometry, radiation):
    """Evaluate the *same* dust heating law used by the trajectory integrator.

    Uses the disk-surface initial position x=rho, y=0, z=H(rho).
    The threshold is inclusive: absorbed >= emitted disallows a dusty launch.
    Each requested radius is assessed individually: no assumption of a
    contiguous launching region or monotonic heating-versus-radius is made.
    """
    rho = float(radius_rg)
    height = float(geometry.at(rho))
    if not np.isfinite(height):
        raise ValueError(f'Nonfinite disk height at launch radius {rho:g} Rg')
    state = np.array([rho, 0., 0., 0., height, 0.], dtype=float)
    # A backend may implement efficient heating-only quadrature. Retain
    # compatibility with independent radiation prescriptions using evaluate().
    if callable(getattr(radiation, 'absorbed_power', None)):
        absorbed = radiation.absorbed_power(state)
    else:
        _, absorbed = radiation.evaluate(state)
    absorbed = float(absorbed)
    emitted = float(radiation.emission_power)
    if not np.isfinite(absorbed) or absorbed < 0:
        raise ValueError(f'Invalid dust absorption power at {rho:g} Rg: {absorbed}')
    if not np.isfinite(emitted) or emitted <= 0:
        raise ValueError(f'Invalid dust emission power: {emitted}')
    return LaunchAssessment(rho, height, absorbed, emitted, absorbed / emitted,
                            absorbed < emitted)
