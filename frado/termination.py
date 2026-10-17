"""Configurable termination policies, independent of disk/radiation dynamics.

The event for actual impact is the *downward* zero crossing of the geometrical
surface separation: g = z-H(rho). It is intentionally not the legacy
1.005*H threshold. A first-excursion near miss is a local minimum of this
same separation, detected via dg/dt changing from negative to positive,
conditional on the cloud having first gone through a separation maximum.
"""
import numpy as np


def clearance(state, geometry):
    rho = np.hypot(state[0], state[2])
    return float(state[4] - geometry.at(rho))


def clearance_rate(state, geometry):
    x, vx, y, vy, _, vz = state
    rho = np.hypot(x, y)
    drho_dt = (x * vx + y * vy) / rho if rho > 0 else 0.0
    return float(vz - geometry.interpolator.derivative()(rho) * drho_dt)


class TerminationEvents:
    def __init__(self, cfg, geometry):
        self.cfg, self.geometry = cfg, geometry

    def impact(self, t, state):
        return clearance(state, self.geometry)

    def max_clearance(self, t, state):
        return clearance_rate(state, self.geometry)

    def min_clearance(self, t, state):
        return clearance_rate(state, self.geometry)

    def is_near_miss(self, state):
        """Close but nonimpacting; threshold is relative to *local* H(rho)."""
        rho = np.hypot(state[0], state[2])
        h = float(self.geometry.at(rho))
        g = float(state[4] - h)
        return h > 0 and g > 0 and g <= self.cfg.near_miss_fraction * h
