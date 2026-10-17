"""Conservative numerical screen and optional targeted ODE-convergence recheck.

A screen is NOT proof of convergence. A recheck tests sensitivity to one
specified refinement and is NOT an independent physical validation.
"""
from dataclasses import dataclass, replace
import numpy as np
from .termination import clearance, clearance_rate


@dataclass(frozen=True)
class ReliabilityAssessment:
    status: str
    flags: tuple[str, ...] = ()
    minimum_clearance_fraction: float | None = None
    normalized_impact_speed: float | None = None
    recheck_reason: str | None = None
    recheck_time_difference_fraction: float | None = None
    note: str = ''


def screen_trajectory(result, geometry, cfg):
    """Flag fragile geometric contacts; do not relabel physical event outcomes.

    Two distinct numerical warnings:
      - a minimum gap near zero or a policy threshold;
      - an impact with unusually small surface-crossing rate.

    Heights are local H(rho), not H at the launch location.
    """
    min_fraction = None
    impact_rate = None
    flags = []
    for event in result.event_records:
        if event.name != 'clearance_minimum':
            continue
        rho = np.hypot(event.state[0], event.state[2])
        h = float(geometry.at(rho))
        if not np.isfinite(h) or h <= 0:
            continue
        fraction = float(clearance(event.state, geometry) / h)
        if min_fraction is None or fraction < min_fraction:
            min_fraction = fraction
        if abs(fraction) <= cfg.grazing_screen_fraction:
            flags.append('MINIMUM_GAP_NEAR_SURFACE')
        if cfg.termination_mode == 'first_excursion' and abs(fraction-cfg.near_miss_fraction) <= cfg.grazing_screen_fraction:
            flags.append('MINIMUM_GAP_NEAR_NEAR_MISS_THRESHOLD')
    if result.termination_reason == 'DISK_IMPACT':
        s = result.final_state
        rho = np.hypot(s[0], s[2])
        h = float(geometry.at(rho))
        if np.isfinite(h) and h > 0:
            impact_rate = float(abs(clearance_rate(s, geometry)) * result.orbital_period_s / h)
            if impact_rate <= cfg.grazing_speed_threshold:
                flags.append('SLOW_SURFACE_CROSSING')
    flags = tuple(dict.fromkeys(flags))
    return ReliabilityAssessment(
        status='SCREEN_FLAGGED' if flags else 'SCREEN_CLEAR',
        flags=flags,
        minimum_clearance_fraction=min_fraction,
        normalized_impact_speed=impact_rate,
        note='Heuristic screening only; SCREEN_CLEAR does not establish numerical convergence.'
    )


def assess_trajectory(result, geometry, cfg, refined_trajectory=None, recheck_error=None):
    """Screen, optionally compare to a second integration with refined ODE controls."""
    if cfg.numerical_check == 'off':
        return ReliabilityAssessment(status='NOT_CHECKED', note='Numerical screening disabled.')
    assessment = screen_trajectory(result, geometry, cfg)
    if cfg.numerical_check != 'recheck' or assessment.status != 'SCREEN_FLAGGED':
        return assessment
    if recheck_error is not None:
        return replace(assessment, status='RECHECK_FAILED', note=str(recheck_error))
    if refined_trajectory is None:
        return assessment
    dt_frac = abs(float(result.time_s[-1]) - float(refined_trajectory.time_s[-1])) / max(
        float(result.orbital_period_s), np.finfo(float).tiny)
    same_reason = refined_trajectory.termination_reason == result.termination_reason
    agree = same_reason and dt_frac <= cfg.recheck_time_tolerance_orbits
    return replace(assessment,
                   status='RECHECK_AGREES' if agree else 'RECHECK_DISAGREES',
                   recheck_reason=refined_trajectory.termination_reason,
                   recheck_time_difference_fraction=dt_frac,
                   note='One tighter ODE integration; other error sources may remain.')
