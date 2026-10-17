"""Event-driven FRADO integration with explicit sublimation and stop policies.

Internal state order: [x, vx, y, vy, z, vz]; positions in Rg,
velocities in Rg/s, times in seconds. This keeps the validated force law.
"""
from dataclasses import dataclass
from typing import Optional
import numpy as np
from scipy.integrate import solve_ivp
from .termination import TerminationEvents, clearance
from .launching import NoDustAtLaunchError

# Event numbers 1, 3, 5 retain their MATLAB diagnostic meanings;
# 2 denotes disk impact and 4 denotes a near miss.
EVENT_IDS = {
    'vertical_turning': 1, 'disk_impact': 2, 'blr_boundary': 3,
    'near_miss': 4, 'escape_speed': 5, 'sublimation': 6,
    'clearance_maximum': 7, 'clearance_minimum': 8,
    'surface_departure': 9,
}


@dataclass(frozen=True)
class EventRecord:
    name: str
    time_s: float
    state: np.ndarray

    @property
    def id(self):
        return EVENT_IDS[self.name]


@dataclass
class Trajectory:
    time_s: np.ndarray
    states: np.ndarray
    dusty: np.ndarray
    # Compatibility representation: (numeric event ID, time in seconds, 6-D state).
    events: list
    sublimation_s: float
    sublimation_state: Optional[np.ndarray]
    orbital_period_s: float
    termination_reason: str
    event_records: list
    initial_radius: float

    @property
    def final_state(self):
        return self.states[-1]


def _event(fn, *, direction=0, terminal=False):
    fn.direction = direction
    fn.terminal = terminal
    return fn


def integrate_trajectory(cfg, scales, geometry, radiation, dynamics, initial_rho,
                         launch_absorbed_power=None):
    z0 = float(geometry.at(initial_rho))
    if not np.isfinite(z0):
        raise ValueError('Launch radius lies outside the disk geometry domain')
    initial_r = np.hypot(initial_rho, z0)
    state0 = np.array([initial_rho, 0., 0., scales.kepler_velocity/np.sqrt(initial_r), z0, 0.])
    orbit = scales.orbital_period_coefficient * initial_r**1.5
    t_end = float(cfg.maximum_orbits * orbit)

    policy = TerminationEvents(cfg, geometry)
    # The segmented state machine only enables a minimum *after* the preceding
    # maximum, so t=0 (cloud launched on the disk) is not a spurious near miss.
    # Track vertical-clearance extrema in either mode; even a tangential
    # disk contact (no sign crossing) should count as an impact.
    phase = 'max'
    begin, state, dusty = 0., state0, True
    departed = False
    # Numerical arming boundary, not a physical alteration to H(rho).
    departure_epsilon = 1e-6 * max(1., z0)
    chunks, records = [], []
    sublimation_time = np.inf
    sublimation_state = None
    reason = 'MAX_TIME_REACHED'

    # No cloud exists to integrate if dust is already sublimated at launch.
    # Sweep runs can pass the pre-screened value and avoid duplicate quadrature.
    absorbed_at_launch = (float(launch_absorbed_power)
                          if launch_absorbed_power is not None
                          else float(radiation.absorbed_power(state0) if callable(getattr(radiation,'absorbed_power',None))
                                     else radiation.evaluate(state0)[1]))
    if not np.isfinite(absorbed_at_launch):
        raise ValueError("Nonfinite dust absorption at launch")
    if absorbed_at_launch >= radiation.emission_power:
        raise NoDustAtLaunchError(initial_rho, absorbed_at_launch,
                                  radiation.emission_power)

    for iteration in range(20_000):  # Guard against nonprogressing event restarts.
        if begin >= t_end:
            break

        def turning(t, s): return s[5]
        def impact(t, s): return policy.impact(t, s)
        def departure(t, s): return policy.impact(t, s) - departure_epsilon
        def downward_from_launch(t, s): return policy.impact(t, s) + departure_epsilon
        def outer(t, s): return np.linalg.norm(s[[0, 2, 4]]) - cfg.outer_blr_radius
        def escape(t, s):
            r = np.linalg.norm(s[[0, 2, 4]])
            return np.linalg.norm(s[[1, 3, 5]]) - scales.escape_velocity/np.sqrt(r)

        names = ['vertical_turning', 'blr_boundary', 'escape_speed']
        events = [_event(turning, direction=-1),
                  _event(outer, direction=1, terminal=True),
                  _event(escape, direction=1)]
        if not departed:
            # At t=0, g=z-H=0 by construction. A disk-impact root at
            # t=0 would be a false *return*. Arm impact detection only once
            # the cloud has demonstrably departed from the surface.
            names.extend(['surface_departure', 'launch_downward'])
            events.extend([_event(departure,direction=1,terminal=True),
                           _event(downward_from_launch,direction=-1,terminal=True)])
        else:
            names.append('disk_impact')
            events.append(_event(impact,direction=-1,terminal=True))
            if phase is not None:
                if phase == 'max':
                    def gap_extremum(t, s): return policy.max_clearance(t, s)
                    names.append('clearance_maximum')
                else:
                    def gap_extremum(t, s): return policy.min_clearance(t, s)
                    names.append('clearance_minimum')
                events.append(_event(gap_extremum, direction=-1 if phase=='max' else 1, terminal=True))

        if dusty:
            def sublimation(t, s):
                return (radiation.absorbed_power(s) if callable(getattr(radiation,'absorbed_power',None))
                        else radiation.evaluate(s)[1]) - radiation.emission_power
            names.append('sublimation')
            events.append(_event(sublimation, direction=1, terminal=True))

        sol = solve_ivp(lambda t,s: dynamics.derivative(t,s,dusty), (begin,t_end), state,
                        method='BDF', atol=cfg.absolute_tolerance,
                        rtol=cfg.relative_tolerance, events=events, dense_output=True,
                        max_step=cfg.max_step_orbit_fraction * orbit)
        if not sol.success:
            raise RuntimeError(f'ODE solver failed at launch R={initial_rho}: {sol.message}')
        if sol.t[-1] <= begin and begin < t_end:
            raise RuntimeError(f'Nonadvancing ODE event at t={begin}, R={initial_rho}')

        chunks.append((begin, float(sol.t[-1]), sol.sol, dusty))
        occurred = {}
        for name, event_times, event_states in zip(names, sol.t_events, sol.y_events):
            for t, s in zip(event_times, event_states):
                t = float(t)
                occurred[name] = (t, s.copy())
                # Internal launch-downward is mapped to the disk-impact code.
                records.append(EventRecord('disk_impact' if name=='launch_downward' else name, t, s.copy()))

        # A physical termination takes precedence over bookkeeping events.
        if 'disk_impact' in occurred or 'launch_downward' in occurred:
            reason = 'DISK_IMPACT'
            break
        if 'blr_boundary' in occurred:
            reason = 'BLR_BOUNDARY'
            break

        # Sublimation is a solver restart, never a terminal event for motion.
        if 'surface_departure' in occurred:
            begin, state = occurred['surface_departure']
            departed = True
            continue

        if 'sublimation' in occurred:
            sublimation_time, sublimation_state = occurred['sublimation']
            begin, state, dusty = sublimation_time, sublimation_state.copy(), False
            continue

        if 'clearance_maximum' in occurred:
            begin, state = occurred['clearance_maximum']
            phase = 'min'
            continue
        if 'clearance_minimum' in occurred:
            begin, state = occurred['clearance_minimum']
            rho_at_min = np.hypot(state[0], state[2])
            h_at_min = float(geometry.at(rho_at_min))
            gap_at_min = float(clearance(state, geometry))
            # Near-tangential physical contact can have no sign-changing root.
            # This tiny numerical tolerance is distinct from the user's
            # (much larger) *physical* near-miss threshold.
            if gap_at_min <= 1e-8 * h_at_min:
                reason = 'DISK_IMPACT'
                records.append(EventRecord('disk_impact', begin, state.copy()))
                break
            if cfg.termination_mode == 'first_excursion' and policy.is_near_miss(state):
                reason = 'NEAR_MISS'
                records.append(EventRecord('near_miss', begin, state.copy()))
                break
            phase = 'max'  # Far from the disk: let it make another excursion.
            continue
        # Normal completion at maximum simulation time.
        break
    else:
        raise RuntimeError('Exceeded maximum event-restart iterations')

    final = chunks[-1][1] if chunks else 0.0
    dt = t_end / (int(cfg.output_samples)-1)
    count = min(int(cfg.output_samples), int(np.floor(final/dt)) + 1)
    times = dt * np.arange(max(count, 1), dtype=float)
    if times[-1] < final - max(1e-9, abs(final)*1e-13):
        times = np.append(times, final)
    else:
        times[-1] = final

    states = np.empty((len(times), 6), dtype=float)
    if chunks:
        for i,(t0,t1,interp,is_dusty) in enumerate(chunks):
            mask = ((times >= t0) if i == 0 else (times > t0)) & (times <= t1)
            if np.any(mask):
                states[mask] = interp(times[mask]).T
    else:
        states[0] = state0

    # Correct phase assignment at the exact sublimation instant.
    dust = times < sublimation_time
    records.sort(key=lambda e: (e.time_s, e.id))
    compat = [(e.id, e.time_s, e.state.copy()) for e in records]
    return Trajectory(times, states, dust, compat, sublimation_time,
                      sublimation_state, orbit, reason, records, float(initial_rho))
