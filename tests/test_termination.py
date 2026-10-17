"""Analytic, inexpensive orbital-clearance tests independent of disk radiation."""
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from frado import Config
from frado.integrator import integrate_trajectory
from frado.termination import clearance,clearance_rate


class ConstantDisk:
    def __init__(self):
        self.radius=np.array([100.,100000.])
        self.height=np.array([1.,1.])
        self.interpolator=PchipInterpolator(self.radius,self.height)
    def at(self,r): return self.interpolator(r)


class FakeRadiation:
    emission_power=100.
    def evaluate(self,s): return np.zeros(3),0.0


class DrivenOscillation:
    """With z(0)=1, z'(0)=0: z-1=A(1-cos(t)) + b*t**2."""
    def __init__(self,b): self.b=b
    def derivative(self,t,s,dusty):
        return np.array([s[1],0.,s[3],0.,s[5],np.cos(t)+2*self.b])


def analytic_trajectory(mode,b,threshold=.1):
    r=1000.
    cfg=Config(initial_radius=r,final_radius=r, maximum_orbits=3,
               output_samples=2500, termination_mode=mode,
               near_miss_fraction=threshold,relative_tolerance=1e-10,
               absolute_tolerance=1e-11)
    disk=ConstantDisk()
    scales=SimpleNamespace(kepler_velocity=0., escape_velocity=1e10,
                           orbital_period_coefficient=2*np.pi/(np.hypot(r,1)**1.5))
    return integrate_trajectory(cfg,scales,disk,FakeRadiation(),DrivenOscillation(b),r)


def test_surface_gradient_uses_radial_motion():
    d=ConstantDisk()
    s=np.array([1000.,1.,0.,0.,1.1,-2.])
    assert np.isclose(clearance(s,d),.1)
    assert np.isclose(clearance_rate(s,d),-2.)


def test_actual_disk_impact_stops_on_downward_crossing():
    traj=analytic_trajectory('disk_impact',b=-.001)
    assert traj.termination_reason=='DISK_IMPACT'
    assert 5.<traj.time_s[-1]<7.
    assert abs(traj.final_state[4]-1.)<1e-5
    assert traj.events[-1][0]==2


def test_first_excursion_stops_on_near_miss():
    traj=analytic_trajectory('first_excursion',b=.001,threshold=.1)
    assert traj.termination_reason=='NEAR_MISS'
    assert 5.5<traj.time_s[-1]<7.
    assert 0.<traj.final_state[4]-1.<.1
    assert any(e.name=='near_miss' for e in traj.event_records)


def test_disk_impact_policy_allows_another_orbit_if_missed():
    traj=analytic_trajectory('disk_impact',b=.001,threshold=.1)
    assert traj.termination_reason=='MAX_TIME_REACHED'
    assert traj.time_s[-1]>15.
    assert not any(e.name=='near_miss' for e in traj.event_records)


def test_first_excursion_does_not_stop_for_far_approach():
    traj=analytic_trajectory('first_excursion',b=.001,threshold=.01)
    assert traj.termination_reason=='MAX_TIME_REACHED'
    assert traj.time_s[-1]>15.


def test_tangential_contact_is_detected_as_impact():
    traj=analytic_trajectory('disk_impact',b=0.)
    assert traj.termination_reason=='DISK_IMPACT'
    assert abs(traj.time_s[-1]-2*np.pi)<.005


def test_invalid_termination_modes_and_thresholds():
    with pytest.raises(ValueError): Config(termination_mode='matlab_legacy')
    with pytest.raises(ValueError): Config(near_miss_fraction=1.5)
