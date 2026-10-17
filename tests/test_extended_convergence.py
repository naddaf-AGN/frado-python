"""Lightweight regression guards for event handling near disk grazing/sublimation.

These independently specified analytic trajectories do not validate the FRADO
radiative force, which is assessed in validation/extended_convergence.py.
"""
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from frado import Config
from frado.integrator import integrate_trajectory
from frado.termination import clearance_rate


class FlatDisk:
    def __init__(self):
        self.interpolator=PchipInterpolator([100.,100_000.],[1.,1.])
    def at(self, r): return self.interpolator(r)


class ToyRadiation:
    def __init__(self, trigger_sublimation=False):
        self.emission_power = 0.5 if trigger_sublimation else 1e5
        self.trigger = trigger_sublimation
    def evaluate(self,state):
        # In test with trigger, sublimation exactly at z=1.5.
        return np.zeros(3), (state[4]-1. if self.trigger else 0.)


class ToyDynamics:
    def __init__(self,b=0.): self.b=b
    def derivative(self,t,s,dusty):
        # z(t) = 2-cos(t)+b*t*t; 100% analytic.
        return np.array([s[1],0.,s[3],0.,s[5],np.cos(t)+2*self.b])


def toy(b=0.,mode='disk_impact',trigger=False,rtol=1e-9,atol=1e-10):
    R=1000.
    config=Config(initial_radius=R,final_radius=R,termination_mode=mode,
                  maximum_orbits=2.,output_samples=1001,
                  relative_tolerance=rtol,absolute_tolerance=atol)
    scale=SimpleNamespace(kepler_velocity=0.,escape_velocity=1e9,
                          orbital_period_coefficient=2*np.pi/(np.hypot(R,1.)**1.5))
    return integrate_trajectory(config,scale,FlatDisk(),ToyRadiation(trigger),ToyDynamics(b),R)


@pytest.mark.parametrize('b,mode,expected',[
    (-1e-4,'disk_impact','DISK_IMPACT'),
    (0.,'disk_impact','DISK_IMPACT'),
    (1e-4,'disk_impact','MAX_TIME_REACHED'),
    (1e-4,'first_excursion','NEAR_MISS'),
])
def test_grazing_classification_is_stable_under_tighter_solver_settings(b,mode,expected):
    for rtol,atol in ((1e-7,1e-8),(1e-10,1e-11)):
        result=toy(b,mode,rtol=rtol,atol=atol)
        assert result.termination_reason==expected
        assert result.time_s[-1]>0.
        if expected=='DISK_IMPACT':
            # Near a double root, finite solver error can leave a tiny negative
            # gap at the g'(t)=0 extremum; tighten tolerance before interpreting it.
            permitted = 1e-5 if rtol >= 1e-7 else 1e-7
            assert abs(result.final_state[4]-1.)<permitted
        if expected=='NEAR_MISS':
            assert 0.<result.final_state[4]-1.<.05
        if expected=='MAX_TIME_REACHED':
            assert np.isclose(result.time_s[-1],4*np.pi,atol=1e-6)


def test_explicit_sublimation_root_and_dust_switch():
    result=toy(b=0.,trigger=True,rtol=1e-10,atol=1e-11)
    assert abs(result.sublimation_s-np.pi/3)<1e-7
    assert result.sublimation_state is not None
    assert abs(result.sublimation_state[4]-1.5)<1e-8
    assert sum(e.name=='sublimation' for e in result.event_records)==1
    assert np.any(result.dusty) and np.any(~result.dusty)
    assert np.all(result.dusty == (result.time_s < result.sublimation_s))
    assert result.termination_reason=='DISK_IMPACT'
    assert abs(result.time_s[-1]-2*np.pi)<1e-6


def test_clearance_rate_includes_radial_disk_slope():
    class SlopingDisk:
        interpolator=PchipInterpolator([100.,1000.,5000.], [2.,11.,51.])
    s=np.array([1000.,3.,0.,0.,15.,-2.])
    # H'(rho)=0.01; drho/dt=3 -> g'(t) = vz - 0.01*3
    assert np.isclose(clearance_rate(s,SlopingDisk()),-2.03,atol=1e-12)
