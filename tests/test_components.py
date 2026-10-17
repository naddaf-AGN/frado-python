import numpy as np
from frado import Config,Simulation
from frado.constants import scales,YEAR
from frado.disk_geometry import build_disk
from frado.dust import emission_power
from frado.shielding import PatchShielding


def test_defaults_and_units():
    c=Config()
    s=scales(c)
    assert len(c.launch_radii)==167
    assert np.isclose(s.velocity_kms_factor,s.gravitational_radius_m/1000)
    assert np.isclose(s.orbital_period_coefficient,2*np.pi*s.gravitational_radius_m/299792458)
    assert YEAR==31536000


def test_disk_interpolation_and_geometry():
    d=build_disk(Config())
    assert d.radius.size==20000
    assert 0<float(d.at(2340))<2340
    assert np.isfinite(d.at(4000))


def test_dust_positive():
    c=Config(); p=emission_power(c,scales(c))
    assert np.isfinite(p) and p>0


def test_shielding_branch():
    c=Config()
    lo,hi,a,b=PatchShielding(c).bounds(2000,100)
    assert np.isclose(lo,1700) and np.isclose(hi,2300)
    assert 0<b<=np.pi


def test_initial_acceleration_dustless():
    sim=Simulation(Config())
    z=float(sim.disk.at(2340))
    state=np.array([2340.,0.,0.,sim.scales.kepler_velocity/np.sqrt(np.hypot(2340,z)),z,0.])
    val=sim.dynamics.derivative(0,state,False)
    assert np.all(np.isfinite(val)) and val[1]<0 and val[5]<0
