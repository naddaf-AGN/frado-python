"""Fast numerical consistency checks using supplied MATLAB trajectory headers."""
from pathlib import Path
import numpy as np
import pytest
from frado import Config, Simulation

REF=Path(__file__).resolve().parents[1]/'validation'/'reference'

@pytest.fixture(scope='module')
def sim():
    return Simulation(Config(bh_mass_log=7.0, mdot_log=-0.4,
                             metallicity=5.0, patch_size=3.0))

@pytest.mark.parametrize('radius',(8370,22990,30860))
def test_matlab_launch_conditions(sim,radius):
    m=np.loadtxt(REF/f'T1500_A3_Z5_M7.0_m-0.4_R{radius}.dat',max_rows=1)
    h=float(sim.disk.at(radius))
    vy=sim.scales.kepler_velocity / np.sqrt(np.hypot(radius,h))*sim.scales.velocity_kms_factor
    assert abs(h-m[3]) < 0.00006  # MATLAB position field printed to 4dp
    assert abs(vy-m[5]) < 0.0051   # MATLAB velocity field printed to 2dp


def test_matlab_event_table_structure():
    a=np.loadtxt(REF/'event_T1500_A3_Z5_M7.0_m-0.4.dat')
    assert a.shape==(2803,16)
    assert a[0,0]==7780 and a[-1,0]==35800
    assert np.allclose(np.diff(a[:,0]),10)
    assert set(np.unique(a[:,7]))=={1.0,3.0,5.0}
    assert set(np.unique(a[:,11]))=={2.0,3.0}
