"""Verify numerical agreement between SciPy and compiled radiation backends."""
from dataclasses import replace
import numpy as np
import pytest
from frado import Config, Simulation
from frado.radiation import Radiation
from frado.constants import scales
from frado.shielding import PatchShielding


def _configured(backend):
    c=Config(bh_mass_log=7, mdot_log=-0.4, metallicity=5.,
             initial_radius=8370, final_radius=8370,
             radiation_backend=backend)
    return Simulation(c)


def test_radiation_backend_invalid():
    with pytest.raises(ValueError, match='radiation_backend'):
        Config(radiation_backend='unknown')


def test_numba_matches_reference_at_different_heights_and_shield_branches():
    pytest.importorskip('numba')
    ref=_configured('scipy')
    fast=_configured('numba')
    for rho,z in [(5000,ref.disk.at(5000)),
                  (8370,ref.disk.at(8370)),
                  (22990,ref.disk.at(22990)),
                  (8370,100), (22990,500),
                  (200,100)]:
        state=np.array([rho,0,0,0,z,0],dtype=float)
        ref_force, ref_heat = ref.radiation.evaluate(state)
        fast_force, fast_heat = fast.radiation.evaluate(state)
        np.testing.assert_allclose(fast_force,ref_force,rtol=1e-8,atol=1e-15)
        np.testing.assert_allclose(fast_heat,ref_heat,rtol=1e-8,atol=1e-15)


def test_numba_respects_shielding_empty_bounds():
    pytest.importorskip('numba')
    fast=_configured('numba')
    # At z=0 the patch spans no radial/azimuthal extent.
    f,heat=fast.radiation.evaluate(np.array([8370.,0,0,0,0,0]))
    np.testing.assert_array_equal(f,np.zeros(3))
    assert heat==0


def test_numba_uses_same_launch_screening():
    pytest.importorskip('numba')
    ref=_configured('scipy')
    fast=_configured('numba')
    for r in (5000, 8370, 22990):
        a=ref.screen_launch(r)
        b=fast.screen_launch(r)
        assert a.status==b.status
        np.testing.assert_allclose(a.heating_ratio,b.heating_ratio,rtol=1e-8)


def test_numba_heating_only_agrees_with_full_radiation():
    pytest.importorskip('numba')
    sim=_configured('numba')
    for r in (5000,8370,22990):
        for z in (sim.disk.at(r), 100.):
            s=np.array([r,0,0,0,z,0],dtype=float)
            full=sim.radiation.evaluate(s)[1]
            heating=sim.radiation.absorbed_power(s)
            np.testing.assert_allclose(heating,full,rtol=1e-12,atol=0)
