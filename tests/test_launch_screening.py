"""Dusty-launch eligibility and per-radius sweep regression tests."""
import csv
import json
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from frado import Config, Simulation, NoDustAtLaunchError
from frado.launching import assess_launch
from frado.integrator import Trajectory, integrate_trajectory


class FlatDisk:
    def __init__(self):
        self.interpolator = PchipInterpolator([100.,100_000.], [1.,1.])

    def at(self, rho):
        return self.interpolator(rho)


class DeterministicRadiation:
    emission_power = 1.0

    def __init__(self, heating):
        self.heating = heating
        self.calls = []

    def evaluate(self, state):
        self.calls.append(np.array(state).copy())
        return np.zeros(3), self.heating(float(state[0]))


def fake_trajectory(rho):
    t = np.array([0., 100.])
    s = np.array([[rho,0.,0.,1e-4,1.,0.], [rho,0.,0.,1e-4,1.,-1e-7]])
    return Trajectory(t, s, np.array([True,True]), [], np.inf, None,
                      1000., 'DISK_IMPACT', [], rho)


def test_thermal_threshold_inclusive_and_correct_surface_position():
    disk = FlatDisk()
    for power, allowed in ((.999999,True), (1.,False), (1.000001,False)):
        rad = DeterministicRadiation(lambda r: power)
        status = assess_launch(500., disk, rad)
        assert status.can_launch == allowed
        assert status.status == ('LAUNCH_ALLOWED' if allowed else 'NO_DUST_AT_LAUNCH')
        np.testing.assert_equal(rad.calls[0], [500.,0.,0.,0.,1.,0.])
        assert status.heating_ratio == power


def test_invalid_heating_not_misclassified_as_physical_no_launch():
    for value in (float('nan'), float('inf'), -1.):
        with pytest.raises(ValueError, match='Invalid dust absorption'):
            assess_launch(500., FlatDisk(), DeterministicRadiation(lambda r: value))


def test_sweep_skips_only_unsurvivable_radii_and_records_each(tmp_path, monkeypatch):
    cfg = Config(initial_radius=500., final_radius=520.,radius_step=5.,
                 numerical_check='off')
    sim = Simulation.__new__(Simulation)
    sim.cfg, sim.disk = cfg, FlatDisk()
    sim.scales = type('Scales', (), {
        'orbital_period_coefficient':2., 'velocity_kms_factor':100.})()
    # Deliberately non-monotonic eligible radii to verify independent checks.
    sim.radiation = DeterministicRadiation(lambda r:
        {500:2.0,505:.8,510:1.0,515:.2,520:1.7}[int(r)])
    integrated = []
    def fake(self, rho, launch_screen=None):
        assert launch_screen is not None
        assert launch_screen.can_launch
        integrated.append(rho)
        return fake_trajectory(rho)
    monkeypatch.setattr(Simulation,'trajectory',fake)
    results = sim.run(tmp_path, on_error='continue', collect_results=True)
    assert integrated == [505.,515.]
    assert len(results) == 2
    with (tmp_path/'launch_screening.csv').open() as f:
        scr = list(csv.DictReader(f))
    assert len(scr) == 5
    assert [s['status'] for s in scr] == [
        'NO_DUST_AT_LAUNCH', 'LAUNCH_ALLOWED', 'NO_DUST_AT_LAUNCH',
        'LAUNCH_ALLOWED', 'NO_DUST_AT_LAUNCH']
    with (tmp_path/'summary.csv').open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 5
    assert [r['termination_reason'] for r in rows] == [
        'NO_DUST_AT_LAUNCH', 'DISK_IMPACT', 'NO_DUST_AT_LAUNCH',
        'DISK_IMPACT', 'NO_DUST_AT_LAUNCH']
    assert all(r['flight_time_years'] == '' and r['trajectory_file'] == ''
               for r in rows if r['termination_reason']=='NO_DUST_AT_LAUNCH')
    assert all(r['launch_heating_ratio'] for r in rows)
    with (tmp_path/'events.csv').open() as f:
        assert len(list(csv.DictReader(f)))==0
    assert len(list((tmp_path/'Trajectory').rglob('*.dat'))) == 2
    report = json.loads((tmp_path/'run_report.json').read_text())
    assert report['attempted_clouds'] == 5
    assert report['eligible_launches'] == 2
    assert report['skipped_no_dust_at_launch'] == 3
    assert report['eligible_requested_radius_min_rg'] == 505.
    assert report['eligible_requested_radius_max_rg'] == 515.


def test_all_radii_no_dust_still_generates_complete_metadata(tmp_path):
    cfg=Config(initial_radius=500., final_radius=510., radius_step=5.)
    sim=Simulation.__new__(Simulation)
    sim.cfg, sim.disk = cfg, FlatDisk()
    sim.scales = type('Scales', (), {'orbital_period_coefficient':2.})()
    sim.radiation = DeterministicRadiation(lambda r:1e6)
    result=sim.run(tmp_path, collect_results=True)
    assert result==[]
    report=json.loads((tmp_path/'run_report.json').read_text())
    assert report['attempted_clouds']==3
    assert report['eligible_launches']==0
    assert report['skipped_no_dust_at_launch']==3
    assert report['eligible_requested_radius_min_rg'] is None
    with (tmp_path/'summary.csv').open() as f:
        assert all(row['termination_reason']=='NO_DUST_AT_LAUNCH'
                   for row in csv.DictReader(f))


def test_direct_integrator_does_not_launch_dust_free_cloud():
    cfg = Config(initial_radius=500, final_radius=500, maximum_orbits=.1)
    rad = DeterministicRadiation(lambda r:2.)
    scales = type('Scales', (), {
        'kepler_velocity':1., 'orbital_period_coefficient':1.,
        'escape_velocity':1.})()
    class NeverNeededDynamics:
        def derivative(self, t, y, dusty):
            raise AssertionError('ODE must not begin for no-launch radius')
    with pytest.raises(NoDustAtLaunchError):
        integrate_trajectory(cfg, scales, FlatDisk(), rad,
                             NeverNeededDynamics(), 500.)
