"""Numerical warning logic, radius specification, streaming and supplementary output."""
import csv
import json
from dataclasses import replace
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from frado import Config
from frado.integrator import EventRecord,Trajectory
from frado.reliability import screen_trajectory, assess_trajectory
from frado.simulation import Simulation
from run import build_parser, parse_config


class FlatDisk:
    def __init__(self):
        self.interpolator=PchipInterpolator([100.,100000.],[1.,1.])
    def at(self,r):return self.interpolator(r)


def trajectory(radius=500., reason='DISK_IMPACT', event=None, end_s=10., vz=0.):
    s0=np.array([radius,0.,0.,0.0001,1.,0.])
    s1=np.array([radius,0.,0.,0.0001,1.+(0. if event is None else 1e-9),vz])
    evs=[event] if event is not None else []
    return Trajectory(np.array([0.,end_s]),np.vstack([s0,s1]),
        np.array([True,True]),[],np.inf,None,100.,reason,evs,radius)


def test_launch_radius_grid_and_validation():
    c=Config(initial_radius=1000.,final_radius=1020.,radius_step=5.)
    np.testing.assert_array_equal(c.launch_radii,[1000,1005,1010,1015,1020])
    with pytest.raises(ValueError): Config(initial_radius=1000.,final_radius=1021.,radius_step=5.).launch_radii


def test_cli_json_config_override_and_dry_run(tmp_path):
    path=tmp_path/'cfg.json'
    path.write_text(json.dumps({'initial_radius':1000.,'final_radius':1020.,'radius_step':5.,
                                 'numerical_check':'recheck'}))
    p=build_parser()
    c=parse_config(p.parse_args(['--config',str(path),'--radius-start','1100','--radius-end','1110',
                                 '--radius-step','5','--numerical-check','screen']))
    np.testing.assert_array_equal(c.launch_radii,[1100.,1105.,1110.])
    assert c.numerical_check=='screen'
    assert parse_config(p.parse_args(['--radius','8370'])).launch_radii.tolist()==[8370.]
    with pytest.raises(ValueError):parse_config(p.parse_args(['--radius-end','8400']))


def test_slow_grazing_contact_flagged():
    cfg=Config(initial_radius=500,final_radius=500,numerical_check='screen')
    obj=trajectory(vz=0.)
    ass=screen_trajectory(obj,FlatDisk(),cfg)
    assert ass.status=='SCREEN_FLAGGED'
    assert 'SLOW_SURFACE_CROSSING' in ass.flags
    assert ass.normalized_impact_speed==0


def test_non_grazing_contact_screen_clear():
    cfg=Config(initial_radius=500,final_radius=500)
    obj=trajectory(vz=-.1)
    assert screen_trajectory(obj,FlatDisk(),cfg).status=='SCREEN_CLEAR'


def test_minimum_and_near_miss_policy_cutoff_flagged():
    disk=FlatDisk()
    cfg=Config(initial_radius=500,final_radius=500,
               termination_mode='first_excursion',near_miss_fraction=.05)
    s=np.array([500.,0.,0.,0.,1.05,0.])
    obj=trajectory(reason='NEAR_MISS', event=EventRecord('clearance_minimum',1.,s))
    ass=screen_trajectory(obj,disk,cfg)
    assert 'MINIMUM_GAP_NEAR_NEAR_MISS_THRESHOLD' in ass.flags


def test_recheck_can_warn_about_different_termination_classification():
    cfg=Config(initial_radius=500,final_radius=500,numerical_check='recheck')
    orig=trajectory(vz=0.)
    ref=trajectory(reason='MAX_TIME_REACHED',vz=0.)
    ass=assess_trajectory(orig,FlatDisk(),cfg,refined_trajectory=ref)
    assert ass.status=='RECHECK_DISAGREES'
    assert ass.recheck_reason=='MAX_TIME_REACHED'


def test_recheck_failure_is_recorded_not_ignored():
    cfg=Config(initial_radius=500,final_radius=500,numerical_check='recheck')
    ass=assess_trajectory(trajectory(),FlatDisk(),cfg,recheck_error=RuntimeError('no progress'))
    assert ass.status=='RECHECK_FAILED'
    assert 'no progress' in ass.note


def test_stream_sweep_writes_every_result_and_keeps_no_arrays(tmp_path,monkeypatch):
    cfg=Config(initial_radius=500,final_radius=510,radius_step=5.,numerical_check='off')
    sc=type('S',(),dict(velocity_kms_factor=100.))()
    sim=Simulation.__new__(Simulation)
    sim.cfg,sim.scales,sim.disk=cfg,sc,FlatDisk()
    # This test stubs the dynamics, so stub the new independent dust screening too.
    sim.radiation=type('Dust',(),{'emission_power':2.,
        'evaluate':lambda self,s:(np.zeros(3),1.)})()
    def fake(self,rho,launch_screen=None):
        if rho==505: raise ValueError('controlled cloud failure')
        return trajectory(radius=rho,vz=-.1)
    monkeypatch.setattr(Simulation,'trajectory',fake)
    updates=[]
    result=sim.run(tmp_path,on_error='continue',collect_results=False,
                   progress=lambda *x:updates.append(x))
    assert result==[]
    with (tmp_path/'summary.csv').open() as f: rows=list(csv.DictReader(f))
    assert len(rows)==3
    assert [r['termination_reason'] for r in rows]==[
        'DISK_IMPACT','NUMERICAL_FAILURE','DISK_IMPACT']
    assert all('max_z_rg_sampled' in r for r in rows)
    assert rows[0]['trajectory_file']
    assert rows[1]['error']=='controlled cloud failure'
    report=json.loads((tmp_path/'run_report.json').read_text())
    assert report['attempted_clouds']==3
    assert report['termination_counts']['NUMERICAL_FAILURE']==1
    assert len(updates)==3
    assert len(list((tmp_path/'Trajectory').rglob('*.dat')))==2
