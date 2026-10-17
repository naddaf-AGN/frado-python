"""Smoke tests for per-run output, event history and optional visualization."""
import csv
import json
import numpy as np
from frado import Config, plot_trajectory
from frado.constants import YEAR
from frado.integrator import Trajectory,EventRecord
from frado.outputs import OutputWriter, save_trajectory


def fake_result():
    a=np.array([500.,0.,0.,0.00002,4.,0.])
    b=np.array([501.,0.,0.,0.00001,3.9,-0.00001])
    c=np.array([502.,0.,0.,0.000005,3.8,-0.00002])
    return Trajectory(
        np.array([0.,YEAR,2*YEAR]), np.stack((a,b,c)),
        np.array([True,False,False]),
        [(6,YEAR,b),(2,2*YEAR,c)], YEAR, b,
        orbital_period_s=10*YEAR,termination_reason='DISK_IMPACT',
        event_records=[EventRecord('sublimation',YEAR,b),EventRecord('disk_impact',2*YEAR,c)],
        initial_radius=500.)


def test_csv_records_status_and_eight_column_trajectory(tmp_path):
    c=Config(initial_radius=500,final_radius=500)
    sc=type('S',(),dict(velocity_kms_factor=100.))()
    t=fake_result()
    writer=OutputWriter(tmp_path,c)
    f=save_trajectory(t,c,sc,500,tmp_path)
    writer.add(t,c,sc)
    data=np.loadtxt(f)
    assert data.shape==(3,8)
    assert np.allclose(data[:,7],[0,1,1])
    with (tmp_path/'summary.csv').open() as h: rows=list(csv.DictReader(h))
    assert rows[0]['termination_reason']=='DISK_IMPACT'
    assert float(rows[0]['sublimation_time_years'])==1.
    with (tmp_path/'events.csv').open() as h: events=list(csv.DictReader(h))
    assert [e['event_name'] for e in events]==['sublimation','disk_impact']
    assert [e['dusty_at_event'] for e in events]==['0','0']
    metadata=json.loads((tmp_path/'metadata.json').read_text())
    assert metadata['configuration']['termination_mode']=='disk_impact'


def test_failure_logging(tmp_path):
    c=Config()
    writer=OutputWriter(tmp_path,c)
    writer.add_failure(999.,c,RuntimeError('test failure'))
    with (tmp_path/'summary.csv').open() as h: row=next(csv.DictReader(h))
    assert row['termination_reason']=='NUMERICAL_FAILURE'
    assert 'test failure' in row['error']


def test_optional_plotting(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    t=fake_result()
    fig,ax=plot_trajectory(t,show_disk=False)
    labels=[line.get_label() for line in ax.lines]
    assert 'Python dusty' in labels and 'Python dust-free' in labels
    fig.savefig(tmp_path/'test.png')
    plt.close(fig)
    assert (tmp_path/'test.png').stat().st_size>100
