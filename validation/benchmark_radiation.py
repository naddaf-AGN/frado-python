"""Reproducible backend cross-check and timed radiation microbenchmark.

Run from project root: python validation/benchmark_radiation.py
Results live in validation/optimization/ (no modification of baseline).
"""
import json
import statistics
import time
from dataclasses import replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from frado import Config, Simulation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'validation' / 'optimization'
OUT.mkdir(exist_ok=True)


def one_backend(config, backend, states, trials=3):
    sim = Simulation(replace(config, radiation_backend=backend))
    # Explicit warm-up includes JIT compilation and SciPy initialization,
    # so the timed figures represent steady-state simulation calls only.
    sim.radiation.evaluate(states[0])
    measurements = []
    for state in states:
        timings = []
        values = []
        for _ in range(trials):
            start = time.perf_counter()
            result = sim.radiation.evaluate(state)
            timings.append(time.perf_counter()-start)
            values.append(np.r_[result[0], result[1]])
        measurements.append({'outputs': values[0].tolist(),
                             'seconds_median': statistics.median(timings),
                             'seconds_runs': timings})
    return measurements


def main():
    cfg=Config(bh_mass_log=7., mdot_log=-0.4, metallicity=5.,
               initial_radius=8370., final_radius=8370.)
    sim=Simulation(cfg)
    # Disk surface/above-surface launch states, plus representative high-z tests.
    points = [(5000., None), (8370., None), (22990., None), (30860., None),
              (8370., 100.), (22990., 500.)]
    states = []
    for rho,z in points:
        if z is None:
            z=float(sim.disk.at(rho))
        states.append(np.array([rho,0,0,0,z,0],dtype=float))
    baseline=one_backend(cfg, 'scipy', states)
    accelerated=one_backend(cfg, 'numba', states)
    checks=[]
    for i, (reference, optimized) in enumerate(zip(baseline, accelerated)):
        a=np.array(reference['outputs'])
        b=np.array(optimized['outputs'])
        abs_difference = np.abs(a-b)
        relative_difference = abs_difference / np.maximum(np.abs(a), 1e-300)
        passed = bool(np.allclose(a, b, rtol=1e-8, atol=1e-15))
        checks.append({'radius':points[i][0], 'height':float(states[i][4]),
                       'baseline':reference, 'optimized':optimized,
                       'max_relative_difference':float(np.max(relative_difference)),
                       'max_absolute_difference':float(np.max(abs_difference)),
                       'speedup':reference['seconds_median']/optimized['seconds_median'],
                       'passed':passed})
    report={'kind':'steady-state callback microbenchmark', 'trials_per_state':3,
            'comparison_rtol':1e-8,'comparison_atol':1e-15,
            'results':checks,'passed':bool(all(c['passed'] for c in checks)),
            'median_speedup':float(statistics.median(c['speedup'] for c in checks))}
    (OUT/'benchmark_results.json').write_text(json.dumps(report,indent=2))
    print('Backend comparison PASS:',report['passed'])
    print('Median per-call speedup:',round(report['median_speedup'],2),'x')
    for c in checks:
        print('R=', c['radius'],'z=',round(c['height'],3),
              'speedup=',round(c['speedup'],2),'x',
              'max rel diff=',format(c['max_relative_difference'],'.3e'),
              'passed=',c['passed'])
    if not report['passed']:
        raise AssertionError('Optimized radiation integrals diverged from reference')

if __name__=='__main__':main()
