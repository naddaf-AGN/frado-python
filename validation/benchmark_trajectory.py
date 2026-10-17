"""Full FRADO trajectory benchmark of reference and accelerated quadrature.

Run from root: python validation/benchmark_trajectory.py
No reliance on MATLAB for equivalence; both solvers use the same FRADO model.
"""
import json
import sys
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from frado import Config, Simulation
from frado.constants import YEAR

OUT=Path(__file__).resolve().parent/'optimization'
OUT.mkdir(exist_ok=True)


def main():
    cfg=Config(initial_radius=8370,final_radius=8370,
               bh_mass_log=7.0, mdot_log=-0.4, metallicity=5.,
               patch_size=3.,maximum_orbits=1.0, output_samples=1001,
               numerical_check='off')
    outputs={}
    for backend in ('scipy','numba'):
        sim=Simulation(replace(cfg, radiation_backend=backend))
        # Do not count first JIT import/compile in the run cost; document it.
        start=time.perf_counter()
        result=sim.trajectory(8370)
        elapsed=time.perf_counter()-start
        outputs[backend]={'result':result, 'seconds':elapsed}
        print(backend,round(elapsed,3),'s',result.termination_reason,
              'flight:',result.time_s[-1]/YEAR,'years',
              'sublimation:',result.sublimation_s/YEAR,'years')
    old=outputs['scipy']['result']; new=outputs['numba']['result']
    # Grids can differ because event termination times differ slightly.
    common_time=np.linspace(0,min(old.time_s[-1],new.time_s[-1]),501)
    a=np.column_stack([np.interp(common_time,old.time_s,old.states[:,j]) for j in range(6)])
    b=np.column_stack([np.interp(common_time,new.time_s,new.states[:,j]) for j in range(6)])
    difference=np.abs(a-b)
    position_max=float(np.max(difference[:,[0,2,4]]))
    velocity_max=float(np.max(difference[:,[1,3,5]]))
    checks={
       'termination_matches':bool(old.termination_reason==new.termination_reason),
       'sublimation_time_abs_seconds':float(abs(old.sublimation_s-new.sublimation_s)),
       'flight_time_abs_seconds':float(abs(old.time_s[-1]-new.time_s[-1])),
       'max_position_difference_rg':position_max,
       'max_velocity_difference_rg_per_second':velocity_max,
       'speedup':outputs['scipy']['seconds']/outputs['numba']['seconds'],
       'times_seconds':{k:v['seconds'] for k,v in outputs.items()},
       'backend_notes':'both solve_ivp BDF; numba uses same adaptive QUADPACK with compiled integrands',
    }
    (OUT/'trajectory_results.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))
    assert checks['termination_matches']
    assert position_max < 1.0, 'Trajectory deviation unexpectedly exceeds 1 Rg'

if __name__=='__main__':main()
