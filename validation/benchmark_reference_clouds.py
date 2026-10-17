"""Run two supplied MATLAB reference radii with *both* Python backends.

This isolates performance/equivalence of the quadrature optimization;
MATLAB-to-Python physical differences were documented in previous releases.
"""
import json, sys, time
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from frado import Config,Simulation
from frado.constants import YEAR


def one(radius):
    cfg=Config(initial_radius=radius, final_radius=radius,
               bh_mass_log=7.0,mdot_log=-.4,metallicity=5.,patch_size=3.,
               maximum_orbits=1.25, output_samples=1001, numerical_check='off')
    runs={}
    for backend in ('scipy','numba'):
        sim=Simulation(replace(cfg, radiation_backend=backend))
        t=time.perf_counter()
        traj=sim.trajectory(radius)
        runs[backend]=(traj,time.perf_counter()-t)
    old, old_seconds=runs['scipy']
    new, new_seconds=runs['numba']
    end=min(old.time_s[-1],new.time_s[-1])
    sample=np.linspace(0,end,501)
    a=np.column_stack([np.interp(sample,old.time_s,old.states[:,j]) for j in range(6)])
    b=np.column_stack([np.interp(sample,new.time_s,new.states[:,j]) for j in range(6)])
    return {
        'radius_rg':radius,
        'termination_reference':old.termination_reason,
        'termination_optimized':new.termination_reason,
        'flight_time_reference_years':old.time_s[-1]/YEAR,
        'flight_time_optimized_years':new.time_s[-1]/YEAR,
        'sublimation_reference_years':None if not np.isfinite(old.sublimation_s) else old.sublimation_s/YEAR,
        'sublimation_optimized_years':None if not np.isfinite(new.sublimation_s) else new.sublimation_s/YEAR,
        'max_position_difference_rg':float(np.max(np.abs(a[:,[0,2,4]]-b[:,[0,2,4]]))),
        'reference_seconds':old_seconds,
        'optimized_seconds':new_seconds,
        'speedup':old_seconds/new_seconds,
    }


def main():
    results=[]
    for radius in (8370.,22990.):
        result=one(radius)
        results.append(result)
        print(f"R {radius:g}: {result['reference_seconds']:.3f} -> {result['optimized_seconds']:.3f} s, "
              f"speedup={result['speedup']:.1f}x, Δr_max={result['max_position_difference_rg']:.3g} Rg",flush=True)
    report={'cases':results,'all_termination_matches':all(r['termination_reference']==r['termination_optimized'] for r in results)}
    out=Path(__file__).parent/'optimization'/'two_cloud_results.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2))
    assert report['all_termination_matches']
    assert all(r['max_position_difference_rg']<.01 for r in results)

if __name__=='__main__':main()
