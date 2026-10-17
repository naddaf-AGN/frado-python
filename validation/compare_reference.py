"""Compare FRADO Python with three supplied MATLAB reference trajectories.

Run from project root: python validation/compare_reference.py
To recompute a single cache: python validation/compare_reference.py --solve 8370
Cached runs are saved in validation/reference_python_R*.npz.
"""
from pathlib import Path
import sys
import json
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from frado import Config,Simulation  # noqa: E402
from frado.constants import YEAR  # noqa: E402

ROOT_REF=ROOT/'validation'/'reference'
HERE=ROOT/'validation'


def _solve_one(radius):
    # Independent process avoids state coupling between adaptive QUADPACK/ODE
    # computations and makes individual wall times reproducible.
    cfg=Config(bh_mass_log=7.,mdot_log=-.4,metallicity=5.,patch_size=3.,
        maximum_orbits=2,output_samples=2001,termination_mode='disk_impact')
    sim=Simulation(cfg)
    start=time.monotonic()
    result=sim.trajectory(radius)
    elapsed=time.monotonic()-start
    np.savez_compressed(HERE/f'reference_python_R{radius}.npz',
        time_s=result.time_s,states=result.states,dusty=result.dusty,
        sublimation_s=result.sublimation_s,
        sublimation_state=(result.sublimation_state if result.sublimation_state is not None else np.full(6,np.nan)),
        termination_reason=result.termination_reason,orbital_period_s=result.orbital_period_s,
        runtime_seconds=elapsed)
    print(f'Computed R={radius} in {elapsed:.2f} s: {result.termination_reason}',flush=True)


def validate():
    cfg=Config(bh_mass_log=7.,mdot_log=-.4,metallicity=5.,patch_size=3.,
        maximum_orbits=2,output_samples=2001,termination_mode='disk_impact')
    sim=Simulation(cfg)
    all_results=[]
    mat_events=np.loadtxt(ROOT_REF/'event_T1500_A3_Z5_M7.0_m-0.4.dat')
    for radius in [8370,22990,30860]:
        ref=np.loadtxt(ROOT_REF/f'T1500_A3_Z5_M7.0_m-0.4_R{radius}.dat')
        cache=HERE/f'reference_python_R{radius}.npz'
        if not cache.exists():
            raise FileNotFoundError(
                f'{cache} missing. Run `python validation/compare_reference.py --solve {radius}` '
                'separately, then rerun this report.')
        with np.load(cache) as f:
            t=f['time_s']/YEAR
            p=f['states']
            py_sub=float(f['sublimation_s'])/YEAR
            reason=str(f['termination_reason'])
            seconds=float(f['runtime_seconds'])
        # Dusty-dynamics comparison *before* the earliest sublimation, if any.
        # MATLAB global dust flag can flip within a rejected adaptive solver step.
        evt=mat_events[mat_events[:,0]==radius][0]
        mat_sub=evt[4]
        cutoff=min(ref[-1,0],t[-1],mat_sub,py_sub)
        if np.isfinite(py_sub): cutoff=max(0,cutoff-.05)
        ref_slice=ref[ref[:,0]<=cutoff]
        vals=np.column_stack((p[:,0],p[:,2],p[:,4],
            p[:,1]*sim.scales.velocity_kms_factor,
            p[:,3]*sim.scales.velocity_kms_factor,
            p[:,5]*sim.scales.velocity_kms_factor))
        py_mat=np.column_stack([np.interp(ref_slice[:,0],t,vals[:,i]) for i in range(6)])
        d=py_mat-ref_slice[:,1:7]
        poserr=np.linalg.norm(d[:,:3],axis=1)
        velerr=np.linalg.norm(d[:,3:],axis=1)
        endrho=np.hypot(p[-1,0],p[-1,2]); h_end=float(sim.disk.at(endrho))
        one=dict(launch_radius_rg=radius,
                 matlab_end_years=float(ref[-1,0]),
                 python_end_years=float(t[-1]),
                 python_termination_reason=reason,
                 python_final_clearance_rg=float(p[-1,4]-h_end),
                 matlab_sublimation_years=(float(mat_sub) if np.isfinite(mat_sub) else None),
                 python_sublimation_years=(float(py_sub) if np.isfinite(py_sub) else None),
                 dusty_max_position_difference_rg=float(poserr.max()),
                 dusty_max_velocity_difference_kms=float(velerr.max()),
                 dusty_comparison_samples=int(len(ref_slice)),
                 runtime_seconds=round(seconds,2))
        all_results.append(one)
        print('R',radius,one,flush=True)
    (HERE/'reference_comparison.json').write_text(json.dumps(all_results,indent=2)+'\n')
    rows=['# FRADO Python compared with supplied MATLAB reference trajectories','',
          'Disk contact follows a downward crossing of z=H(rho). The MATLAB reference uses different event semantics.',
          'Final times need not match when the event or sublimation semantics differ.','',
          '| Launch [Rg] | MATLAB end [yr] | Python end [yr] | Python status | Python final z-H [Rg] | Max dusty Δpos [Rg] | Max dusty Δvel [km/s] |',
          '|---:|---:|---:|---|---:|---:|---:|']
    for v in all_results:
        rows.append(f'| {v["launch_radius_rg"]:,} | {v["matlab_end_years"]:.6f} | '
                    f'{v["python_end_years"]:.6f} | {v["python_termination_reason"]} | '
                    f'{v["python_final_clearance_rg"]:.8g} | '
                    f'{v["dusty_max_position_difference_rg"]:.4f} | '
                    f'{v["dusty_max_velocity_difference_kms"]:.4f} |')
    rows += ['','These are the same three supplied cases, **not** a full re-run of the 2,803-cloud MATLAB grid.',
             'The pre-sublimation comparison is intentionally limited; different state-transition handling',
             'in MATLAB/Python prevents interpreting full-flight discrepancies as purely integration error.',
             'A suite of artificial analytic trajectories independently tests disk impact, near misses,',
             'physical grazing contact and continuation into additional orbits.',
             'MATLAB event IDs 2 and 4 do not have identical semantics in the Python stopping policies.',
             '']
    (HERE/'REFERENCE_COMPARISON.md').write_text('\n'.join(rows))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--solve':
        _solve_one(int(sys.argv[2]))
    else:
        validate()
