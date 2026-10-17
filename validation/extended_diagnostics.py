"""Independent checks at stored FRADO sublimation states; no ODE rerun.

Run after validation/extended_convergence.py --suite full.
"""
from pathlib import Path
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from frado import Config,Simulation

HERE=ROOT/'validation'/'extended'
results=json.loads((HERE/'results.json').read_text())
records=[]
for row in results['real_runs']:
    if row['sublimation_years'] is None: continue
    r=row['radius_rg']; label=row['case']
    with np.load(HERE/f'R{r}_{label}.npz') as data:
        state=data['sublimation_state'];sub=float(data['sublimation_s'])
        times=data['time_s'];dusty=data['dusty']
    cfg=Config(bh_mass_log=7.,mdot_log=-.4,metallicity=5.,patch_size=3.,
               maximum_orbits=2.,quad_epsabs=row['quad_epsabs'],quad_epsrel=row['quad_epsrel'])
    sim=Simulation(cfg)
    radiation,absorbed=sim.radiation.evaluate(state)
    f=absorbed/sim.radiation.emission_power-1
    a_before=sim.dynamics.acceleration(state,True)
    a_after=sim.dynamics.acceleration(state,False)
    discontinuity=float(np.linalg.norm(a_before-a_after-radiation))
    mask_consistent=bool(np.array_equal(dusty,times<sub))
    records.append({'radius_rg':r,'case':label,
        'threshold_relative_residual':float(f),
        'acceleration_jump_residual_rg_s2':discontinuity,
        'dust_state_tags_consistent':mask_consistent,
        'sublimation_years':row['sublimation_years'],
        'sublimation_radius_rg':float(np.hypot(state[0],state[2]))})
(HERE/'sublimation_diagnostics.json').write_text(json.dumps(records,indent=2)+'\n')
for rec in records:
    print(rec)
