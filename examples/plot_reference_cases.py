"""Plot three precomputed Python and reference MATLAB trajectories (no solver run).

From package root: python examples/plot_reference_cases.py
Result: figures/validation_three_trajectories.png
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib.pyplot as plt
from frado import Config, Simulation, plot_trajectory
from frado.integrator import Trajectory

sim=Simulation(Config(bh_mass_log=7.,mdot_log=-.4,metallicity=5.,patch_size=3.))
fig, axes=plt.subplots(1,3,figsize=(18,5))
for ax,r in zip(axes,(8370,22990,30860)):
    cache=Path(f'validation/reference_python_R{r}.npz')
    with np.load(cache) as data:
        sub_state=data['sublimation_state']
        sub_state=None if not np.all(np.isfinite(sub_state)) else sub_state
        traj=Trajectory(time_s=data['time_s'],states=data['states'],dusty=data['dusty'],
                        events=[],sublimation_s=float(data['sublimation_s']),
                        sublimation_state=sub_state,
                        orbital_period_s=float(data['orbital_period_s']),
                        termination_reason=str(data['termination_reason']),
                        event_records=[],initial_radius=float(r))
    ref=Path(f'validation/reference/T1500_A3_Z5_M7.0_m-0.4_R{r}.dat')
    plot_trajectory(traj,disk=sim.disk,matlab_path=ref,ax=ax)
Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/validation_three_trajectories.png',dpi=150,bbox_inches='tight')
print('Saved figures/validation_three_trajectories.png')
