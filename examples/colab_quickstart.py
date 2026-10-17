"""Run from package root: python examples/colab_quickstart.py

Produces an R=8370 trajectory and two easy-to-review files:
OutPut/summary.csv, figures/trajectory_8370.png.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from frado import Config, Simulation, plot_trajectory

cfg=Config(bh_mass_log=7.0,mdot_log=-0.4,metallicity=5.0,patch_size=3,
           initial_radius=8370,final_radius=8370,
           termination_mode='disk_impact', maximum_orbits=10,output_samples=5000)
sim=Simulation(cfg)
launch,result,path=sim.run(output_dir='OutPut')[0]
print(f'R={launch:g} Rg: {result.termination_reason}; time={result.time_s[-1]/31536000:.4f} yr')
print('Trajectory saved:',path)
reference=Path('validation/reference/T1500_A3_Z5_M7.0_m-0.4_R8370.dat')
fig,ax=plot_trajectory(result,disk=sim.disk,matlab_path=reference if reference.exists() else None)
Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/trajectory_8370.png',dpi=170,bbox_inches='tight')
print('Comparison plot: figures/trajectory_8370.png')
