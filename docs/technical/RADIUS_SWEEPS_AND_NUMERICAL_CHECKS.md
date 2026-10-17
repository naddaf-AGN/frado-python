# Launch-radius sweeps and numerical reliability

## Goal

A user specifies one physical model and a **range of cloud-launching radii**. Each radius produces an individual trajectory file (unless disabled), a row in `summary.csv`, and zero or more rows in `events.csv`. No plotting is performed by the physics engine. Launch radii use the disk's `Rg` units, **not metres**. The endpoints are inclusive, and the step must divide the interval to within floating-point precision.

## Example (3 clouds)

```bash
python run.py --radius-start 8370 --radius-end 8390 --radius-step 10 \
  --bh-mass-log 7 --mdot-log -0.4 --metallicity 5 --patch-size 3 \
  --max-orbits 2 --samples 2001 --numerical-check screen \
  --output-dir OutPut/example_sweep
python examples/plot_sweep_summary.py OutPut/example_sweep
```

For **one cloud**, substitute `--radius 8370`. For a large job, use `--dry-run` to confirm grid size before executing. `examples/reference_sweep.json` shows the equivalent reproducible JSON configuration, invoked by `--config examples/reference_sweep.json`.

On Google Colab, use `!python run.py ...` and `!python examples/plot_sweep_summary.py ...` in cells after extracting the ZIP and changing to its folder.

### Python API

```python
from frado import Config, Simulation
cfg = Config(bh_mass_log=7, mdot_log=-0.4, metallicity=5,
             initial_radius=8370, final_radius=8390, radius_step=10,
             maximum_orbits=2, output_samples=2001,
             numerical_check='screen', termination_mode='disk_impact')
sim = Simulation(cfg)
sim.run('OutPut/my_sweep', on_error='continue', collect_results=False)
```

`collect_results=False` writes the data incrementally without storing every cloud trajectory in RAM. For `collect_results=True`, the returned list contains `(radius, trajectory, path)` tuples. The CLI always uses streaming.

## Outputs

```
OutPut/my_sweep/
    metadata.json       # full input parameters, launch radii, unit definitions
    run_report.json     # attempted clouds + counts by termination and numerical status
    summary.csv         # one row per attempted cloud, including failures
    events.csv          # one row per detected event per cloud
    Trajectory/
        <model-label>/
            <model-label>_R8370.dat
            <model-label>_R8380.dat
            <model-label>_R8390.dat
    sweep_summary.png   # generated only when plotting script is explicitly used
```

Trajectory column meanings are in `metadata.json` and the header of each `.dat` file: time (years), `x,y,z` (`Rg`), `vx,vy,vz` (km/s), `sublimated_0_or_1`. Dust sublimation is an event: status 0 before and 1 after. The model saves all requested trajectories, not only long-lived ones.

`summary.csv` includes initial height, orbital and flight times, termination reason, final coordinates and velocities, sublimation time/location, whether escape speed was crossed, first vertical turning time, sampled extrema, and a numerical warning label. Extrema explicitly ending in `_sampled` are **sample-grid estimates**, not exact continuous maxima. `events.csv` holds all detected turning/clearance/sublimation/escape/impact events, including events on multiple orbital excursions. Each event row records its own state and time.

**IMPORTANT:** an event record is not necessarily a reason for termination; e.g., `escape_speed` is a diagnostic crossing. `MAX_TIME_REACHED` does not mean the cloud physically escaped or returned. `BLR_BOUNDARY` denotes crossing the chosen BLR radius.

## Numerical screening around grazing impacts

Adaptive ODE integration can shift apparent contact times substantially near a tangent touch; a very shallow impact and a genuine close near miss may be impossible to distinguish at the selected tolerance. The software supports **configurable numerical monitoring**, without altering the physical force law or termination criteria.

- `numerical_check='off'`: no screen.
- `numerical_check='screen'` (**default**): compute inexpensive heuristic warning indicators.
- `numerical_check='recheck'`: when a warning is raised, integrate that cloud again with smaller ODE tolerances and maximum solver step; compare outcome and flight time. Only flagged clouds are repeated.

Screen conditions:

1. For a local clearance minimum, relative gap `g/H(rho)` is within `grazing_screen_fraction` of zero, or (in first-excursion mode) within the same distance of the configured near-miss threshold. Here `g=z-H(rho)`.
2. For an actual impact, the dimensionless crossing rate `|dg/dt| * orbital_period / H(rho)` is below `grazing_speed_threshold`. The derivative includes the disk radial gradient and cloud radial velocity, so it detects shallow **surface** encounters rather than simply small vertical velocities.

Defaults: `grazing_screen_fraction=0.001`, `grazing_speed_threshold=0.3`; these are transparent **heuristics**, not physical constants or proven error bounds.

Status meanings:

| Status | Interpretation |
|---|---|
| `NOT_CHECKED` | Screen disabled. |
| `SCREEN_CLEAR` | No configured warning criterion triggered; this is **not** convergence proof. |
| `SCREEN_FLAGGED` | An ambiguous near-surface feature requires a tighter numerical check. |
| `RECHECK_AGREES` | A second run has the same termination reason, with flight-time difference below `recheck_time_tolerance_orbits` (default 0.005 orbital periods). |
| `RECHECK_DISAGREES` | Recheck changed termination reason or flight time beyond threshold. Inspect before scientific interpretation. |
| `RECHECK_FAILED` | The refined integration failed; original result remains saved but unverified. |

Default refinement: `relative_tolerance` and `absolute_tolerance` multiplied by `0.1`; `max_step_orbit_fraction` multiplied by `0.5`. **Quadrature settings remain unchanged during the targeted ODE recheck**, so this test does not bound radiation-integral error. For high-stakes scientific use, separately refine quadrature, solver method and maximum step, and check for convergence over multiple runs.

### Numerical example: a shallow impact

For the reference physical parameters and launch radius **35,750 Rg**, the baseline solver produced a trajectory termination time **60.471246711 yr** and flagged its shallow impact. A tighter ODE recheck also classified a disk impact but returned a terminal-time shift of **0.008027 orbital periods**, exceeding the default acceptance threshold: `RECHECK_DISAGREES`. This is **exactly** the ambiguity our screen is intended to expose. See `validation/extended/REPORT.md` for tolerance sweeps and a detailed discussion.

The recheck result does **not** silently replace the baseline outcome; both the decision and numerical warning are recorded separately.

## Safety and reproducibility

- Use **different output directories for different runs**. `summary.csv`, `events.csv`, `metadata.json` and `run_report.json` are recreated when a run starts. The code does not protect previous runs against reuse of an output path.
- The default `output_samples=500000` in the Python API can generate extremely large trajectory files. Explicitly choose `output_samples=1000–5000` (or a scientifically justified resolution) before running a sweep. Sampling affects saved curve resolution and sampled extrema, **not** the ODE solver tolerances.
- MATLAB filenames encode some parameters to only 1 decimal place; `metadata.json` is the authoritative source for exact input values.
- Large launch-radius grids can be computationally expensive. Choose the optional Numba backend to accelerate radiation integrals; automatic parallel execution is not supported.
- The reference comparisons in `validation/` involve selected trajectories; they do not validate all 2,803 supplied MATLAB event records.
