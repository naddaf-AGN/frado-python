# Dust-survival screening at launch

## Rationale

A dusty cloud cannot be launched in the model from a disk surface where its grains have already sublimated. The launch screen excludes such cases while retaining every requested radius in the output tables.

## Criterion

For every requested cylindrical radius `rho_0`, compute the disk surface
`z_0=H(rho_0)` and evaluate the **same patch-shielded, fitted-opacity dust
absorption integral** used by the trajectory model at
`(x,y,z)=(rho_0,0,H(rho_0))`. Compare the absorbed power `P_abs(rho_0)` to
the configured dust thermal emission power `P_emit(T_sub)` computed from
`k_average_tot.dat`.

- `P_abs < P_emit`: `LAUNCH_ALLOWED` — **dust survival** at launch.
- `P_abs >= P_emit`: `NO_DUST_AT_LAUNCH` — the requested radius is skipped,
  including exact equality.
- Invalid or nonfinite powers: `SCREENING_FAILED`, not incorrectly classified
  as a physical no-launch; with `--on-error continue`, the other radii proceed.

**Interpretation:** Survival of dust is a *necessary*, but not a sufficient,
condition for successful dynamical lifting from the disk. This filter does NOT
compute a separate acceleration-based lift-off criterion. During integration,
dust can still sublimate later at a finite flight time.

The screening is performed for **each radius separately** rather than
assuming all allowed launching radii form one continuous interval. The launch
scan uses the actual user-specified spacing; the report's minimum and maximum
eligible radii are bounds of screened samples, NOT inferred physical thresholds.

## Outputs

- `launch_screening.csv`: one row for every requested radius, with physical
  disk height, status, absorbed power, emitted power, ratio, and error column.
- `summary.csv`: one row for every requested radius; rejected entries have
  `termination_reason=NO_DUST_AT_LAUNCH`, `launch_status=NO_DUST_AT_LAUNCH`,
  heating diagnostic, and blank flight/velocity/trajectory fields.
- `events.csv`: **no** fabricated sublimation event for skipped clouds.
- `Trajectory/`: only launched and successfully integrated cloud trajectories.
- `run_report.json`: counts of accepted/rejected/failed checks and extrema of
  accepted sampled launching radii.

A run where all radii are rejected is a valid run. It still writes all reports.

## Usage

```sh
python run.py --config examples/launch_screening_demo.json --output-dir OutPut/screen_test
```

The bundled demonstration uses `M=10^7 M_sun`, Eddington ratio `10^-0.4`,
`Z=5`, `T_sub=1500 K`, patch size 3, and two launch radii. With this model,
`5000 Rg` fails dust survival (power ratio ~5.925) whereas `10000 Rg` passes
(~0.358). The demonstration uses a deliberately short maximum integration
time for fast execution; its trajectory termination is **not** a science-case
orbital solution.

Single-radius checks from Python:

```python
from frado import Config, Simulation
cfg = Config(bh_mass_log=7, mdot_log=-0.4,
             initial_radius=5000, final_radius=10000, radius_step=5000)
sim = Simulation(cfg)
for rho in cfg.launch_radii:
    s = sim.screen_launch(rho)
    print(rho, s.status, s.heating_ratio)
```

A direct call to `sim.trajectory(rho)` on a rejected radius raises
`NoDustAtLaunchError` instead of manufacturing a gravitational-only cloud.

## Numerical interpretation

When `P_abs/P_emit` is very close to 1, quadrature accuracy can affect the
binary classification. The output preserves the exact ratio so boundary cases
can be reviewed and retested at tighter quadrature tolerances. This screen does **not** claim to locate the exact continuous inner dust-survival radius.
