# FRADO Python 1.0.0 — technical reference

## Physics and numerical flow

1. `Config` defines the black-hole mass, Eddington ratio, dust-to-gas scaling (via metallicity), shielding, numerical tolerances and launching-radius grid.
2. The built-in fitted disk surface determines launch height `H(rho)` and disk-impact geometry.
3. The shipped `frado/k_average_tot.dat` supplies absorption opacity for computing the sublimation-temperature dust emission. Seven-term fitted opacities (total and absorption) are used in 2-D disk-integrated acceleration/heating, not interpolated from the table on every step.
4. At **each** requested radius, check absorbed vs emitted dust power. If heating >= emission, record `NO_DUST_AT_LAUNCH` and skip integrating that radius. Survival is necessary but not sufficient for launching.
5. Integrate the 3-D position/velocity system; detect and restart at sublimation, after which radiation pressure is disabled. Record event history and numerical quality flags.

Internal positions use `Rg=GM/c²`, time uses seconds, and velocities use `Rg/s`; exported times use years and velocities use km/s.

## Essential configuration keys

| Key | Meaning | Default |
| --- | --- | --- |
| `bh_mass_log` | log10(M_BH/M_sun), 6–10 | 8.5 |
| `mdot_log` | log10(Eddington ratio), −2–0 | −1.7 |
| `metallicity` | Dust-scaled metallicity (solar) | 5 |
| `sublimation_temperature` | Dust destruction threshold (K) | 1500 |
| `patch_size` | Patch shielding parameter | 3 |
| `initial_radius`, `final_radius`, `radius_step` | Launch grid (Rg), endpoint included only if step divides interval | 2340, 4000, 10 |
| `termination_mode` | `disk_impact` or `first_excursion` | `disk_impact` |
| `near_miss_fraction` | Maximum gap/H for `first_excursion` near miss | 0.05 |
| `maximum_orbits` | Maximum number of initial orbital periods | 5 |
| `output_samples` | Output samples per trajectory | 500000 |
| `radiation_backend` | `scipy` or `numba` | `scipy` |
| `numerical_check` | `off`, `screen` or `recheck` | `screen` |
| `relative_tolerance`, `absolute_tolerance` | ODE tolerances | 1e-7, 1e-6 |
| `max_step_orbit_fraction` | Maximum ODE step / orbital period | 0.02 |

All valid keys appear in `frado/config.py` and `python run.py --help`. The bundled optimized example overrides the default backend and output sample count.

**Termination meanings:** `DISK_IMPACT` = downward disk crossing or tangential contact; `NEAR_MISS` = first close clearance minimum (only in `first_excursion` mode); `BLR_BOUNDARY` and `MAX_TIME_REACHED` = other stopping limits. The `disk_impact` mode *can* continue into subsequent excursions after a miss. Grazing classifications can depend on tolerances; interpret numerical flags.

## Outputs

| File | Contents |
| --- | --- |
| `summary.csv` | one record per requested radius, including no-launch and failure states, launching/ending/sublimation quantities, reliability diagnostics |
| `launch_screening.csv` | initial radius, disk height, heating/emission powers and ratio, launching status |
| `events.csv` | time, position, velocity, event name/ID and dust status per event |
| `Trajectory/.../*.dat` | numeric columns: `time_years x_rg y_rg z_rg vx_kms vy_kms vz_kms sublimated_0_or_1` |
| `metadata.json` | complete run configuration, units and output column definitions |
| `run_report.json` | aggregate outcome and screening counts |

Use a **unique output directory per run**. Set `--on-error raise` for fail-fast debugging, `--on-error continue` to record individual failures, and `--skip-trajectories` for summaries without trajectory files.

## Extending the physics

`frado/shielding.py` encapsulates the patch prescription, while `frado/radiation.py` and `frado/radiation_numba.py` expose reference and compiled quadratures. `frado/simulation.py` accepts a shielding factory. Modified physics needs **independent scientific validation**—a faster computation alone is not proof of physical equivalence. Custom-opacity upload and multiprocessing are planned, not implemented here.
