# SciPy and Numba radiation integration backends

## Available backends

The `scipy` backend is the default. The optional `numba` backend compiles integrands to accelerate SciPy's adaptive integration; each run's `metadata.json` records the selected backend.

## Implementation

The `scipy` implementation of `Radiation.evaluate` integrates the cylindrical radiation force, vertical radiation force, and absorbed dust power using three independent SciPy `dblquad` calls. For each integrand sample, a Python callback evaluates seven-term analytical opacity fits, a geometrical denominator, and the shielding patch limits.

The `numba` backend uses a **Numba-compiled low-level callback** called directly from SciPy QUADPACK through `LowLevelCallable`, avoiding most Python-per-sample overhead. SciPy `nquad` uses the same nested R/phi adaptive quadrature order as `dblquad` (inner phi, outer R), and uses the *same* quadrature tolerances. The seven-term total/absorption opacity fits, force coefficients, shielding, integrator, disk geometry, and user-facing stopping conditions are unchanged.

`Radiation.absorbed_power(state)`: the compiled backend integrates only the absorption kernel when screening dusty launches or determining sublimation. This saves two unnecessary force integrals for those checks. In the SciPy backend, the full `evaluate` path is used. Third-party radiation implementations without `absorbed_power` can use `evaluate`.

The backend comparison does not use an interpolation grid, modified opacity fit, or relaxed quadrature tolerance. The `numba` backend is optional. Installing `requirements-optimized.txt` enables the additional dependency; running the SciPy backend does not require Numba.

## User interface

```bash
python run.py --config examples/optimized_sweep.json --output-dir OutPut/optimized_01
# Or override the backend for any existing JSON config:
python run.py --config examples/reference_sweep.json --radiation-backend numba --output-dir OutPut/optimized_02
# Select the SciPy backend explicitly:
python run.py --config examples/reference_sweep.json --radiation-backend scipy --output-dir OutPut/reference_02
```

Python:

```python
from frado import Config, Simulation
cfg = Config(bh_mass_log=7.0, mdot_log=-0.4,
             initial_radius=8370, final_radius=8390, radius_step=10,
             radiation_backend='numba')
sim = Simulation(cfg)
sim.run('OutPut/example', on_error='continue', collect_results=False)
```

For Colab, install dependencies with `!pip -q install -r requirements-optimized.txt` before running.

## Selected benchmark measurements

The stored, runnable benchmarks are `validation/benchmark_radiation.py`, `validation/benchmark_trajectory.py`, and `validation/benchmark_reference_clouds.py`. Their machine-readable measurements are under `validation/optimization/`. Timing values are **machine-, CPU-load-, SciPy-, and Numba-version-dependent**, and initial Numba compilation was excluded from the steady-state microbenchmark.

| Test | Reference SciPy | Compiled Numba | Speedup |
| --- | ---: | ---: | ---: |
| Radiation evaluations at six representative states (median of per-state ratios) | varies | varies | ~36× |
| Full trajectory, R=8,370 Rg | 12.60 s | 0.31 s | ~41× |
| Full trajectory, R=22,990 Rg | 18.94 s | 0.63 s | ~30× |

SciPy-vs-Numba **same-physics** trajectories agreed to ~1.0×10⁻⁶ Rg maximum position difference for R=8,370, ~4.2×10⁻⁷ Rg for R=22,990; termination classifications matched, and sublimation time differed by ~0.0003 s for the former. The six quadrature comparisons were within 6.4×10⁻¹² relative difference at the tested points. These are *not* universal guarantees for all parameter choices, and do *not* eliminate numerical uncertainty in grazing impacts.

The complete R=30,860 trajectory was not included in the full backend-to-backend benchmark: its reference solve was unusually slow during the test window. Its **initial and selected raised-height force calculations** passed the quadrature comparison. More challenging grazing and long-duration runs still require dedicated convergence checks.

**43 automated tests passed**, including regression comparisons of heating-only integration, launch-screening status, and radiation at multiple states. The optional optimized demonstration computed three trajectories (R=8,370, 8,380, 8,390 Rg) and saved `summary.csv`, `events.csv`, `launch_screening.csv`, and the trajectory files.

## Reproduction

```bash
python -m pytest -q
python validation/benchmark_radiation.py
python validation/benchmark_trajectory.py
python validation/benchmark_reference_clouds.py
python run.py --config examples/optimized_sweep.json --output-dir OutPut/demo
```

**Limits:** one CPU / interpreter benchmark, with a moderate sample of six force positions and two fully solved clouds; not a comprehensive parameter-grid validation, nor a formal comparison to original MATLAB across the entire launch-radius range.
