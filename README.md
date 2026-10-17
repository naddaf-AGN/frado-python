# FRADO Python

**Version 1.0.0** · [Scientific validation](docs/VALIDATION.md) · [Google Colab quickstart](docs/QUICKSTART.md)

FRADO Python is a scientific implementation of the **2.5D Failed Radiatively Accelerated Dusty Outflow (FRADO) model**, simulating three-dimensional trajectories of dusty clouds above accretion disks in active galactic nuclei (AGN). It provides dust heating/sublimation, disk radiation pressure, configurable trajectory termination, launch-radius sweeps, and optional Numba acceleration of radiation integrals.

**Scope:** This implementation has numerical regression tests and selected comparisons with reference calculations. It has **not** been validated across the entire AGN parameter space. See [Validation and limitations](docs/VALIDATION.md) before interpreting scientific results.

## Run and test

Use Python 3.11 or a compatible Python environment. Install required packages:

```bash
python -m pip install -r requirements-optimized.txt
python -m pytest -q
```

Run a short demonstrative sweep (Numba backend):

```bash
python run.py --config examples/optimized_sweep.json --dry-run
python run.py --config examples/optimized_sweep.json --output-dir OutPut/demo
python examples/plot_sweep_summary.py OutPut/demo
```

The example solves three launch radii (8370, 8380, 8390 gravitational radii). Initial Numba compilation may add startup time. For a single cloud:

```bash
python run.py --radius 8370 --bh-mass-log 7 --mdot-log -0.4 \
  --metallicity 5 --patch-size 3 --max-orbits 1.25 \
  --radiation-backend numba --output-dir OutPut/one_cloud
```

**No local Python installation is necessary:** the [Google Colab instructions](docs/QUICKSTART.md) run the code in an online notebook.

## Parameters and numerical methods

Example configuration (`examples/optimized_sweep.json`):

```json
{
  "bh_mass_log": 7.0,
  "mdot_log": -0.4,
  "metallicity": 5.0,
  "patch_size": 3.0,
  "initial_radius": 8370.0,
  "final_radius": 8390.0,
  "radius_step": 10.0,
  "radiation_backend": "numba",
  "termination_mode": "disk_impact",
  "numerical_check": "screen",
  "maximum_orbits": 1.25,
  "output_samples": 1001
}
```

- The `scipy` radiation backend is the default; `numba` uses compiled callbacks with the same fitted opacities and adaptive quadrature.
- `disk_impact` follows a trajectory to disk contact or another stopping boundary; `first_excursion` can also stop at a configured near-miss clearance minimum.
- A launch-radius sample is skipped when the dust cannot survive at its disk surface (`NO_DUST_AT_LAUNCH`). Dust survival alone **does not** prove successful radiative launching.
- The fitted disk prescriptions require `6 ≤ bh_mass_log ≤ 10` and `−2 ≤ mdot_log ≤ 0`; supported launch radii span 100–100,000 `Rg`. These input bounds do not imply uniform numerical or physical validity across the range.
- The built-in dust opacity file is used as supplied. User-defined opacity input and parallelized sweeps are not supported in this version.

See the [parameter and output reference](docs/REFERENCE.md) and the [numerical-method notes](docs/technical/) for details.

## Output

Each run writes `summary.csv`, `events.csv`, `launch_screening.csv`, `metadata.json`, and `run_report.json`. Optional trajectory files contain eight columns:

`time [yr], x [Rg], y [Rg], z [Rg], vx [km/s], vy [km/s], vz [km/s], sublimated [0/1]`.

Use a fresh output directory for each run; files in an existing output directory may be overwritten.

## Scientific references

The physics of the 2.5D FRADO model and its applications are described in:

- **Naddaf, Czerny & Szczerba (2021)**, *ApJ*, **920**, 30 — 2.5D FRADO dynamics and geometry. [DOI: 10.3847/1538-4357/ac139d](https://doi.org/10.3847/1538-4357/ac139d).
- **Naddaf & Czerny (2022)**, *A&A*, **663**, A77 — low-ionization emission-line modelling. [DOI: 10.1051/0004-6361/202142806](https://doi.org/10.1051/0004-6361/202142806).
- **Naddaf et al. (2025)**, *A&A*, **702**, L13 — quasar main-sequence applications. [DOI: 10.1051/0004-6361/202556852](https://doi.org/10.1051/0004-6361/202556852).

Please cite **FRADO Python** using [`CITATION.cff`](CITATION.cff) together with the relevant model paper(s). The software release DOI will be added after Zenodo archiving.

## Documentation

- [Google Colab quickstart](docs/QUICKSTART.md)
- [Physical and numerical reference](docs/REFERENCE.md)
- [Validation and numerical limitations](docs/VALIDATION.md)
- [Reference comparisons](validation/REFERENCE_COMPARISON.md)
- [Convergence and near-grazing study](validation/extended/REPORT.md)

**Author:** Mohammad Hassan Naddaf · [ORCID](https://orcid.org/0000-0002-7604-9594).

**Licensing:** A software license has not yet been assigned; public redistribution awaits confirmation of rights for the bundled opacity and reference data.
