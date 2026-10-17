# FRADO Python — convergence and grazing impact validation

## Scope

- Numerically integrate the FRADO force law without altering it.
- Compare ODE tolerances, maximum solver steps, and quadrature tolerances independently where possible.
- Test actual crossings, exact tangencies, and nonimpacting near misses with an independent analytic model.
- Distinguish a numerical convergence test from any MATLAB-model agreement.

## Main numerical findings

- Controlled grazing experiments: **40/44** expected classifications match; **4** are tolerance-sensitive failures. In particular, `b=10^-9` and `10^-8` are genuine near misses but the normal numerical settings can mislabel them as impacts. Tight settings recover the correct distinction.
- The sublimating R=8,370 Rg trajectory is comparatively stable; ODE and quadrature sensitivity should be read separately in the tables below.
- The near-grazing R=35,750 Rg case is MUCH more sensitive: the standard-to-strict terminal-time difference is measured in months, even when all cases ultimately classify the event as disk impact. This is a warning against relying on one solver-tolerance setting near tangencies.
- The validation comparisons vary numerical settings rather than physical prescriptions.

### Sublimation-event residual checks

See `sublimation_diagnostics.json` for consistency checks at each stored sublimation state: absorbed/emitted power ratio, acceleration jump, and saved dust-status labels. The diagnostic must be run after the main sweep with `python validation/extended_diagnostics.py`. It uses the matching quadrature settings for each run.

### Caution about grazing-impact times

The R=35,750 Rg cloud reaches the disk with extremely small relative surface velocity. Under the standard solver, the return time is about 60.47125 years; under the strictest solver it is about 59.86454 years (a difference of roughly 221 days). The shift decreases when tolerances are tightened, but full convergence is not proven. The original MATLAB event records rounded near-zero vertical speeds; this alone cannot prove the reference orbit is precisely tangent.

## Analytic near-grazing tests

Expected classifications: **40/44** matched.

This test **intentionally includes coefficients extremely close to an exact tangency**, where classification can become tolerance-dependent.

The analytic surface clearance is g(t) = 1 − cos(t) + b t². Negative b crosses the surface, b = 0 touches it tangentially, and positive b has a genuine near miss. All cases launch at g=0 with zero vertical speed.

| b | mode | tier | expected | obtained | abs. event-time error [s] |
|---:|---|---|---|---|---:|
| -1.0e-03 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 1.63e-05 |
| -1.0e-03 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 6.06e-08 |
| -1.0e-03 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 1.63e-05 |
| -1.0e-03 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 6.06e-08 |
| -1.0e-04 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 5.17e-05 |
| -1.0e-04 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 1.88e-07 |
| -1.0e-04 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 5.17e-05 |
| -1.0e-04 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 1.88e-07 |
| -1.0e-06 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 0.00889 |
| -1.0e-06 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 1.86e-06 |
| -1.0e-06 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 0.00889 |
| -1.0e-06 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 1.86e-06 |
| -1.0e-08 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 0.000889 |
| -1.0e-08 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 0.000889 |
| -1.0e-08 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 0.000889 |
| -1.0e-08 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 0.000889 |
| -1.0e-09 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 0.000281 |
| -1.0e-09 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 0.000281 |
| -1.0e-09 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 0.000281 |
| -1.0e-09 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 0.000281 |
| 0.0e+00 | disk_impact | reference | DISK_IMPACT | DISK_IMPACT | 2.83e-07 |
| 0.0e+00 | disk_impact | tight | DISK_IMPACT | DISK_IMPACT | 1.01e-11 |
| 0.0e+00 | first_excursion | reference | DISK_IMPACT | DISK_IMPACT | 2.83e-07 |
| 0.0e+00 | first_excursion | tight | DISK_IMPACT | DISK_IMPACT | 1.01e-11 |
| 1.0e-09 | disk_impact | reference | MAX_TIME_REACHED | DISK_IMPACT | 12.6 |
| 1.0e-09 | disk_impact | tight | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-09 | first_excursion | reference | NEAR_MISS | DISK_IMPACT | 2.83e-07 |
| 1.0e-09 | first_excursion | tight | NEAR_MISS | NEAR_MISS | 1.01e-11 |
| 1.0e-08 | disk_impact | reference | MAX_TIME_REACHED | DISK_IMPACT | 12.6 |
| 1.0e-08 | disk_impact | tight | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-08 | first_excursion | reference | NEAR_MISS | DISK_IMPACT | 2.83e-07 |
| 1.0e-08 | first_excursion | tight | NEAR_MISS | NEAR_MISS | 1.01e-11 |
| 1.0e-06 | disk_impact | reference | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-06 | disk_impact | tight | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-06 | first_excursion | reference | NEAR_MISS | NEAR_MISS | 2.83e-07 |
| 1.0e-06 | first_excursion | tight | NEAR_MISS | NEAR_MISS | 1e-11 |
| 1.0e-04 | disk_impact | reference | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-04 | disk_impact | tight | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-04 | first_excursion | reference | NEAR_MISS | NEAR_MISS | 2.85e-07 |
| 1.0e-04 | first_excursion | tight | NEAR_MISS | NEAR_MISS | 3.41e-12 |
| 1.0e-03 | disk_impact | reference | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-03 | disk_impact | tight | MAX_TIME_REACHED | MAX_TIME_REACHED | 0 |
| 1.0e-03 | first_excursion | reference | NEAR_MISS | NEAR_MISS | 3.03e-07 |
| 1.0e-03 | first_excursion | tight | NEAR_MISS | NEAR_MISS | 5.78e-11 |

## Full-physics FRADO solver sweeps

Unless otherwise stated, M = 10^7 Msun, mdot = 10^-0.4, Z = 5, Tsub=1500 K, patch size=3, disk_impact mode, 8,001 output samples, max orbits=2.

| Rlaunch [Rg] | solver | last time [yr] | sublimation [yr] | impact clearance [Rg] | status | run [s] |
|---:|---|---:|---:|---:|---|---:|
| 8,370 | baseline | 5.819792927 | 1.732858752 | 4.83e-13 | DISK_IMPACT | 13.3 |
| 8,370 | ode_only | 5.819775025 | 1.732858936 | 4.83e-13 | DISK_IMPACT | 14.9 |
| 8,370 | step_only | 5.819791216 | 1.732858752 | 3.55e-13 | DISK_IMPACT | 12.2 |
| 8,370 | tight_ode | 5.819775023 | 1.732858936 | 2.13e-13 | DISK_IMPACT | 16.9 |
| 8,370 | very_tight_ode | 5.819771719 | 1.732858844 | -1.56e-13 | DISK_IMPACT | 19.5 |
| 8,370 | tighter_quadrature | 5.819792846 | 1.732858706 | 3.41e-13 | DISK_IMPACT | 27.5 |
| 22,990 | baseline | 33.592238307 | — | 0 | DISK_IMPACT | 18.9 |
| 22,990 | tight_ode | 33.591557242 | — | 2.84e-14 | DISK_IMPACT | 23.9 |
| 35,750 | baseline | 60.471246711 | — | 0 | DISK_IMPACT | 14.6 |
| 35,750 | tight_ode | 59.938871367 | — | 0 | DISK_IMPACT | 18.3 |
| 35,750 | very_tight_ode | 59.868667197 | — | 0 | DISK_IMPACT | 25.5 |
| 35,750 | ultra_tight_ode | 59.864537209 | — | 0 | DISK_IMPACT | 41.1 |

### Cross-configuration differences (relative to baseline)

| Launch [Rg] | experiment | Δend [days] | Δsublimation [days] | Max Δposition [Rg] | Max Δspeed [km/s] |
|---:|---|---:|---:|---:|---:|
| 8,370 | ode_only | -0.00653423 | +6.7371e-05 | 0.125188 | 0.037827 |
| 8,370 | step_only | -0.000624508 | +0 | 0.0197908 | 0.00774364 |
| 8,370 | tight_ode | -0.006535 | +6.7371e-05 | 0.125312 | 0.0379073 |
| 8,370 | very_tight_ode | -0.00774112 | +3.3735e-05 | 0.147699 | 0.0446876 |
| 8,370 | tighter_quadrature | -2.98853e-05 | -1.67547e-05 | 0.000256185 | 6.57739e-05 |
| 22,990 | tight_ode | -0.248589 | — | 0.226194 | 0.017528 |
| 35,750 | tight_ode | -194.317 | — | 1.42929 | 0.0623621 |
| 35,750 | very_tight_ode | -219.942 | — | 1.6266 | 0.0708419 |
| 35,750 | ultra_tight_ode | -221.449 | — | 1.63663 | 0.0712766 |

## Important caveats

1. Analytic near-grazing cases validate event-detection logic, not the full FRADO force law.
2. MAT-file near-zero rounded final vertical velocity is only a *candidate* for grazing.
3. `tight_ode` and `very_tight_ode` tighten both ODE tolerances and maximum step. `ode_only` and `step_only` separate the two effects.
4. Event times are checked at the solver root; plotted output sampling does not set the accuracy.
5. Position/speed differences use PCHIP interpolation of output samples at shared physical times; they include a small output-interpolation error, not purely solver error.
6. A near tangency can change outcome when parameters change, and exactly tangent contacts are intrinsically sensitive to numerical error.
7. The Python runs do not imply full-grid FRADO validation or physical model correctness.

Reproduce: `python validation/extended_convergence.py --suite full` (recompute every run),
or add `--resume` to reuse saved runs whose settings match.
Synthetic-only: `python validation/extended_convergence.py --suite synthetic`.


## Figures generated from the saved outputs

- `grazing_impact_sensitivity.png`: termination-time sensitivity near grazing.
- `sublimation_time_sensitivity.png`: differences in sublimation times.
- `analytic_grazing_classification.png`: false-impact classifications at ordinary tolerances in an exactly solvable toy model.

Run `python validation/plot_extended.py` without re-solving the ODEs.
