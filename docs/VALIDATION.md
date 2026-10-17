# Scientific verification and limitations

## Verification provided

- **Automated regression tests:** 43 tests passed in the prepared environment; rerun `python -m pytest -q` on the checked-out code. These cover reference initial conditions, dust-survival screening, integration events, grazing tests, file output, and consistency of the radiation backends.
- **Reference output checks:** for launch radii 8370, 22990 and 30860 `Rg`, the initial disk heights and Keplerian velocities are compared with supplied MATLAB outputs. Reference files are in `validation/reference/`; selected comparisons are documented in [`validation/REFERENCE_COMPARISON.md`](../validation/REFERENCE_COMPARISON.md). These are **not** complete end-to-end agreement tests with MATLAB.
- **Backend comparison:** at six sample locations, SciPy and Numba radiation integrals agreed within approximately 6.4×10⁻¹² relative difference. Selected trajectories at 8370 and 22990 `Rg` agreed within approximately 10⁻⁶ `Rg` in position and shared termination outcomes under tested settings.
- **Solver convergence:** synthetic analytic grazing examples and selected full-physics trajectories are analyzed in [`validation/extended/REPORT.md`](../validation/extended/REPORT.md). Close tangencies can remain numerically ambiguous.

## Performance observations

In one measured environment, compiled Numba callbacks accelerated selected steady-state radiation evaluations by about 36× and two tested trajectories by about 30–41×. These benchmarks are not universal performance guarantees. Reproduction scripts and records appear in `validation/` and [radiation backend notes](technical/RADIATION_BACKENDS.md).

## Interpretation constraints

1. No full grid-scale convergence proof or exhaustive cross-implementation validation is provided for all AGN parameters.
2. Grazing disk impacts and near misses may change termination times or classifications when solver tolerances and maximum step are varied. The `screen` and `recheck` options flag possible uncertainty; they are **not confidence intervals**.
3. The dust-survival condition is a necessary, **not sufficient**, criterion for radiative launch.
4. The fixed opacity fits, disk geometry and shielding are specific model assumptions, not observationally established universal laws.
5. A maximum integration time is a numerical stopping criterion; it is not evidence that a cloud physically escaped or returned.
6. Changing solver backends or parameters should be accompanied by convergence checks for the scientific cases of interest.

These limitations must be retained in publications or applications based on the code.
