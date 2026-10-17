# FRADO Python compared with supplied MATLAB reference trajectories

Disk contact follows a downward crossing of z=H(rho). The MATLAB reference uses different event semantics.
Final times need not match when the event or sublimation semantics differ.

| Launch [Rg] | MATLAB end [yr] | Python end [yr] | Python status | Python final z-H [Rg] | Max dusty Δpos [Rg] | Max dusty Δvel [km/s] |
|---:|---:|---:|---|---:|---:|---:|
| 8,370 | 5.763116 | 5.819793 | DISK_IMPACT | 4.8316906e-13 | 0.0090 | 0.0110 |
| 22,990 | 33.617466 | 33.592238 | DISK_IMPACT | 0 | 1.0715 | 0.1013 |
| 30,860 | 52.730011 | 52.726855 | DISK_IMPACT | 0 | 0.8362 | 0.0314 |

These are the same three supplied cases, **not** a full re-run of the 2,803-cloud MATLAB grid.
The pre-sublimation comparison is intentionally limited; different state-transition handling
in MATLAB/Python prevents interpreting full-flight discrepancies as purely integration error.
A suite of artificial analytic trajectories independently tests disk impact, near misses,
physical grazing contact and continuation into additional orbits.
MATLAB event IDs 2 and 4 do not have identical semantics in the Python stopping policies.
