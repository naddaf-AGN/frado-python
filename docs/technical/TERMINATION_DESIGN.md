# Trajectory termination and event semantics

## Disk geometry

The disk surface is `H(rho)` from `disk_geometry.py`, where `rho=sqrt(x*x+y*y)`. Define the clearance `g=z-H(rho)`. The code computes `H` using the *current* cloud radius at each event evaluation, avoiding the stale `global rho` in the original MATLAB event function.

- **Impact crossing:** a **downward** zero of `g` (SciPy event `direction=-1`). A tangential minimum of `g` sufficiently close to zero (`g <= 1e-8*H`, with one-sided numerical tolerance) is also counted as a touch. This tolerance is for contact detection only; it is not a physical disk thickness.
- **Cloud departure:** the original initial condition is `z=H`, `vz=0`; an impact event armed at that exact instant would incorrectly terminate the simulation. The code first requires the cloud to reach `g >= 1e-6 * max(1,H_launch)` to arm its *return* event. If it moves down into the disk before departing, the initial motion is terminated just below the surface at `g=-1e-6 * max(1,H_launch)` to avoid a false solver root at the initial equality. Thus the arm offset is strictly numerical, not a physical continuation model.
- **Clearance turning points:** the derivative is computed as `dg/dt = vz - H'(rho) (x*vx+y*vy)/rho`. The disk's derivative comes from the same PCHIP interpolation. A clearance maximum (`dg/dt` crossing + to −) arms the next clearance minimum (`dg/dt` crossing − to +). The fact that `vz=0` does **not** generally imply minimum separation from a radially varying disk surface.
- **Near miss:** with `termination_mode='first_excursion'`, a local clearance minimum of strictly positive `g` terminates only if `g <= near_miss_fraction*H(rho)`. The default fraction 0.05 is a **configurable modeling choice**, not a derived physical thickness. The result is `NEAR_MISS`, not `DISK_IMPACT`. Otherwise the calculation resumes and another orbit/excursion is possible.
- **Continue after miss:** with `termination_mode='disk_impact'`, even a very small positive local separation is *not* a collision and the integrator continues, subject to physical/numerical stop limits.

## Other stopping conditions

- `BLR_BOUNDARY`: outward crossing of spherical radius equal to `outer_blr_radius`.
- `MAX_TIME_REACHED`: integration exhausted `maximum_orbits * initial_orbital_period`. The default is **5 orbits**, to avoid an unexpectedly expensive calculation when a cloud never intersects the disk.
- `NUMERICAL_FAILURE`: an exception in the solver, force quadrature or disk evaluation; with `on_error='continue'` a failure record is saved and the next cloud is attempted.
- `escape_speed`: upward crossing of the local escape-speed threshold is **diagnostic only** and does not terminate the trajectory (as in the MATLAB model).
- `vertical_turning`: downward zero of `vz`, diagnostic only.
- `sublimation`: radiative absorbed power crosses the thermal emission threshold from below, causing a clean dusty-to-dust-free solver restart. It never marks trajectory termination by itself.

## Numerical limitations / future validation

Numerical thresholds for contact and departure should be studied under tolerance and step-size convergence, especially for grazing orbits. The `max_step_orbit_fraction=0.02` step cap is a practical safeguard against bypassing short-duration surface crossings and may affect runtime; it is not a guarantee of robust event discovery for every conceivable orbit. New choices must be benchmarked on representative models before scientific adoption.

The documented `first_excursion` policy does **not** universally terminate after one orbit: if the cloud stays far from the disk at the first clearance minimum, it continues. This is deliberate and differs from stopping indiscriminately at the first vertical velocity reversal.
