"""FRADO physical model and streaming launch-radius sweep controller."""
from dataclasses import dataclass, replace
from .config import Config
from .constants import scales
from .disk_geometry import build_disk
from .dust import emission_power
from .shielding import PatchShielding
from .radiation import Radiation
from .dynamics import Dynamics
from .integrator import integrate_trajectory
from .reliability import screen_trajectory, assess_trajectory
from .outputs import save_trajectory, OutputWriter
from .launching import assess_launch, NoDustAtLaunchError


@dataclass
class Simulation:
    cfg: Config
    shielding_factory: object = PatchShielding

    def __post_init__(self):
        self.scales = scales(self.cfg)
        self.disk = build_disk(self.cfg)
        self.radiation = Radiation(self.cfg, self.scales, self.shielding_factory(self.cfg))
        self.radiation.emission_power = emission_power(self.cfg, self.scales)
        self.dynamics = Dynamics(self.scales, self.radiation)

    def screen_launch(self, rho):
        """Check dusty launching at one radius without integrating motion."""
        return assess_launch(rho, self.disk, self.radiation)

    def trajectory(self, rho, launch_screen=None):
        """Integrate one dusty trajectory; raise if dust cannot survive at launch."""
        return integrate_trajectory(
            self.cfg, self.scales, self.disk, self.radiation, self.dynamics, rho,
            launch_absorbed_power=(launch_screen.absorbed_power
                                   if launch_screen is not None else None))

    def _assess(self, result, rho):
        if self.cfg.numerical_check != 'recheck':
            return assess_trajectory(result, self.disk, self.cfg)
        basic = screen_trajectory(result, self.disk, self.cfg)
        if basic.status != 'SCREEN_FLAGGED':
            return basic
        refined_cfg = replace(
            self.cfg,
            relative_tolerance=self.cfg.relative_tolerance * self.cfg.recheck_tolerance_factor,
            absolute_tolerance=self.cfg.absolute_tolerance * self.cfg.recheck_tolerance_factor,
            max_step_orbit_fraction=self.cfg.max_step_orbit_fraction * self.cfg.recheck_step_factor)
        try:
            # Keep quadrature/forces fixed to isolate ODE sensitivity.
            refined = integrate_trajectory(refined_cfg, self.scales, self.disk,
                                           self.radiation, self.dynamics, rho)
            return assess_trajectory(result, self.disk, self.cfg,
                                     refined_trajectory=refined)
        except (ArithmeticError, RuntimeError, ValueError, OverflowError, FloatingPointError) as err:
            return assess_trajectory(result, self.disk, self.cfg,
                                     recheck_error=err)

    def run(self, output_dir='OutPut', on_error='raise', save_trajectories=True,
            collect_results=True, progress=None):
        """Screen every requested radius, then integrate only dusty-launchable ones.

        Each disallowed launch is represented in summary.csv and
        launch_screening.csv, never by a fake dust-free trajectory.

        Args:
          on_error: 'continue' records a NUMERICAL_FAILURE and proceeds.
          collect_results: False avoids accumulating trajectories in RAM.
          progress: optional callable(rho, reason, index, count), called once/radius.

        Returns: list of (launch_rho, trajectory, path) when
          collect_results=True; otherwise an empty list. All results are saved
          independently via OutputWriter in both cases.
        """
        if on_error not in ('raise', 'continue'):
            raise ValueError('on_error must be raise or continue')
        writer = OutputWriter(output_dir, self.cfg)
        results = []
        try:
            for i, rho in enumerate(self.cfg.launch_radii, start=1):
                rho = float(rho)
                reason = 'NUMERICAL_FAILURE'
                launch_screen = None
                try:
                    launch_screen = self.screen_launch(rho)
                    writer.add_launch_screening(launch_screen)
                    if not launch_screen.can_launch:
                        writer.add_no_launch(launch_screen, self.cfg, self.scales)
                        reason = 'NO_DUST_AT_LAUNCH'
                        continue
                    result = self.trajectory(rho, launch_screen=launch_screen)
                    assessment = self._assess(result, rho)
                    path = (save_trajectory(result, self.cfg, self.scales, rho, output_dir)
                            if save_trajectories else None)
                    writer.add(result, self.cfg, self.scales, assessment, path,
                               launch_screen=launch_screen)
                    reason = result.termination_reason
                    if collect_results:
                        results.append((rho, result, path))
                except NoDustAtLaunchError as err:
                    # Defensive fallback if a caller supplies inconsistent launch
                    # data; never write a fabricated cloud trajectory.
                    if launch_screen is not None:
                        from .launching import LaunchAssessment
                        actual = LaunchAssessment(
                            rho, launch_screen.height_rg, err.absorbed_power,
                            err.emission_power, err.absorbed_power/err.emission_power, False)
                        writer.add_no_launch(actual, self.cfg, self.scales)
                    else:
                        writer.add_launch_screening_error(rho, err)
                        writer.add_failure(rho, self.cfg, err)
                    reason = 'NO_DUST_AT_LAUNCH'
                except (ArithmeticError, RuntimeError, ValueError, OverflowError, FloatingPointError) as err:
                    if launch_screen is None:
                        writer.add_launch_screening_error(rho, err)
                    writer.add_failure(rho, self.cfg, err)
                    if on_error == 'raise':
                        raise
                finally:
                    if progress is not None:
                        progress(rho, reason, i, len(self.cfg.launch_radii))
        finally:
            writer.finalize()
        return results
