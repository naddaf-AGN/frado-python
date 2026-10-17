"""Streaming, self-describing FRADO outputs for single clouds or radius sweeps.

One summary record is written per *attempted* launch radius, even on failure.
Events use a long/tidy representation to accommodate multiple orbits.
"""
from dataclasses import asdict
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import csv
import json
import numpy as np
from .constants import YEAR

LAUNCH_SCREENING_FIELDS = ['launch_radius_rg', 'launch_height_rg', 'status',
                           'absorbed_power', 'emission_power', 'heating_to_emission_ratio', 'error']

EVENT_FIELDS = ['launch_radius_rg', 'event_id', 'event_name', 'time_years',
                'x_rg', 'y_rg', 'z_rg', 'rho_rg', 'vx_kms', 'vy_kms', 'vz_kms',
                'dusty_at_event']
SUMMARY_FIELDS = [
    'launch_radius_rg', 'launch_height_rg', 'orbital_period_years',
    'flight_time_years', 'termination_reason', 'termination_mode',
    'end_rho_rg', 'end_z_rg', 'end_vz_kms', 'sublimation_time_years',
    'sublimation_rho_rg', 'sublimation_z_rg', 'number_of_events', 'error',
    # Supplementary physical quantities. Maxima are based on output samples.
    'launch_vy_kms', 'end_x_rg', 'end_y_rg', 'end_vx_kms', 'end_vy_kms',
    'max_z_rg_sampled', 'max_rho_rg_sampled', 'max_spherical_radius_rg_sampled',
    'maximum_speed_kms_sampled', 'first_vertical_turning_time_years',
    'escape_speed_crossed', 'dust_sublimated', 'trajectory_file',
    # Numerical quality labels do not silently change physical classification.
    'numerical_status', 'numerical_flags', 'minimum_gap_fraction_at_turning',
    'normalized_impact_speed', 'recheck_reason', 'recheck_time_difference_orbits',
    'numerical_note',
    # Screening quantities: blank physical trajectory columns mean no launch.
    'launch_status', 'launch_absorbed_power', 'dust_emission_power',
    'launch_heating_ratio'
]
TRAJECTORY_COLUMNS = ['time_years','x_rg','y_rg','z_rg',
                      'vx_kms','vy_kms','vz_kms','sublimated_0_or_1']


def model_label(cfg):
    return (f'T{cfg.sublimation_temperature:4.0f}_A{cfg.patch_size:g}'
            f'_Z{cfg.metallicity:g}_M{cfg.bh_mass_log:.1f}_m{cfg.mdot_log:.1f}')


def save_trajectory(result, cfg, scales, initial_rho, root='OutPut'):
    """Eight-column MATLAB-compatible trajectory, every requested cloud.

    Values are written with higher precision than the original MATLAB text
    so users can perform accurate subsequent analyses without rerunning FRADO.
    """
    label = model_label(cfg)
    folder = Path(root) / 'Trajectory' / label
    folder.mkdir(parents=True, exist_ok=True)
    filename = folder / f'{label}_R{initial_rho:.12g}.dat'
    p = result.states
    values = np.column_stack((result.time_s/YEAR, p[:,0],p[:,2],p[:,4],
        p[:,1]*scales.velocity_kms_factor,p[:,3]*scales.velocity_kms_factor,
        p[:,5]*scales.velocity_kms_factor,(~result.dusty).astype(int)))
    np.savetxt(filename, values, fmt=['%.10e']*7+['%d'],
               header=' '.join(TRAJECTORY_COLUMNS))
    return filename


def event_rows(result, scales):
    for rec in result.event_records:
        s = rec.state
        yield dict(launch_radius_rg=result.initial_radius, event_id=rec.id,
                   event_name=rec.name, time_years=rec.time_s/YEAR,
                   x_rg=s[0], y_rg=s[2], z_rg=s[4],
                   rho_rg=np.hypot(s[0], s[2]),
                   vx_kms=s[1]*scales.velocity_kms_factor,
                   vy_kms=s[3]*scales.velocity_kms_factor,
                   vz_kms=s[5]*scales.velocity_kms_factor,
                   dusty_at_event=int(rec.time_s < result.sublimation_s))


def summary_row(result, cfg, scales, assessment=None, path=None, launch_screen=None):
    s = result.final_state
    ss = result.sublimation_state
    p = result.states
    turn = next((e.time_s/YEAR for e in result.event_records
                 if e.name == 'vertical_turning'), '')
    velocity = p[:,[1,3,5]] * scales.velocity_kms_factor
    row = dict(
        launch_radius_rg=result.initial_radius,
        launch_height_rg=p[0,4],
        orbital_period_years=result.orbital_period_s/YEAR,
        flight_time_years=result.time_s[-1]/YEAR,
        termination_reason=result.termination_reason,
        termination_mode=cfg.termination_mode,
        end_rho_rg=np.hypot(s[0],s[2]), end_z_rg=s[4],
        end_vz_kms=s[5]*scales.velocity_kms_factor,
        sublimation_time_years=result.sublimation_s/YEAR if np.isfinite(result.sublimation_s) else '',
        sublimation_rho_rg=np.hypot(ss[0],ss[2]) if ss is not None else '',
        sublimation_z_rg=ss[4] if ss is not None else '',
        number_of_events=len(result.event_records), error='',
        launch_vy_kms=p[0,3]*scales.velocity_kms_factor,
        end_x_rg=s[0], end_y_rg=s[2],
        end_vx_kms=s[1]*scales.velocity_kms_factor,
        end_vy_kms=s[3]*scales.velocity_kms_factor,
        max_z_rg_sampled=np.max(p[:,4]),
        max_rho_rg_sampled=np.max(np.hypot(p[:,0],p[:,2])),
        max_spherical_radius_rg_sampled=np.max(np.linalg.norm(p[:,[0,2,4]],axis=1)),
        maximum_speed_kms_sampled=np.max(np.linalg.norm(velocity,axis=1)),
        first_vertical_turning_time_years=turn,
        escape_speed_crossed=int(any(e.name == 'escape_speed' for e in result.event_records)),
        dust_sublimated=int(np.isfinite(result.sublimation_s)),
        trajectory_file=str(path) if path is not None else '',
        numerical_status=assessment.status if assessment is not None else 'NOT_CHECKED',
        numerical_flags=';'.join(assessment.flags) if assessment else '',
        minimum_gap_fraction_at_turning=(assessment.minimum_clearance_fraction
                                        if assessment and assessment.minimum_clearance_fraction is not None else ''),
        normalized_impact_speed=(assessment.normalized_impact_speed
                                if assessment and assessment.normalized_impact_speed is not None else ''),
        recheck_reason=(assessment.recheck_reason or '') if assessment else '',
        recheck_time_difference_orbits=(assessment.recheck_time_difference_fraction
                                       if assessment and assessment.recheck_time_difference_fraction is not None else ''),
        numerical_note=assessment.note if assessment else '',
        launch_status='LAUNCHED',
        launch_absorbed_power=(launch_screen.absorbed_power if launch_screen else ''),
        dust_emission_power=(launch_screen.emission_power if launch_screen else ''),
        launch_heating_ratio=(launch_screen.heating_ratio if launch_screen else '')
    )
    return row


class OutputWriter:
    """One row per attempted cloud; never mixes different runs by accident.

    OutputWriter recreates summary/events under the specified output path.
    Use a unique directory for each set of physical parameter values.
    """
    def __init__(self, root, cfg):
        self.directory = Path(root)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.summary_path = self.directory / 'summary.csv'
        self.events_path = self.directory / 'events.csv'
        self.launch_screening_path = self.directory / 'launch_screening.csv'
        for path, fields in ((self.summary_path, SUMMARY_FIELDS),
                             (self.events_path, EVENT_FIELDS),
                             (self.launch_screening_path, LAUNCH_SCREENING_FIELDS)):
            with path.open('w', newline='') as fh:
                csv.DictWriter(fh, fieldnames=fields).writeheader()
        self.started = datetime.now(timezone.utc)
        self.outcomes = Counter()
        self.statuses = Counter()
        self.nattempted = 0
        self.launch_statuses = Counter()
        self.eligible_radius_min = None
        self.eligible_radius_max = None
        metadata = {
            'package_version':'1.0.0',
            'started_utc':self.started.isoformat(),
            'configuration':asdict(cfg),
            'launch_radii_rg':[float(x) for x in cfg.launch_radii],
            'number_of_requested_clouds':len(cfg.launch_radii),
            'trajectory_columns':TRAJECTORY_COLUMNS,
            'summary_columns':SUMMARY_FIELDS,
            'event_columns':EVENT_FIELDS,
            'launch_screening_columns':LAUNCH_SCREENING_FIELDS,
            'units': {'radius':'Rg','height':'Rg','time':'years','velocity':'km/s'},
            'notes': [
                'Physical DISK_IMPACT means a downward crossing or tangential contact with z=H(rho).',
                'NEAR_MISS means positive local minimum of z-H(rho) within the user threshold.',
                'BLR_BOUNDARY and MAX_TIME_REACHED are computational stopping limits.',
                'Sublimation is a one-way event and does not end a trajectory.',
                'Numerical flags are heuristic, not statistical confidence intervals.',
                'RECHECK_AGREES records agreement with one tighter integration, not a proof.',
                'Sampled maximum diagnostics depend on output_samples.',
                'Events may contain multiple vertical turns for multi-orbit runs.',
                'Unlike MATLAB 16-column events, these are one row per detected event.',
                'Before each cloud, dust survival at z=H(rho) is tested using the same absorption law.',
                'NO_DUST_AT_LAUNCH means absorbed power >= dust emission power; no trajectory is computed.',
                'Allowed launch radii are individually screened; no monotonic/contiguous range is assumed.',
                'Dust survival is necessary but not sufficient for a dynamical cloud launch.'
            ]
        }
        (self.directory/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')

    def add_launch_screening(self, screen):
        with self.launch_screening_path.open('a', newline='') as fh:
            csv.DictWriter(fh, fieldnames=LAUNCH_SCREENING_FIELDS).writerow(dict(
                launch_radius_rg=screen.radius_rg, launch_height_rg=screen.height_rg,
                status=screen.status, absorbed_power=screen.absorbed_power,
                emission_power=screen.emission_power,
                heating_to_emission_ratio=screen.heating_ratio, error=''))
        self.launch_statuses[screen.status] += 1
        if screen.can_launch:
            r = screen.radius_rg
            self.eligible_radius_min = (r if self.eligible_radius_min is None
                                        else min(r, self.eligible_radius_min))
            self.eligible_radius_max = (r if self.eligible_radius_max is None
                                        else max(r, self.eligible_radius_max))

    def add_launch_screening_error(self, rho, error):
        with self.launch_screening_path.open('a', newline='') as fh:
            csv.DictWriter(fh, fieldnames=LAUNCH_SCREENING_FIELDS).writerow(dict(
                launch_radius_rg=rho, status='SCREENING_FAILED', error=str(error)))
        self.launch_statuses['SCREENING_FAILED'] += 1

    def add_no_launch(self, screen, cfg, scales):
        """Record non-launch with no fabricated position/velocity trajectory."""
        r, h = screen.radius_rg, screen.height_rg
        orbital = scales.orbital_period_coefficient * np.hypot(r, h)**1.5 / YEAR
        row = dict.fromkeys(SUMMARY_FIELDS, '')
        row.update(launch_radius_rg=r, launch_height_rg=h,
                   orbital_period_years=orbital,
                   termination_reason='NO_DUST_AT_LAUNCH',
                   termination_mode=cfg.termination_mode,
                   launch_status='NO_DUST_AT_LAUNCH',
                   launch_absorbed_power=screen.absorbed_power,
                   dust_emission_power=screen.emission_power,
                   launch_heating_ratio=screen.heating_ratio,
                   number_of_events=0, numerical_status='NOT_APPLICABLE')
        with self.summary_path.open('a', newline='') as fh:
            csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS).writerow(row)
        self.nattempted += 1
        self.outcomes['NO_DUST_AT_LAUNCH'] += 1

    def add(self, result, cfg, scales, assessment=None, path=None, launch_screen=None):
        with self.summary_path.open('a',newline='') as fh:
            csv.DictWriter(fh,fieldnames=SUMMARY_FIELDS).writerow(
                summary_row(result,cfg,scales,assessment,path,launch_screen))
        with self.events_path.open('a',newline='') as fh:
            writer=csv.DictWriter(fh,fieldnames=EVENT_FIELDS)
            writer.writerows(event_rows(result, scales))
        self.nattempted += 1
        self.outcomes[result.termination_reason] += 1
        self.statuses[assessment.status if assessment else 'NOT_CHECKED'] += 1

    def add_failure(self, initial_rho, cfg, error):
        row=dict.fromkeys(SUMMARY_FIELDS,'')
        row.update(launch_radius_rg=initial_rho,termination_reason='NUMERICAL_FAILURE',
                   termination_mode=cfg.termination_mode,error=str(error),
                   numerical_status='NOT_CHECKED', launch_status='FAILED')
        with self.summary_path.open('a',newline='') as fh:
            csv.DictWriter(fh,fieldnames=SUMMARY_FIELDS).writerow(row)
        self.nattempted += 1
        self.outcomes['NUMERICAL_FAILURE'] += 1

    def finalize(self):
        report = {
            'attempted_clouds':self.nattempted,
            'screening_counts':dict(self.launch_statuses),
            'eligible_launches':self.launch_statuses['LAUNCH_ALLOWED'],
            'skipped_no_dust_at_launch':self.launch_statuses['NO_DUST_AT_LAUNCH'],
            'eligible_requested_radius_min_rg':self.eligible_radius_min,
            'eligible_requested_radius_max_rg':self.eligible_radius_max,
            'eligible_radius_note':('Extrema of individually screened requested radii; '
                                    'not a proven continuous launching interval.'),
            'termination_counts':dict(self.outcomes),
            'numerical_status_counts':dict(self.statuses),
            'finished_utc':datetime.now(timezone.utc).isoformat(),
            'paths':{'summary':'summary.csv','events':'events.csv',
                     'launch_screening':'launch_screening.csv',
                     'metadata':'metadata.json',
                     'trajectories':'Trajectory/'}
        }
        (self.directory/'run_report.json').write_text(json.dumps(report,indent=2)+'\n')
        return report
