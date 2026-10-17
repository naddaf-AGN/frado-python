"""CLI: FRADO single trajectory or complete, streamed radius sweep.

For reproducible runs, optionally provide a JSON configuration file:
    python run.py --config examples/reference_sweep.json --dry-run
    python run.py --config examples/reference_sweep.json --output-dir OutPut/my_run
"""
import argparse
import json
from dataclasses import fields, replace
from pathlib import Path
from frado import Config, Simulation


def build_parser():
    p = argparse.ArgumentParser(description='FRADO Python 1.0.0 — single cloud or launch-radius sweep')
    p.add_argument('--config', type=Path, help='JSON object of Config parameter values')
    radii = p.add_mutually_exclusive_group()
    radii.add_argument('--radius', type=float, help='single launch radius [Rg]')
    radii.add_argument('--radius-start', type=float, help='first launch radius [Rg]')
    p.add_argument('--radius-end', type=float, help='last launch radius [Rg]; requires --radius-start')
    p.add_argument('--radius-step', type=float, help='spacing between launch radii [Rg]')
    p.add_argument('--bh-mass-log', type=float)
    p.add_argument('--mdot-log', type=float)
    p.add_argument('--metallicity', type=float)
    p.add_argument('--patch-size', type=float)
    p.add_argument('--termination', choices=['disk_impact','first_excursion'])
    p.add_argument('--near-miss-fraction', type=float)
    p.add_argument('--max-orbits', type=float)
    p.add_argument('--samples', type=int)
    p.add_argument('--radiation-backend', choices=['scipy','numba'],
                   help='SciPy reference or compiled Numba/QUADPACK quadrature')
    p.add_argument('--numerical-check', choices=['off', 'screen', 'recheck'])
    p.add_argument('--output-dir', type=Path, default=Path('OutPut/frado_run'))
    p.add_argument('--on-error', choices=['raise','continue'], default='continue')
    p.add_argument('--skip-trajectories', action='store_true', help='write summary/events only')
    p.add_argument('--dry-run', action='store_true', help='print configuration and cloud count; do not integrate')
    return p


def parse_config(args):
    params = {}
    if args.config is not None:
        data = json.loads(args.config.read_text())
        if not isinstance(data, dict):
            raise ValueError('Configuration JSON must contain an object')
        fields_known = {f.name for f in fields(Config)}
        unknown = set(data)-fields_known
        if unknown:
            raise ValueError(f'Unknown Config entries: {sorted(unknown)}')
        params.update(data)
    overrides = {
        'bh_mass_log':args.bh_mass_log, 'mdot_log':args.mdot_log,
        'metallicity':args.metallicity, 'patch_size':args.patch_size,
        'termination_mode':args.termination, 'near_miss_fraction':args.near_miss_fraction,
        'maximum_orbits':args.max_orbits, 'output_samples':args.samples,
        'numerical_check':args.numerical_check,
        'radiation_backend':args.radiation_backend
    }
    params.update({k:v for k,v in overrides.items() if v is not None})
    if args.radius is not None:
        if args.radius_end is not None or args.radius_step is not None:
            raise ValueError('--radius cannot be combined with --radius-end / --radius-step')
        params.update(initial_radius=args.radius, final_radius=args.radius)
    elif args.radius_start is not None:
        params['initial_radius'] = args.radius_start
        params['final_radius'] = args.radius_start if args.radius_end is None else args.radius_end
    elif args.radius_end is not None:
        raise ValueError('--radius-end requires --radius-start')
    if args.radius_step is not None:
        params['radius_step'] = args.radius_step
    # A safe, explicit single-cloud CLI default if nothing else is provided.
    if args.config is None and args.radius is None and args.radius_start is None:
        params.update(initial_radius=8370., final_radius=8370.)
    return Config(**params)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg = parse_config(args)
        count = len(cfg.launch_radii)
    except (ValueError, TypeError, OSError) as err:
        parser.error(str(err))
    print(f'FRADO: {count} cloud(s), {cfg.initial_radius:g}–{cfg.final_radius:g} Rg, '
          f'step {cfg.radius_step:g} Rg; termination={cfg.termination_mode}; '
          f'numerical_check={cfg.numerical_check}; radiation={cfg.radiation_backend}', flush=True)
    if args.dry_run:
        print(json.dumps({'clouds':count,'launch_radii_rg':cfg.launch_radii.tolist(),
                          'output_dir':str(args.output_dir)}, indent=2))
        return
    sim = Simulation(cfg)
    def progress(rho, reason, index, total):
        print(f'[{index}/{total}] R={rho:g} Rg → {reason}', flush=True)
    sim.run(args.output_dir, on_error=args.on_error,
            save_trajectories=not args.skip_trajectories,
            collect_results=False, progress=progress)
    print('Saved:', args.output_dir / 'summary.csv',
          args.output_dir / 'events.csv',
          args.output_dir / 'launch_screening.csv',
          args.output_dir / 'metadata.json',
          args.output_dir / 'run_report.json', sep='\n  ')


if __name__ == '__main__':
    main()
