"""Optional sweep summary plot, no FRADO integration required.

Usage: python examples/plot_sweep_summary.py OutPut/demo
"""
import argparse
import csv
from pathlib import Path
import matplotlib.pyplot as plt


def plot_summary(folder: Path):
    with (folder/'summary.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError('No clouds found in summary.csv')
    color_for = {'DISK_IMPACT':'tab:blue', 'NEAR_MISS':'tab:orange',
                 'BLR_BOUNDARY':'tab:green', 'MAX_TIME_REACHED':'tab:purple',
                 'NUMERICAL_FAILURE':'tab:red'}
    fig, ax = plt.subplots(figsize=(9,5))
    for reason, color in color_for.items():
        group = [r for r in rows if r['termination_reason']==reason and r['flight_time_years']]
        if group:
            ax.scatter([float(r['launch_radius_rg']) for r in group],
                       [float(r['flight_time_years']) for r in group],
                       label=f'{reason} ({len(group)})',color=color,s=25)
    flags = [r for r in rows if r['numerical_status'] in ('SCREEN_FLAGGED','RECHECK_DISAGREES','RECHECK_FAILED')
             and r['flight_time_years']]
    if flags:
        ax.scatter([float(r['launch_radius_rg']) for r in flags],
                   [float(r['flight_time_years']) for r in flags],
                   facecolors='none',edgecolors='black',s=75,
                   label=f'Numerically flagged ({len(flags)})')
    ax.set(xlabel='Launch radius (Rg)',ylabel='Flight time (years)',
           title='FRADO radius-sweep outcomes')
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    target=folder/'sweep_summary.png'
    fig.savefig(target,dpi=160)
    plt.close(fig)
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir',type=Path)
    args=parser.parse_args()
    print(plot_summary(args.output_dir))
