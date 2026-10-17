"""Create standalone figures from saved extended-convergence outputs; no ODE solve."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parent/'extended'
data=json.loads((OUT/'results.json').read_text())

# Numerical sensitivity of a near-grazing REAL cloud.
runs=[r for r in data['real_runs'] if r['radius_rg']==35750]
if runs:
    fig,ax=plt.subplots(figsize=(9,4.5))
    xs=np.arange(len(runs))
    ax.plot(xs,[r['end_years'] for r in runs],marker='o')
    ax.set_xticks(xs,[r['case'].replace('_','\n') for r in runs])
    ax.set_ylabel('Time to disk impact [years]')
    ax.set_title('Near-grazing FRADO trajectory: Rlaunch = 35,750 Rg')
    ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(OUT/'grazing_impact_sensitivity.png',dpi=160)
    plt.close(fig)

# The sublimation timestamp is obtained through ODE root finding, not by
# selecting the nearest sampled output row.
runs=[r for r in data['real_runs'] if r['radius_rg']==8370]
if runs:
    baseline=next(r for r in runs if r['case']=='baseline')
    fig,ax=plt.subplots(figsize=(9,4.5))
    xs=np.arange(len(runs))
    ax.plot(xs,[(r['sublimation_years']-baseline['sublimation_years'])*365*86400
                 for r in runs],marker='o')
    ax.set_xticks(xs,[r['case'].replace('_','\n') for r in runs])
    ax.set_ylabel('Sublimation-time shift relative to baseline [seconds]')
    ax.set_title('Sublimation-event numerical sensitivity: Rlaunch = 8,370 Rg')
    ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(OUT/'sublimation_time_sensitivity.png',dpi=160)
    plt.close(fig)

# Analytic sensitivity diagram for *near-grazing positive b*, showing
# where a normal solver falsely reports impact instead of near miss.
rows=[r for r in data['synthetic'] if r['b']>0 and r['mode']=='first_excursion']
if rows:
    by_tier={name:sorted((r for r in rows if r['tier']==name),key=lambda r:r['b'])
             for name in ('reference','tight')}
    fig,ax=plt.subplots(figsize=(9,4.5))
    for tier, items in by_tier.items():
        y=[1 if r['actual']=='NEAR_MISS' else 0 for r in items]
        ax.plot([r['b'] for r in items],y,'o-',label=tier)
    ax.set_xscale('log');ax.set_yticks([0,1],['False impact','Correct near miss'])
    ax.set_xlabel('Positive grazing control parameter b')
    ax.set_title('Artificial grazing test: impact/near-miss numerical ambiguity')
    ax.legend();ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(OUT/'analytic_grazing_classification.png',dpi=160)
    plt.close(fig)
print('Figures generated in',OUT)
