"""Optional Matplotlib view; completely separate from FRADO force calculations."""
from pathlib import Path
import numpy as np
from .constants import YEAR


def plot_trajectory(result, disk=None, matlab_path=None, ax=None,
                    show_sublimation=True, show_disk=True):
    """Plot z versus rho with blue=dusty, red=dust-free.

    Solid = Python; dashed = supplied MATLAB reference (when provided).
    Returns (fig, ax). No plotting is done unless this function is called.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(8,5))
    else:
        fig = ax.figure

    p, t = result.states, result.time_s
    rho = np.hypot(p[:,0],p[:,2]); z=p[:,4]

    if matlab_path is not None:
        mat=np.loadtxt(Path(matlab_path), ndmin=2)
        rm=np.hypot(mat[:,1],mat[:,2]); zm=mat[:,3]
        for is_dusty,color,name in ((True,'tab:blue','MATLAB dusty'),(False,'tab:red','MATLAB dust-free')):
            sel=(mat[:,7]==0)==is_dusty
            if np.any(sel):
                ax.plot(np.where(sel,rm,np.nan),np.where(sel,zm,np.nan),
                        '--',color=color,label=name,lw=1.6)

    if show_sublimation and np.isfinite(result.sublimation_s) and result.sublimation_state is not None:
        mask=t < result.sublimation_s
        st=result.sublimation_state
        rsub=np.hypot(st[0],st[2]); zsub=st[4]
        # Appending exact sublimation position gives both segments a shared end.
        r0=np.concatenate([rho[mask],[rsub]])
        z0=np.concatenate([z[mask],[zsub]])
        r1=np.concatenate([[rsub],rho[~mask]])
        z1=np.concatenate([[zsub],z[~mask]])
        if len(r0)>1: ax.plot(r0,z0,color='tab:blue',label='Python dusty',lw=2)
        if len(r1)>1: ax.plot(r1,z1,color='tab:red',label='Python dust-free',lw=2)
        ax.scatter([rsub],[zsub],color='black',s=30,zorder=6,label='Sublimation')
    else:
        ax.plot(rho,z,color='tab:blue',label='Python dusty',lw=2)

    if show_disk and disk is not None:
        lo,hi=max(float(np.min(rho))*0.95,float(disk.radius[0])), min(float(np.max(rho))*1.05,float(disk.radius[-1]))
        if hi>lo:
            rgrid=np.linspace(lo,hi,300)
            ax.plot(rgrid,disk.at(rgrid),color='dimgray',ls=':',lw=1.5,label='Disk surface H(rho)')
    ax.scatter([rho[0]],[z[0]],color='green',zorder=7,s=35,label='Launch')
    ax.scatter([rho[-1]],[z[-1]],marker='x',s=55,color='black',zorder=7,label=result.termination_reason)
    ax.set(xlabel='Cylindrical radius (Rg)',ylabel='Height (Rg)',
           title=f'FRADO: launch R={result.initial_radius:g} Rg')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    return fig,ax
