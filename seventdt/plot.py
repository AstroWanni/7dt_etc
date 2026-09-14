"""
Plotting utilities for 7DTC results.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import patheffects

from .constants import FILTER_WAVE, BROAD_BAND_COLORS


def _get_filter_colors(filters: list) -> dict:
    colors = dict(BROAD_BAND_COLORS)
    mfilters = [f for f in filters if f.startswith("m")]
    if mfilters:
        waves = np.array([FILTER_WAVE[f][0] for f in mfilters])
        norm = mcolors.Normalize(vmin=waves.min(), vmax=waves.max())
        cmap = plt.get_cmap("jet")
        for f, w in zip(mfilters, waves):
            colors[f] = cmap(norm(w))
    return colors


def plot_exptime_per_filter(
    result: dict,
    magnitude: float,
    snr: float,
    warn_threshold: float = 3600.0,
    save_path: str = None,
):
    """
    Filter vs required exposure time plot (point source).

    Parameters
    ----------
    result        : dict returned by calculate_exptime
    magnitude     : AB magnitude (for info box)
    snr           : target SNR (for info box)
    warn_threshold: exposure time [s] above which ⚠ is shown
    save_path     : file path to save the figure (optional)
    """
    filt_list = [f for f in result if f in FILTER_WAVE]
    colors = _get_filter_colors(filt_list)
    cond = result[list(result.keys())[0]]

    fig, ax = plt.subplots(figsize=(14, 6))

    for f in filt_list:
        w, hw = FILTER_WAVE[f]
        t = result[f]['t_exp']
        c = colors.get(f, 'gray')
        ax.errorbar(w, t, xerr=hw/2, fmt='o', color=c, markersize=7,
                    elinewidth=2.5, capsize=4, capthick=2,
                    markeredgecolor='black', markeredgewidth=0.8)
        label = f"{t:.0f}s"
        if t > warn_threshold:
            label = "⚠ " + label
        txt = ax.text(w, t * 1.12, label, ha='center', va='bottom',
                      fontsize=12, color=c, fontweight='bold')
        txt.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])
        ftxt = ax.text(w, t / 1.12, f, ha='center', va='top',
                       fontsize=10, color=c, fontweight='bold')
        ftxt.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])

    cond_text = (
        f"mag = {magnitude:.2f}  |  SNR = {snr:.1f}\n"
        f"moon_phase = {cond['moon_phase']:.2f}  |  "
        f"moon_sep = {cond['moon_separation']:.1f}°  |  "
        f"hour_since_sunset = {cond['hour_since_sunset']:.1f} h  |  "
        f"seeing = {cond['seeing']:.1f}\""
    )
    ax.text(0.01, 0.97, cond_text, transform=ax.transAxes,
            fontsize=10.5, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9))

    ax.set_yscale('log')
    ax.set_xlabel("Wavelength (nm)", fontsize=18)
    ax.set_ylabel("Required exposure time (s)", fontsize=18)
    ax.set_title("7DT Required Exposure Time per Filter", fontsize=19.5)
    ax.tick_params(axis='both', labelsize=15)
    ax.set_xlim(350, 970)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()
