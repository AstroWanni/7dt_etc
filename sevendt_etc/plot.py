"""Plots for ExposureCalculator (called as calc.plot_*)."""
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from . import constants as C
from .conditions import resolve_conditions


def filter_colors(filters):
    """Broad bands: fixed colors; medium bands: rainbow over 400-875 nm (same color for a band in every plot)."""
    norm, cmap = mcolors.Normalize(*C.MEDIUM_RANGE), plt.get_cmap('rainbow')
    return {f: C.BROAD_COLORS.get(f) or cmap(norm(C.wavelength(f))) for f in filters}


def _unit(mode):
    return 'mag' if mode == 'point' else 'mag/arcsec$^2$'


def _finish(calc, fig, ax, cond, mode, area, title, save_path, return_fig=False):
    """Show the figure once; return it only if requested (a returned Figure is displayed again by Jupyter)."""
    ax.text(0.02, 0.98, calc.conditions_text(cond, mode, area), transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', alpha=0.8))
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=8, ncol=2, loc='lower right')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    if return_fig:
        return fig
    plt.show()
    return None


def plot_snr_vs_exptime(calc, filter, value, mode='point', area=1.0, exptime_range=(1, 1000), n_points=300,
                        n_stack=100, save_path=None, return_fig=False, **cond_kw):
    """
    SNR vs exposure time.
    line    : single exposure, exptime_range [s]
    scatter : N x frame_time frames stacked (N = 1..n_stack), SNR = sqrt(N) x one-frame SNR
    """
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    T = np.logspace(*np.log10(exptime_range), n_points)
    n = np.arange(1, n_stack + 1)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        S = 10 ** (-0.4 * (value - p['zp_e_s'])) * A
        ax.plot(T, calc.snr(S, p, npix, T, n_frame=1), color=col[f], label=f)
        ax.scatter(n * calc.frame_time, calc.snr(S, p, npix, n * calc.frame_time, n_frame=n), color=col[f], s=10)
    ax.plot([], [], 'k-', label='single exposure')
    ax.scatter([], [], color='k', s=10, label=f'stacked {calc.frame_time:g} s frames')
    ax.set(xscale='log', yscale='log', xlabel='exposure time [s] (stacked: N x %g s)' % calc.frame_time, ylabel='SNR')
    return _finish(calc, fig, ax, cond, mode, area, f'SNR vs exposure, {mode} {value} {_unit(mode)}', save_path, return_fig)


def plot_exptime_vs_snr(calc, filter, value, mode='point', area=1.0, snr_range=(1, 100), n_points=60, save_path=None, return_fig=False, **cond_kw):
    """Exposure needed vs target SNR. line: single exposure; scatter: stacked total N x frame_time."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    snrs = np.logspace(*np.log10(snr_range), n_points)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        S = 10 ** (-0.4 * (value - p['zp_e_s'])) * A
        ax.plot(snrs, [calc.exptime_single(S, p, npix, s) for s in snrs], color=col[f], label=f)
        ax.scatter(snrs, [calc.n_frame_stack(S, p, npix, s) * calc.frame_time for s in snrs], color=col[f], s=10)
    ax.plot([], [], 'k-', label='single exposure')
    ax.scatter([], [], color='k', s=10, label=f'stacked {calc.frame_time:g} s frames')
    ax.set(xscale='log', yscale='log', xlabel='SNR', ylabel='exposure time [s]')
    return _finish(calc, fig, ax, cond, mode, area, f'Exposure vs SNR, {mode} {value} {_unit(mode)}', save_path, return_fig)


def plot_limit_vs_exptime(calc, filter, snr=5.0, mode='point', area=1.0, exptime_range=(1, 1000), n_points=300,
                          save_path=None, return_fig=False, **cond_kw):
    """Limiting magnitude (point) or surface brightness (extended) vs single-exposure time."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    T = np.logspace(*np.log10(exptime_range), n_points)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        ax.plot(T, p['zp_e_s'] - 2.5 * np.log10(calc.limit_flux(p, npix, T, snr, n_frame=1) / A), color=col[f], label=f)
    ylabel = f'{snr:g}$\\sigma$ limiting mag [AB]' if mode == 'point' else f'{snr:g}$\\sigma$ limiting SB [mag/arcsec$^2$]'
    ax.set(xscale='log', xlabel='single exposure time [s]', ylabel=ylabel)
    return _finish(calc, fig, ax, cond, mode, area, f'Depth vs single exposure ({mode})', save_path, return_fig)


def plot_noise_budget(calc, filt, value, mode='point', area=1.0, exptime_range=(1, 1000), n_points=300,
                      save_path=None, return_fig=False, **cond_kw):
    """Variance fractions (source, sky, dark, read, systematic) and SNR vs single-exposure time for one filter."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    p = calc.band_params(calc.filter_list(filt)[0], cond)
    T = np.logspace(*np.log10(exptime_range), n_points)
    S = 10 ** (-0.4 * (value - p['zp_e_s'])) * A
    comp = calc.var_components(S, p, npix, T, n_frame=1)
    tot = sum(comp.values())
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.stackplot(T, *[c / tot for c in comp.values()], labels=list(comp), alpha=0.8)
    ax2 = ax.twinx()
    ax2.plot(T, S * T / np.sqrt(tot), 'k-', lw=1.2)
    ax2.set(yscale='log', ylabel='SNR (black line)')
    ax.set(xscale='log', xlim=exptime_range, ylim=(0, 1), xlabel='single exposure time [s]', ylabel='variance fraction')
    return _finish(calc, fig, ax, cond, mode, area, f'Noise budget {filt} (single exposure), {mode} {value} {_unit(mode)}', save_path, return_fig)
