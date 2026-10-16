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


def _finish(calc, fig, ax, cond, mode, area, title, save_path):
    ax.text(0.02, 0.98, calc.conditions_text(cond, mode, area), transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', alpha=0.8))
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=8, ncol=2, loc='lower right')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()
    return fig


def plot_snr_vs_exptime(calc, filter, value, mode='point', area=1.0, exptime_range=(1, 3600), n_points=300,
                        show_sqrt_t=True, save_path=None, **cond_kw):
    """SNR vs total exposure. Dotted (point only): UL5 sqrt(t) scaling without read/source noise."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    T = np.logspace(*np.log10(exptime_range), n_points)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        ax.plot(T, calc.snr(10 ** (-0.4 * (value - p['zp_e_s'])) * A, p, npix, T), color=col[f], label=f)
        if show_sqrt_t and mode == 'point':
            ax.plot(T, calc.snr_sqrt_t(value, p, T), color=col[f], ls=':', lw=0.8)
    ax.set(xscale='log', yscale='log', xlabel='total exposure time [s]', ylabel='SNR')
    note = ' (dotted: sqrt(t) scaling)' if show_sqrt_t and mode == 'point' else ''
    return _finish(calc, fig, ax, cond, mode, area, f'SNR vs exposure, {mode} {value} {_unit(mode)}{note}', save_path)


def plot_exptime_vs_snr(calc, filter, value, mode='point', area=1.0, snr_range=(1, 100), n_points=60, save_path=None, **cond_kw):
    """Total exposure needed vs target SNR."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    snrs = np.logspace(*np.log10(snr_range), n_points)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        S = 10 ** (-0.4 * (value - p['zp_e_s'])) * A
        ax.plot(snrs, [calc.exptime(S, p, npix, s)[0] for s in snrs], color=col[f], label=f)
    ax.set(xscale='log', yscale='log', xlabel='SNR', ylabel='total exposure time [s]')
    return _finish(calc, fig, ax, cond, mode, area, f'Exposure vs SNR, {mode} {value} {_unit(mode)}', save_path)


def plot_limit_vs_exptime(calc, filter, snr=5.0, mode='point', area=1.0, exptime_range=(1, 3600), n_points=300,
                          save_path=None, **cond_kw):
    """Limiting magnitude (point) or surface brightness (extended) vs total exposure."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    T = np.logspace(*np.log10(exptime_range), n_points)
    filters = calc.filter_list(filter)
    col = filter_colors(filters)
    fig, ax = plt.subplots(figsize=(8, 6))
    for f in filters:
        p = calc.band_params(f, cond)
        ax.plot(T, p['zp_e_s'] - 2.5 * np.log10(calc.limit_flux(p, npix, T, snr) / A), color=col[f], label=f)
    ylabel = f'{snr:g}$\\sigma$ limiting mag [AB]' if mode == 'point' else f'{snr:g}$\\sigma$ limiting SB [mag/arcsec$^2$]'
    ax.set(xscale='log', xlabel='total exposure time [s]', ylabel=ylabel)
    return _finish(calc, fig, ax, cond, mode, area, f'Depth vs exposure ({mode})', save_path)


def plot_noise_budget(calc, filt, value, mode='point', area=1.0, exptime_range=(1, 3600), n_points=300,
                      save_path=None, **cond_kw):
    """Variance fractions (source, sky, dark, read, systematic) and SNR vs total exposure for one filter."""
    cond = resolve_conditions(**cond_kw)
    npix, A = calc.geometry(mode, area)
    p = calc.band_params(calc.filter_list(filt)[0], cond)
    T = np.logspace(*np.log10(exptime_range), n_points)
    S = 10 ** (-0.4 * (value - p['zp_e_s'])) * A
    comp = calc.var_components(S, p, npix, T)
    tot = sum(comp.values())
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.stackplot(T, *[c / tot for c in comp.values()], labels=list(comp), alpha=0.8)
    ax2 = ax.twinx()
    ax2.plot(T, S * T / np.sqrt(tot), 'k-', lw=1.2)
    ax2.set(yscale='log', ylabel='SNR (black line)')
    ax.set(xscale='log', xlim=exptime_range, ylim=(0, 1), xlabel='total exposure time [s]', ylabel='variance fraction')
    return _finish(calc, fig, ax, cond, mode, area, f'Noise budget {filt}, {mode} {value} {_unit(mode)}', save_path)
