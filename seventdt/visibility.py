"""Night plan for one target + exposure time per filter at the best (lowest airmass) time."""
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.time import Time
from matplotlib import patheffects
from matplotlib.ticker import LogLocator
from . import constants as C
from .calculator import ExposureCalculator
from .conditions import get_observer, resolve_conditions
from .plot import filter_colors


def night_times(date, n_points=300):
    """Astronomical evening -> morning twilight of the night starting on local date 'YYYY-MM-DD'."""
    obs = get_observer()
    noon = Time(f'{date}T16:00:00')   # ~local noon (UTC-4)
    eve = obs.twilight_evening_astronomical(noon, which='next')
    mor = obs.twilight_morning_astronomical(eve, which='next')
    return Time(np.linspace(eve.mjd, mor.mjd, n_points), format='mjd'), eve, mor


def altitudes(ra, dec, times):
    """(target altitude, Moon altitude) [deg] at the 7DT site."""
    frame = AltAz(obstime=times, location=get_observer().location)
    return (SkyCoord(ra, dec, unit='deg').transform_to(frame).alt.deg,
            get_body('moon', times).transform_to(frame).alt.deg)


def observe_and_etc(ra, dec, date, value, snr=5.0, mode='point', area=100.0, seeing=2.0, aperture=C.DEFAULT_APERTURE,
                    bands=C.MEDIUM_BANDS, calc=None, min_alt=30.0, warn_threshold=3600.0, figsize=(13, 12),
                    save_path=None, verbose=True):
    """
    ra, dec        : target [deg]
    date           : local night start date 'YYYY-MM-DD'
    value          : AB mag (point) or surface brightness [mag/arcsec^2] (extended)
    snr            : target SNR (default 5)
    mode, area     : 'point' or 'extended'; area [arcsec^2] for extended
    seeing         : ["]
    aperture       : 10 (default) or 5 ["] (ignored if calc is given)
    bands          : filters to plot (default m400-m875, 25 nm steps)
    min_alt        : altitude [deg] shaded as unobservable
    warn_threshold : exposure [s] above which the label gets a warning mark
    Returns (result dict per filter, info dict)
    """
    calc = calc or ExposureCalculator(aperture=aperture)
    bands = calc.filter_list(bands)
    times, eve, mor = night_times(date)
    t_alt, m_alt = altitudes(ra, dec, times)
    i = int(np.argmax(t_alt))
    best, best_alt = times[i], float(t_alt[i])
    res = calc.calculate_exptime(bands, value, snr=snr, mode=mode, area=area, ra=ra, dec=dec, obsdate=best,
                                 seeing=seeing, verbose=False)
    cond = resolve_conditions(ra=ra, dec=dec, obsdate=best, seeing=seeing)
    up = t_alt >= min_alt
    info = {'date': date, 'best_time_utc': best.iso, 'best_alt': best_alt, 'airmass': cond['airmass'],
            'moon_phase': cond['moon_phase'], 'moon_sep': cond['moon_separation'], 'seeing': seeing,
            'twilight_utc': (eve.iso, mor.iso),
            'hours_above_min_alt': float(up.sum() * (mor.mjd - eve.mjd) * 24 / len(times)),
            'window_utc': (times[up][0].iso, times[up][-1].iso) if up.any() else None}
    vunit_txt = 'mag' if mode == 'point' else f'mag/arcsec^2, area {area:g} arcsec^2'
    geom_txt = f'{calc.aperture:g}" aperture' if mode == 'point' else f'area {area:g} arcsec^2'
    if verbose:
        print(f"Target RA, Dec = {ra:.4f}, {dec:.4f} | night of {date}")
        print(f"Astronomical twilight (UTC): {eve.iso[:16]} -> {mor.iso[:16]}")
        print(f"Above {min_alt:g} deg (UTC): " + (f"{info['window_utc'][0][:16]} -> {info['window_utc'][1][:16]} "
              f"({info['hours_above_min_alt']:.1f} h)" if up.any() else 'never'))
        print(f"Best time (lowest airmass): {best.iso[:19]} UTC, alt = {best_alt:.1f} deg, airmass = {cond['airmass']:.2f}")
        print(f"Moon phase = {cond['moon_phase']:.2f}, Moon sep = {cond['moon_separation']:.1f} deg, seeing = {seeing:g}\"")
        print(f"{mode} source, {value} {vunit_txt.split(',')[0]}, {geom_txt}, SNR = {snr:g}, frame = {calc.frame_time:g} s\n")
        print(f"{'filter':6} {'T[s]':>10}  frames")
        for f in bands:
            r = res[f]
            print(f"{f:6} {r['exptime']:10.1f}  {r['frames']}" + ('  !' if r['exptime'] > warn_threshold else ''))
    # ---- figure ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [1, 1.15]})
    tp = times.datetime
    ax1.plot(tp, t_alt, color='steelblue', lw=2, label='Target')
    ax1.plot(tp, m_alt, color='orange', lw=1.5, ls='--', label='Moon')
    ax1.axhspan(0, min_alt, color='gray', alpha=0.12, label=f'alt < {min_alt:g}°')
    ax1.axvline(best.datetime, color='red', ls='--', lw=1.5, label='Lowest airmass')
    ax1.scatter([best.datetime], [best_alt], color='red', s=80, zorder=5)
    ax1.set(ylim=(0, 90), xlim=(tp[0], tp[-1]))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax1.tick_params(labelsize=13)
    ax1.set_xlabel('UTC', fontsize=15)
    ax1.set_ylabel('Altitude (deg)', fontsize=15)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=12, loc='upper right')
    am = ax1.secondary_yaxis('right', functions=(lambda a: 1 / np.sin(np.radians(np.clip(a, 5, 90))),
                                                 lambda x: np.degrees(np.arcsin(1 / np.clip(x, 1, 11.5)))))
    am.set_yticks([1.0, 1.1, 1.3, 1.5, 2.0, 3.0])
    am.set_ylabel('Airmass', fontsize=15)
    am.tick_params(labelsize=13)
    ax1.text(0.01, 0.97, f"Night of {date}   (twilight {eve.iso[11:16]} - {mor.iso[11:16]} UTC)\n"
             f"Best: {best.iso[:16]} UTC   alt {best_alt:.1f}°   airmass {cond['airmass']:.2f}\n"
             f"Moon phase {cond['moon_phase']:.2f}   Moon sep {cond['moon_separation']:.1f}°   seeing {seeing:g}\"",
             transform=ax1.transAxes, va='top', fontsize=12, bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.9))
    ax1.set_title(f'RA = {ra:.3f}°, Dec = {dec:.3f}°', fontsize=16)
    col = filter_colors(bands)
    Ts = []
    for f in bands:
        w, T = C.wavelength(f), res[f]['exptime']
        if not np.isfinite(T):
            continue
        Ts.append(T)
        ax2.errorbar(w, T, xerr=C.BAND_WIDTH.get(f, 25.0) / 2, fmt='o', color=col[f], ms=8, elinewidth=3, capsize=4,
                     mec='black', mew=0.8)
        t_ = ax2.text(w, T * 1.15, ('⚠ ' if T > warn_threshold else '') + f'{T:.0f}s', ha='center', va='bottom',
                      fontsize=11, color=col[f], fontweight='bold')
        t_.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])
    waves = [C.wavelength(f) for f in bands]
    ax2.axhline(calc.frame_time, color='gray', ls=':', lw=1)
    ax2.text(max(waves) + 13, calc.frame_time * 1.05, f'{calc.frame_time:g} s = 1 frame', color='gray', fontsize=10,
             va='bottom', ha='right')
    ax2.set_yscale('log')
    if Ts:
        ax2.set_ylim(min(min(Ts), calc.frame_time) / 1.5, max(max(Ts), calc.frame_time) * 2.5)
    ax2.set_xlim(min(waves) - 15, max(waves) + 15)
    ax2.set_xticks(waves)
    ax2.set_xticklabels(bands, rotation=45, fontsize=11)
    ax2.tick_params(axis='y', labelsize=13)
    ax2.set_xlabel('Filter', fontsize=15)
    ax2.set_ylabel(f'Exposure for SNR = {snr:g} (s)', fontsize=15)
    ax2.yaxis.set_major_locator(LogLocator(base=10))
    ax2.grid(axis='y', which='major', alpha=0.4)
    ax2.grid(axis='y', which='minor', visible=False)
    vunit = 'mag' if mode == 'point' else 'mag/arcsec$^2$'
    geom = f'{calc.aperture:g}" aperture' if mode == 'point' else f'area {area:g} arcsec$^2$'
    ax2.set_title(f'{mode} source  {value} {vunit},  {geom},  SNR = {snr:g}  (frame {calc.frame_time:g} s)', fontsize=15)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()
    return res, info
