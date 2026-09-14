"""
Visibility and scheduling utilities for 7DT.

observe_and_etc()
    Given a target (ra, dec) and a date, computes the altitude curve for the night,
    finds the optimal observation time, and runs the ETC at that time.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib import patheffects
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, get_body

from .core import get_observer
from .etc import calculate_exptime
from .constants import FILTER_WAVE, BROAD_BAND_COLORS, APERTURE_AREA, APERTURE_DIAMETER


# ----------------------------------------------------------------
# Night / altitude helpers
# ----------------------------------------------------------------

def get_night_times(observer, date: str, n_points: int = 300):
    """
    Return time array spanning astronomical twilight for the given date.

    Parameters
    ----------
    observer : astroplan.Observer
    date     : 'YYYY-MM-DD' (local date at the site)
    n_points : number of time steps

    Returns
    -------
    times, eve_twilight, mor_twilight : Time arrays
    """
    # UTC 16:00 ≈ local noon at the 7DT site (UTC-4)
    noon = Time(f"{date}T16:00:00")
    eve_twilight = observer.twilight_evening_astronomical(noon, which='next')
    mor_twilight = observer.twilight_morning_astronomical(eve_twilight, which='next')
    times = Time(np.linspace(eve_twilight.mjd, mor_twilight.mjd, n_points), format='mjd')
    return times, eve_twilight, mor_twilight


def get_altitudes(observer, ra: float, dec: float, times: Time) -> np.ndarray:
    """Target altitude [deg] at each time step."""
    target = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    altaz = target.transform_to(AltAz(obstime=times, location=observer.location))
    return altaz.alt.deg


def get_moon_altitudes(observer, times: Time) -> np.ndarray:
    """Moon altitude [deg] at each time step."""
    moon = get_body('moon', times)
    altaz = moon.transform_to(AltAz(obstime=times, location=observer.location))
    return altaz.alt.deg


def find_optimal_time(times: Time, altitudes: np.ndarray, min_alt: float = 30.0):
    """
    Return the time and altitude of the best observation window.

    Prefers the highest-altitude moment within the min_alt constraint.
    Falls back to the global maximum if the target never reaches min_alt.
    """
    mask = altitudes >= min_alt
    idx = int(np.argmax(altitudes * mask)) if mask.any() else int(np.argmax(altitudes))
    return times[idx], float(altitudes[idx])


# ----------------------------------------------------------------
# Plot helper
# ----------------------------------------------------------------

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


# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

def observe_and_etc(
    ra: float,
    dec: float,
    date: str,
    magnitude: float,
    snr: float,
    filters: str | list = 'all',
    seeing: float = None,
    min_alt: float = 30.0,
    warn_threshold: float = 3600.0,
    save_path: str = None,
) -> dict:
    """
    Altitude-time plot for one night + ETC at the optimal observation time.

    Parameters
    ----------
    ra, dec       : target coordinates [deg]
    date          : observation date 'YYYY-MM-DD' (local date at the site)
    magnitude     : source AB magnitude
    snr           : target SNR
    filters       : filter name, list, or 'all'
    seeing        : seeing FWHM [arcsec] (default: 2.0")
    min_alt       : minimum altitude [deg] for optimal-time search
    warn_threshold: exposure time [s] above which ⚠ is shown on the plot
    save_path     : file path to save the figure (optional)

    Returns
    -------
    dict  — ETC result at the optimal time (same structure as calculate_exptime)
    """
    observer = get_observer()

    # Night times and altitudes
    times, eve_twilight, mor_twilight = get_night_times(observer, date)
    target_alts = get_altitudes(observer, ra, dec, times)
    moon_alts   = get_moon_altitudes(observer, times)

    # Optimal time
    opt_time, opt_alt = find_optimal_time(times, target_alts, min_alt)
    sunset = observer.sun_set_time(opt_time, which='previous', horizon=-18 * u.deg)
    hour_since_sunset = (opt_time.mjd - sunset.mjd) * 24

    # ETC at optimal time
    result = calculate_exptime(
        magnitude=magnitude, snr=snr, filters=filters,
        ra=ra, dec=dec, obsdate=opt_time, seeing=seeing, verbose=False,
    )

    cond = result[list(result.keys())[0]]
    filt_list = [f for f in result if f in FILTER_WAVE]
    colors = _get_filter_colors(filt_list)

    # ---- Figure ----
    fig, (ax_alt, ax_etc) = plt.subplots(2, 1, figsize=(13, 12))

    # Upper: altitude vs time
    ax_alt.plot(times.datetime, target_alts, color='steelblue', lw=2, label='Target')
    ax_alt.plot(times.datetime, moon_alts, color='orange', lw=1.5, ls='--', label='Moon')
    ax_alt.axhline(min_alt, color='gray', ls=':', lw=1, alpha=0.7,
                   label=f'Min altitude ({min_alt}°)')
    ax_alt.axvline(opt_time.datetime, color='red', ls='--', lw=1.5, label='Optimal time')
    ax_alt.scatter([opt_time.datetime], [opt_alt], color='red', zorder=5, s=80)
    ax_alt.axvspan(eve_twilight.datetime, mor_twilight.datetime, alpha=0.05, color='navy')

    info = (
        f"Optimal time : {opt_time.iso} UTC\n"
        f"Target altitude : {opt_alt:.1f}°\n"
        f"Moon phase : {cond['moon_phase']:.2f}  |  "
        f"Moon sep : {cond['moon_separation']:.1f}°  |  "
        f"Hour since sunset : {hour_since_sunset:.1f} h  |  "
        f"Seeing : {cond['seeing']:.1f}\""
    )
    ax_alt.text(0.01, 0.97, info, transform=ax_alt.transAxes,
                fontsize=13.5, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9))

    ax_alt.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax_alt.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    plt.setp(ax_alt.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=15)
    ax_alt.tick_params(axis='y', labelsize=15)
    ax_alt.set_xlabel('UTC', fontsize=18)
    ax_alt.set_ylabel('Altitude (deg)', fontsize=18)
    ax_alt.set_title(f'Altitude vs Time  |  RA={ra:.3f}°, Dec={dec:.3f}°  |  {date}',
                     fontsize=19.5)
    ax_alt.set_ylim(0, 90)
    ax_alt.grid(alpha=0.3)
    ax_alt.legend(fontsize=15, loc='upper right')

    # Lower: ETC
    for f in filt_list:
        w, hw = FILTER_WAVE[f]
        t = result[f]['t_exp']
        c = colors.get(f, 'gray')
        ax_etc.errorbar(w, t, xerr=hw/2, fmt='o', color=c, markersize=7,
                        elinewidth=2.5, capsize=4, capthick=2,
                        markeredgecolor='black', markeredgewidth=0.8)
        label = f"{t:.0f}s"
        if t > warn_threshold:
            label = "⚠ " + label
        txt = ax_etc.text(w, t * 1.12, label, ha='center', va='bottom',
                          fontsize=12, color=c, fontweight='bold')
        txt.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])
        ftxt = ax_etc.text(w, t / 1.12, f, ha='center', va='top',
                           fontsize=10, color=c, fontweight='bold')
        ftxt.set_path_effects([patheffects.withStroke(linewidth=2, foreground='black')])

    ax_etc.set_yscale('log')
    ax_etc.set_xlabel('Wavelength (nm)', fontsize=18)
    ax_etc.set_ylabel('Required exposure time (s)', fontsize=18)
    ax_etc.set_title(f'ETC at Optimal Time  |  mag={magnitude:.2f}, SNR={snr:.1f}',
                     fontsize=19.5)
    ax_etc.tick_params(axis='both', labelsize=15)
    ax_etc.set_xlim(350, 970)
    ax_etc.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()

    return result
