"""Observing conditions at the 7DT site: Moon phase, Moon separation, airmass."""
import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time
from .constants import OBSERVER_ELEVATION, OBSERVER_LATITUDE, OBSERVER_LONGITUDE, OBSERVER_NAME, OBSERVER_TIMEZONE

DEFAULT_CONDITIONS = {'moon_phase': 0.5, 'moon_separation': 90.0, 'airmass': 1.2, 'seeing': 2.0}


def get_location():
    return EarthLocation.from_geodetic(lat=OBSERVER_LATITUDE, lon=OBSERVER_LONGITUDE, height=OBSERVER_ELEVATION)


def get_observer():
    return Observer(location=get_location(), name=OBSERVER_NAME, timezone=OBSERVER_TIMEZONE)


def moon_phase(obsdate):
    """Illuminated fraction of the Moon (0 = new, 1 = full)."""
    return float(get_observer().moon_illumination(Time(obsdate)))


def moon_separation(ra, dec, obsdate):
    """Moon-target separation [deg]."""
    moon = get_body('moon', Time(obsdate))   # geocentric (GCRS) Moon position
    target = SkyCoord(ra, dec, unit=(u.deg, u.deg)).transform_to(moon.frame)
    return float(moon.separation(target).deg)


def airmass(ra, dec, obsdate):
    """sec(z) of the target at the 7DT site; NaN if below the horizon."""
    alt = SkyCoord(ra, dec, unit=(u.deg, u.deg)).transform_to(AltAz(obstime=Time(obsdate), location=get_location())).alt.deg
    return float(1 / np.sin(np.radians(alt))) if alt > 0 else np.nan


def resolve_conditions(ra=None, dec=None, obsdate=None, moon_phase=None, moon_separation=None, airmass=None, seeing=None):
    """
    Conditions for the depth model.
      - ra, dec, obsdate given : Moon phase, Moon separation and airmass are computed (given values are overridden)
      - only obsdate given     : Moon phase is computed
      - ra, dec without obsdate: obsdate = now
      - seeing is never computed; missing values fall back to DEFAULT_CONDITIONS
    Returns dict(moon_phase, moon_separation, airmass, seeing, ra, dec, obsdate)
    """
    has_target = ra is not None and dec is not None
    if has_target and obsdate is None:
        obsdate = Time.now()
    if obsdate is not None:
        obsdate = Time(obsdate)
        moon_phase = _moon_phase(obsdate)
        if has_target:
            moon_separation = _moon_separation(ra, dec, obsdate)
            airmass = _airmass(ra, dec, obsdate)
            if not np.isfinite(airmass):
                print(f'WARNING: target below the horizon at {obsdate.iso} (airmass = NaN)')
    given = {'moon_phase': moon_phase, 'moon_separation': moon_separation, 'airmass': airmass, 'seeing': seeing}
    cond = {k: float(DEFAULT_CONDITIONS[k] if v is None else v) for k, v in given.items()}
    cond.update(ra=ra, dec=dec, obsdate=obsdate)
    return cond


_moon_phase, _moon_separation, _airmass = moon_phase, moon_separation, airmass
