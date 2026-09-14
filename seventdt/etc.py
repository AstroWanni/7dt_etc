"""
High-level ETC functions for point source.

    from seventdt import calculate_exptime, calculate_snr
"""
from astropy.time import Time
from .core import ExposureCalculator


def _make_calc() -> ExposureCalculator:
    calc = ExposureCalculator()
    calc.load_model()
    return calc


def calculate_exptime(
    magnitude: float,
    snr: float,
    filters: str | list = 'all',
    ra: float = None,
    dec: float = None,
    obsdate: str | Time = None,
    moon_phase: float = None,
    moon_separation: float = None,
    hour_since_sunset: float = None,
    seeing: float = None,
    verbose: bool = True,
) -> dict:
    """
    Required exposure time for a point source.

    Parameters
    ----------
    magnitude  : AB magnitude of the source
    snr        : target SNR
    filters    : filter name, list of names, or 'all'
    ra, dec    : target coordinates [deg] — used to compute moon conditions
    obsdate    : observation time (ISO string or astropy Time)
    moon_phase, moon_separation, hour_since_sunset : manual condition override
    seeing     : seeing FWHM [arcsec] (default: 2.0")
    verbose    : print computed conditions

    Returns
    -------
    dict  keyed by filter name, each value contains:
        t_exp, magnitude, snr, ul5_ref,
        moon_phase, moon_separation, hour_since_sunset, seeing
    """
    if isinstance(obsdate, str):
        obsdate = Time(obsdate)
    return _make_calc().calculate_exptime(
        filter=filters, magnitude=magnitude, snr=snr,
        ra=ra, dec=dec, obsdate=obsdate,
        moon_phase=moon_phase, moon_separation=moon_separation,
        hour_since_sunset=hour_since_sunset, seeing=seeing,
        verbose=verbose,
    )


def calculate_snr(
    magnitude: float,
    exptime: float,
    filters: str | list = 'all',
    ra: float = None,
    dec: float = None,
    obsdate: str | Time = None,
    moon_phase: float = None,
    moon_separation: float = None,
    hour_since_sunset: float = None,
    seeing: float = None,
    verbose: bool = True,
) -> dict:
    """
    Expected SNR for a point source in a given exposure time.

    Parameters
    ----------
    magnitude : AB magnitude of the source
    exptime   : exposure time [s]
    (rest identical to calculate_exptime)

    Returns
    -------
    dict  keyed by filter name, each value contains:
        snr, magnitude, exptime, ul5_exptime, filter,
        moon_phase, moon_separation, hour_since_sunset, seeing
    """
    if isinstance(obsdate, str):
        obsdate = Time(obsdate)
    return _make_calc().calculate_snr(
        filter=filters, magnitude=magnitude, exptime=exptime,
        ra=ra, dec=dec, obsdate=obsdate,
        moon_phase=moon_phase, moon_separation=moon_separation,
        hour_since_sunset=hour_since_sunset, seeing=seeing,
        verbose=verbose,
    )
