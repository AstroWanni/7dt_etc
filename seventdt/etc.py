"""
Function interface (one call, no calculator object).

    from seventdt import calculate_exptime, calculate_snr, calculate_limit
    calculate_exptime(19.5, snr=10, filters='all', ra=30.97, dec=4.25, obsdate='2026-09-15T07:00:00')
    calculate_limit(300, snr=3, mode='extended', area=100, moon_phase=0.1, moon_separation=120, airmass=1.2)
"""
from functools import lru_cache
from . import constants as C
from .calculator import ExposureCalculator


@lru_cache(maxsize=None)
def _calc(aperture):
    return ExposureCalculator(aperture=aperture)


def calculate_exptime(value, snr=5.0, filters='all', mode='point', area=1.0, aperture=C.DEFAULT_APERTURE,
                      verbose=True, **conditions):
    """Total exposure [s] to reach 'snr'. value: AB mag (point) or mag/arcsec^2 (extended)."""
    return _calc(aperture).calculate_exptime(filters, value, snr=snr, mode=mode, area=area, verbose=verbose, **conditions)


def calculate_snr(value, exptime, filters='all', mode='point', area=1.0, aperture=C.DEFAULT_APERTURE,
                  verbose=True, **conditions):
    """SNR for total exposure 'exptime' [s]."""
    return _calc(aperture).calculate_snr(filters, value, exptime, mode=mode, area=area, verbose=verbose, **conditions)


def calculate_limit(exptime, snr=5.0, filters='all', mode='point', area=1.0, aperture=C.DEFAULT_APERTURE,
                    verbose=True, **conditions):
    """Limiting magnitude (point) or surface brightness (extended) for total exposure 'exptime' [s]."""
    return _calc(aperture).calculate_limit(filters, exptime, snr=snr, mode=mode, area=area, verbose=verbose, **conditions)
