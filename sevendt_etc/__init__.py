"""
7dt_etc — 7DT exposure time calculator with an empirical noise model
(point and extended sources; GAIN 2750, 1x1 binning).
"""
from .calculator import ExposureCalculator, load_depth_model
from .conditions import get_observer, resolve_conditions
from .constants import DEFAULT_BANDS, EXTRA_BANDS, MEDIUM_BANDS
from .etc import calculate_exptime, calculate_limit, calculate_snr
from .visibility import observe_and_etc

__version__ = '1.0.0'
__all__ = ['ExposureCalculator', 'load_depth_model', 'get_observer', 'resolve_conditions', 'calculate_exptime',
           'calculate_snr', 'calculate_limit', 'observe_and_etc', 'DEFAULT_BANDS', 'EXTRA_BANDS', 'MEDIUM_BANDS']
