"""7DT site, detector and filter constants (GAIN 2750, 1x1 binning)."""
import re

# Site (El Sauce, Chile)
OBSERVER_NAME = '7DT'
OBSERVER_LATITUDE = -30.4704      # [deg]
OBSERVER_LONGITUDE = -70.7804     # [deg]
OBSERVER_ELEVATION = 1580         # [m]
OBSERVER_TIMEZONE = 'America/Santiago'

# Detector (measured from 2026 bias / dark / flat, GAIN 2750, 1x1)
GAIN_MODE = 2750
EGAIN = 0.256          # [e-/ADU] (flat photon transfer curve)
READ_NOISE = 1.50      # [e-] (bias, RMS)
DARK_CURRENT = 0.0005  # [e-/s/pix] (100 s master dark)
PIXSCALE = 0.505       # [arcsec/pix]
FRAME_TIME = 100.0     # [s] single-frame exposure

# Depth models shipped with the package (models/depth_model_ul5_<n>.json)
APERTURE_MODELS = {10: 'depth_model_ul5_5.json', 5: 'depth_model_ul5_4.json'}   # diameter ["] -> file
DEFAULT_APERTURE = 10

# Filters
BROAD_BANDS = ['g', 'r', 'i']
MEDIUM_BANDS = [f'm{w}' for w in range(400, 876, 25)]
DEFAULT_BANDS = BROAD_BANDS + MEDIUM_BANDS          # returned by filter='all'
EXTRA_BANDS = ['z', 'm438', 'm483', 'm512', 'm534', 'm561', 'm586', 'm615', 'm640', 'm661', 'm769w', 'm832w']
EXTRA_BANDS_NOTE = 'few training frames or a different bandwidth; use with caution'
BAND_WIDTH = {**{f: 25.0 for f in MEDIUM_BANDS}, 'g': 130.0, 'r': 130.0, 'i': 130.0, 'z': 100.0}   # [nm]
BROAD_WAVE = {'g': 475, 'r': 620, 'i': 760, 'z': 900}   # [nm]
BROAD_COLORS = {'g': 'green', 'r': 'red', 'i': 'gold', 'z': 'k'}
MEDIUM_RANGE = (400, 875)   # [nm] fixed rainbow color scale for medium bands


def wavelength(filt):
    """Central wavelength [nm] of a filter name (m400 -> 400, g -> 475)."""
    m = re.match(r'm(\d+)', filt)
    return int(m.group(1)) if m else BROAD_WAVE.get(filt)


def band_key(filt):
    """Sort key: broad bands first (g, r, i, z), then medium bands by wavelength."""
    order = ['g', 'r', 'i', 'z']
    return (0, order.index(filt), 0) if filt in order else (1, 0, wavelength(filt) or 10 ** 9)
