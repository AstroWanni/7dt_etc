import numpy as np
from seventdt import DEFAULT_BANDS, load_depth_model
from conftest import COND

KEYS = {'intercept', 'coefficients', 'features', 'exptime', 'aperture', 'zp_e_s'}


def test_models_load():
    for ap in (10, 5):
        m = load_depth_model(ap)
        assert m, f'empty model for {ap}"'
        for f, d in m.items():
            assert KEYS <= set(d), f'{ap}" {f}: missing {KEYS - set(d)}'
            assert d['aperture'] == ap
            assert len(d['coefficients']) == len(d['features'])


def test_default_bands_present(calc):
    missing = [f for f in DEFAULT_BANDS if f not in calc.depth]
    assert not missing, f'default bands missing from the 10" model: {missing}'


def test_sky_consistent_between_apertures(calc, calc5):
    """Sky is a property of the sky, not of the aperture: both models should agree at typical conditions."""
    from seventdt import resolve_conditions
    cond = resolve_conditions(**COND)
    d = [calc.band_params(f, cond)['mu_sky'] - calc5.band_params(f, cond)['mu_sky'] for f in calc.default_bands]
    assert abs(np.nanmedian(d)) < 0.2
