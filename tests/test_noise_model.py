import numpy as np
import pytest
from conftest import COND, DARK


def test_sky_limited_sqrt_t(calc):
    """Long exposures of a faint source in whole frames: SNR grows as sqrt(T)."""
    r = calc.calculate_snr('m500', 25.0, [900, 3600], verbose=False, **COND)
    ratio = r['m500[1]']['snr'] / r['m500[0]']['snr']
    assert ratio == pytest.approx(2.0, rel=1e-3)


def test_read_noise_limited_linear(calc):
    """Very short exposures of a faint source in a dark sky: SNR grows ~linearly with T."""
    r = calc.calculate_snr('m400', 25.0, [0.1, 0.2], verbose=False, **DARK)
    ratio = r['m400[1]']['snr'] / r['m400[0]']['snr']
    assert 1.9 < ratio <= 2.0


def test_bright_source_photon_noise(calc):
    """Bright source: SNR ~ sqrt(source counts), far below the sqrt(t)-scaled UL5 extrapolation."""
    r = calc.calculate_snr('m500', 12.0, 100, verbose=False, **COND)['m500']
    assert r['frac_source'] > 0.9
    assert r['snr'] < r['snr_sqrt_t']


@pytest.mark.parametrize('mode,area,value', [('point', 1.0, 19.5), ('extended', 100.0, 24.0)])
@pytest.mark.parametrize('snr', [3, 10, 50])
def test_exptime_roundtrip(calc, mode, area, value, snr):
    t = calc.calculate_exptime('m600', value, snr=snr, mode=mode, area=area, verbose=False, **COND)['m600']['exptime']
    s = calc.calculate_snr('m600', value, t, mode=mode, area=area, verbose=False, **COND)['m600']['snr']
    assert s == pytest.approx(snr, rel=1e-6)


@pytest.mark.parametrize('mode,area', [('point', 1.0), ('extended', 100.0)])
def test_limit_roundtrip(calc, mode, area):
    lim = calc.calculate_limit('m700', 300, snr=5, mode=mode, area=area, verbose=False, **COND)['m700']['limit']
    s = calc.calculate_snr('m700', lim, 300, mode=mode, area=area, verbose=False, **COND)['m700']['snr']
    assert s == pytest.approx(5.0, rel=1e-6)


def test_extended_area_scaling(calc):
    """Sky-limited extended source: 4x area -> 2x SNR."""
    r1 = calc.calculate_snr('m500', 26.0, 3600, mode='extended', area=25, verbose=False, **COND)['m500']['snr']
    r4 = calc.calculate_snr('m500', 26.0, 3600, mode='extended', area=100, verbose=False, **COND)['m500']['snr']
    assert r4 / r1 == pytest.approx(2.0, rel=0.02)


def test_variance_fractions_sum_to_one(calc):
    r = calc.calculate_snr('m500', 20.0, 150, verbose=False, **COND)['m500']
    assert sum(r[k] for k in r if k.startswith('frac_')) == pytest.approx(1.0)


def test_fainter_needs_longer(calc):
    t = [calc.calculate_exptime('m500', m, snr=5, verbose=False, **COND)['m500']['exptime'] for m in (18, 19, 20)]
    assert np.all(np.diff(t) > 0)
