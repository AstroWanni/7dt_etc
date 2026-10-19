import pytest


@pytest.mark.parametrize('T,n', [(0.5, 1), (50, 1), (100, 1), (101, 2), (300, 3), (350, 4)])
def test_n_frame(calc, T, n):
    assert int(calc.n_frame(T)) == n


def test_extra_read_noise_frame(calc):
    """101 s = 2 frames: the extra readout lowers the SNR of a faint source (plain co-add)."""
    from conftest import COND
    r = calc.calculate_snr('m500', 21.0, [100, 101], verbose=False, **COND)
    assert r['m500[1]']['snr'] < r['m500[0]']['snr']


def test_exptime_result_keys(calc):
    from conftest import COND
    r = calc.calculate_exptime('m500', 20.0, snr=10, verbose=False, **COND)['m500']
    assert r['n_frame'] == int(calc.n_frame(r['exptime']))
    assert 'frames' not in r and r['skysig_adu'] == pytest.approx(r['skysig_e'] / 0.256)
