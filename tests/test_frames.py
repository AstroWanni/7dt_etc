import pytest


@pytest.mark.parametrize('T,k,r,n,text', [
    (50, 0, 50, 1, '1x50.0s'),
    (100, 1, 0, 1, '1x100s'),
    (101, 1, 1, 2, '1x100s + 1x1.0s'),
    (300, 3, 0, 3, '3x100s'),
    (350, 3, 50, 4, '3x100s + 1x50.0s'),
])
def test_split(calc, T, k, r, n, text):
    kk, rr = calc.split(T)
    assert (float(kk), float(rr)) == pytest.approx((k, r))
    assert int(calc.n_frame(T)) == n
    assert calc.frames_text(T) == text


def test_extra_read_noise_frame(calc):
    """101 s = 100 s + 1 s: the extra readout lowers the SNR of a faint source (plain co-add)."""
    from conftest import COND
    r = calc.calculate_snr('m500', 21.0, [100, 101], verbose=False, **COND)
    assert r['m500[1]']['snr'] < r['m500[0]']['snr']
