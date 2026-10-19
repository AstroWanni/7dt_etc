import pytest
from sevendt_etc import DEFAULT_BANDS, EXTRA_BANDS


def test_all_is_default_only(calc):
    bands = calc.filter_list('all')
    assert bands == [f for f in DEFAULT_BANDS if f in calc.depth]
    assert not set(bands) & set(EXTRA_BANDS)


def test_extra_band_warns(calc):
    extra = [f for f in EXTRA_BANDS if f in calc.depth]
    if not extra:
        pytest.skip('no extra bands in the model')
    with pytest.warns(UserWarning):
        calc.filter_list(extra[0])


def test_unknown_band_raises(calc):
    with pytest.raises(KeyError):
        calc.filter_list('m999')


def test_available_superset(calc):
    assert set(calc.filter_list('all')) <= set(calc.filter_list('available'))
