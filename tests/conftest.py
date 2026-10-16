import pytest
from seventdt import ExposureCalculator

# Fixed conditions (no ra/dec/obsdate -> no ephemeris calculation, tests run offline)
COND = dict(moon_phase=0.5, moon_separation=90.0, airmass=1.2, seeing=2.0)
DARK = dict(moon_phase=0.0, moon_separation=120.0, airmass=1.1, seeing=2.0)


@pytest.fixture(scope='session')
def calc():
    return ExposureCalculator()


@pytest.fixture(scope='session')
def calc5():
    return ExposureCalculator(aperture=5)
