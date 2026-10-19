"""
sevendt_etc usage examples (run from the repo root: python example.py).
See examples/7dt_etc_demo.ipynb for the full walkthrough.
"""
from sevendt_etc import ExposureCalculator, observe_and_etc

calc = ExposureCalculator()
cond = dict(moon_phase=0.2, moon_separation=120, airmass=1.2, seeing=2.0)

# 1. Exposure time per band for a 19.5 mag point source at SNR = 10 (main function)
res = calc.calculate_exptime('all', 19.5, snr=10, **cond)

# 2. Recommended time on a night + exposure per band (text report + figure)
res, info = observe_and_etc(ra=30.968, dec=-4.245, date='2026-09-14', value=19.5, snr=5,
                            save_path='observe_and_etc.png', report_path='observe_and_etc.txt')
