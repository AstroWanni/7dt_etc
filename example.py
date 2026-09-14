"""
7DTC usage examples.
Run from the repo root: python example.py
"""
from astropy.time import Time
import seventdt as dtc
from seventdt.plot import plot_exptime_per_filter
from seventdt.visibility import observe_and_etc

# 1. Point source — exposure time
result = dtc.calculate_exptime(
    magnitude=19.5, snr=10, filters='all',
    seeing=1.5, ra=30.968, dec=4.245,
    obsdate=Time('2026-09-14T04:00:00'),
)
for filt, r in result.items():
    print(f"  {filt:6}  t_exp = {r['t_exp']:8.1f} s   UL5@100s = {r['ul5_ref']:.2f}")

# 2. Point source — SNR
result_snr = dtc.calculate_snr(
    magnitude=19.5, exptime=300,
    filters=['g', 'r', 'i'], seeing=1.5,
)
for filt, r in result_snr.items():
    print(f"  {filt:6}  SNR = {r['snr']:.1f}")

# 3. Point source — plot
result_all = dtc.calculate_exptime(magnitude=19.5, snr=10, filters='all', seeing=2.0)
plot_exptime_per_filter(result_all, magnitude=19.5, snr=10,
                        save_path='exptime_point_source.png')

# 4. Visibility + ETC
observe_and_etc(
    ra=30.968, dec=4.245,
    date='2026-09-14',
    magnitude=19.5, snr=10,
    filters='all', seeing=1.5,
    save_path='visibility_etc.png',
)
