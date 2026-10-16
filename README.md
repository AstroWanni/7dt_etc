# 7dt_etc — 7DT Exposure Time Calculator

Exposure time calculator for the **7-Dimensional Telescope (7DT)** built on an empirical
noise model of archived 7DT frames. Supports **point sources** and **extended sources
(surface brightness)**, short and long exposures, and multi-frame observations.

Developed by **GWUniverse @ Seoul National University**.

> **v0.3.0** replaces the previous depth-only calculator. See [CHANGELOG](CHANGELOG.md) for breaking changes.

---

## Installation

```bash
git clone https://github.com/AstroWanni/7dt_etc.git
cd 7dt_etc
pip install -e .
```

Requirements: `numpy`, `astropy`, `astroplan`, `matplotlib`

---

## Quick Start

```python
from seventdt import observe_and_etc

# Point source: exposure time to reach SNR = 5 for a 19.5 mag source
res, info = observe_and_etc(ra=30.968, dec=4.245, date='2026-09-14', value=19.5, snr=5)

# Extended source: 24.0 mag/arcsec^2 measured over 100 arcsec^2
res, info = observe_and_etc(ra=30.968, dec=4.245, date='2026-09-14', value=24.0, snr=5,
                            mode='extended', area=100)
```

`observe_and_etc` finds the best time of the night for the target (lowest airmass), computes the
observing conditions at that time, and returns the exposure time needed in each filter.

---

## Output of `observe_and_etc`

### Figure

![observe_and_etc example](docs/observe_and_etc_point.png)

**Top — target visibility on the given night**

| Element | Meaning |
|---|---|
| Blue line | Target altitude between astronomical evening and morning twilight |
| Orange dashed line | Moon altitude |
| Gray band | Altitude below `min_alt` (default 30°) |
| Red dashed line / dot | Best time = lowest airmass (highest altitude) of the night |
| Right axis | Airmass |
| Text box | Night date, twilight times, best time with altitude and airmass, Moon phase, Moon separation, seeing |

**Bottom — exposure time per filter at the best time**

| Element | Meaning |
|---|---|
| Points | Total exposure time needed to reach the requested SNR (m400–m875, 25 nm steps) |
| Horizontal bar | Filter bandwidth (25 nm) |
| Label above point | Total exposure time; ⚠ marks times longer than `warn_threshold` (default 1 h) |
| Dotted line | One frame (100 s) |
| Title | Source type, brightness, aperture (point) or area (extended), SNR |

### Text output

```text
Target RA, Dec = 30.9680, 4.2450 | night of 2026-09-14
Astronomical twilight (UTC): 2026-09-14 23:54 -> 2026-09-15 09:21
Above 30 deg (UTC): 2026-09-15 03:48 -> 2026-09-15 09:21 (5.6 h)
Best time (lowest airmass): 2026-09-15 07:10:50 UTC, alt = 55.2 deg, airmass = 1.22
Moon phase = 0.18, Moon sep = 162.4 deg, seeing = 2"
point source, 19.5 mag, 10" aperture, SNR = 5, frame = 100 s

filter       T[s]  frames
m400        478.4  4x100s + 1x78.4s
m425        197.8  1x100s + 1x97.8s
m450        177.5  1x100s + 1x77.5s
m475        173.6  1x100s + 1x73.6s
m500        178.6  1x100s + 1x78.6s
m525        234.1  2x100s + 1x34.1s
m550        434.4  4x100s + 1x34.4s
m575        349.6  3x100s + 1x49.6s
m600        455.7  4x100s + 1x55.7s
m625        736.8  7x100s + 1x36.8s
m650        677.1  6x100s + 1x77.1s
m675        679.3  6x100s + 1x79.3s
m700        966.4  9x100s + 1x66.4s
m725       2069.2  20x100s + 1x69.2s
m750       2875.1  28x100s + 1x75.1s
m775       6089.4  60x100s + 1x89.4s  !
m800      12790.1  127x100s + 1x90.1s  !
m825      14072.4  140x100s + 1x72.4s  !
m850      27394.5  273x100s + 1x94.5s  !
m875      56027.7  560x100s + 1x27.7s  !
```


| Line | Meaning |
|---|---|
| Astronomical twilight | Start and end of the usable night (Sun below −18°) |
| Above 30 deg | Time window when the target is above `min_alt`, and its length |
| Best time | Time of lowest airmass, used for the calculation below |
| Moon / seeing | Conditions used for the calculation |
| `T[s]` | Total exposure time to reach the requested SNR |
| `frames` | How the total is split: `k x 100s + 1 x r s` (100 s frames plus one remainder frame) |
| `!` | Longer than `warn_threshold` |

### Return values

```python
res, info = observe_and_etc(...)
res['m500']['exptime']   # total exposure time [s]
res['m500']['frames']    # e.g. '2x100s + 1x48.0s'
info['best_time_utc'], info['airmass'], info['window_utc']
```

---

## Other functions

```python
from seventdt import ExposureCalculator

calc = ExposureCalculator()              # 10" aperture model (default); ExposureCalculator(aperture=5) for 5"
cond = dict(ra=30.968, dec=4.245, obsdate='2026-09-15T07:10:00', seeing=2.0)

calc.calculate_snr('all', 19.5, exptime=300, **cond)                               # SNR for a given exposure
calc.calculate_exptime('all', 19.5, snr=10, **cond)                                # exposure for a given SNR
calc.calculate_limit('all', exptime=300, snr=5, **cond)                            # limiting magnitude
calc.calculate_limit('all', exptime=300, snr=3, mode='extended', area=100, **cond) # limiting surface brightness

calc.plot_snr_vs_exptime('all', 19.5, **cond)
calc.plot_exptime_vs_snr('all', 19.5, **cond)
calc.plot_limit_vs_exptime('all', snr=3, mode='extended', area=100, **cond)
calc.plot_noise_budget('m400', 24.0, mode='extended', area=100, **cond)
```

### Inputs

| Argument | Description |
|---|---|
| `filter` | Filter name, list, or `'all'` (g, r, i, m400–m875 in 25 nm steps) |
| `value` | AB magnitude (`mode='point'`) or surface brightness in mag/arcsec² (`mode='extended'`) |
| `exptime` | Total exposure time [s]; split into 100 s frames plus one remainder frame |
| `snr` | Target signal-to-noise ratio |
| `mode` | `'point'` (default) or `'extended'` |
| `area` | Measurement area for extended sources [arcsec²] |
| `ra`, `dec`, `obsdate` | Target and time; Moon phase, Moon separation and airmass are computed from them |
| `moon_phase`, `moon_separation`, `airmass` | Set conditions directly instead of `ra`/`dec`/`obsdate` |
| `seeing` | Seeing FWHM [arcsec] (default 2.0) |

---

## Model

For each filter, the 5σ limiting magnitude of a 100 s frame is predicted from observing conditions:

```
UL5 = b0 + b1·moon_phase + b2·moon_sep + b3·airmass + b4·seeing
```

UL5 is converted into the per-pixel background (sky + dark), and the SNR is computed with the
read noise counted once per frame:

```
SNR = S·T / sqrt( S·T + Npix · (B·T + N·RN²) )
```

| Symbol | Point source | Extended source |
|---|---|---|
| S [e⁻/s] | 10^(−0.4 (m − ZP)) | 10^(−0.4 (μ − ZP)) · A |
| Npix | pixels in the 10″ (or 5″) aperture | A / pixscale² |

`B` sky + dark rate per pixel, `T` total exposure, `N` number of frames, `RN` = 1.50 e⁻,
dark = 0.0005 e⁻/s/pix, pixscale = 0.505″/pix.

Read noise makes short exposures less efficient than the √t scaling of the previous version,
and the SNR of bright sources now includes the source photon noise.

---

## Scope and limitations

- **Camera setting**: gain 2750 (high gain), 1×1 binning only.
- **Training data**: 2024–2026 100 s frames. Predictions are most reliable for typical conditions
  (airmass 1.2–1.7, average Moon). Very bright Moon and airmass > 1.7 are less accurate.
- **Filters**: `'all'` returns g, r, i and m400–m875 (25 nm steps). Filters with few frames or a different
  bandwidth (z, m438, m483, m512, m534, m561, m586, m615, m640, m661, m769w, m832w) are kept in the model
  files but excluded by default; request them by name to use them (a warning is printed).
- **Extended sources over large areas**: pixel noise is assumed independent, so limits over very large areas
  are optimistic (flat-field and background residuals are not included).
