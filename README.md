# 7dt_etc — 7DT Empirical Exposure Time Calculator

Empirical exposure time calculator for the **7-Dimensional Telescope (7DT)**, fitted on
5σ limiting magnitudes (UL5) of archived 100 s frames.

Developed by **GWUniverse @ Seoul National University**.

---

## Model

For each of 24 filters (g, r, i, z, m400–m875), UL5 at the reference exposure time
(100 s) is predicted from observing conditions via linear regression:

```
UL5_ref = b0 + b1·moon_phase + b2·moon_sep + b3·hour_since_sunset + b4·seeing
```

Exposure time and SNR are then computed from:

```
t_exp = t_ref · 10^((mag − UL5_ref) / 1.25) · (SNR / 5)²
SNR   = 5 · 10^(−0.4 · (mag − UL5_t))        where UL5_t = UL5_ref + 1.25·log10(t / t_ref)
```

UL5 was measured using **SExtractor fixed-aperture photometry (10″ diameter)**.

---

## Installation

```bash
git clone https://github.com/AstroWanni/7dt_etc.git
cd 7dt_etc
pip install -e .
```

### Requirements

`numpy`, `astropy`, `astroplan`, `matplotlib`, `tqdm`  
`scikit-learn` is only needed to refit the depth model (`ExposureCalculator.fit_depth_model`).

---

## Quick Start

```python
from seventdt.visibility import observe_and_etc

result = observe_and_etc(
    ra=30.968, dec=4.245,
    date='2026-09-14',
    magnitude=19.5, snr=10,
    filters='all', seeing=1.5,
)
```

---

## Observing Conditions

| Input | Computed |
|---|---|
| `ra` + `dec` + `obsdate` | moon_phase, moon_separation, hour_since_sunset |
| `obsdate` only | moon_phase, hour_since_sunset |
| nothing | DEFAULT_CONDITIONS (moon_phase 0.5, moon_sep 90°, 5 h after sunset) |

`seeing` is never computed automatically — pass it explicitly.  
Any condition can also be given directly (`moon_phase=`, `moon_separation=`, etc.).

---

## Package Structure

```
seventdt/
├── core.py          # ExposureCalculator (original standalone calculator)
├── constants.py     # instrument constants, filter wavelengths
├── etc.py           # calculate_exptime / calculate_snr
├── plot.py          # plot_exptime_per_filter
├── visibility.py    # observe_and_etc (altitude curve + ETC)
└── depth_model_*.json
```

---

## Refitting the Model

```python
from seventdt import ExposureCalculator
calc = ExposureCalculator()
calc.update_depth_model(depth_table='my_depth_data.csv')
```

Input table must contain: `filter, exptime, ul5, seeing, ra_center, dec_center, mjd`.  
`scikit-learn` is required for fitting.
