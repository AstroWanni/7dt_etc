# Changelog

## [1.0.0] — unreleased (pre-release updates are kept under 1.0.0)

First public release. The calculator is built on an empirical noise model of archived 7DT frames.
Earlier versions (0.x) were development snapshots and are not supported.

### Added
- Package renamed `seventdt` → `sevendt_etc` (`import sevendt_etc`).
- Text report for `calculate_exptime` / `observe_and_etc` (detector, source, conditions, observation plan,
  per-band sky [mag/arcsec², e⁻/s/pix], skysig [e⁻, ADU], exposure time and number of frames);
  `observe_and_etc(..., report_path=...)` saves it, `info['report']` returns it.
- Plot functions do not return the Figure by default (each plot is shown once in Jupyter); `return_fig=True` returns it.
- **Noise model**: SNR = S·T / sqrt(S·T + Npix·(B·T + N·RN²)), with the per-pixel background B derived from the
  depth model and the read noise counted once per frame (RN = 1.50 e⁻, dark = 0.0005 e⁻/s/pix, gain 0.256 e⁻/ADU,
  measured from 2026 bias / dark / flat frames, GAIN 2750, 1×1).
- **Extended sources**: `mode='extended'` with surface brightness [mag/arcsec²] and measurement `area` [arcsec²].
- **Multi-frame exposures**: totals longer than 100 s are taken as N = ceil(T / 100 s) frames;
  the number of frames (`n_frame`) is reported in every result.
- **New depth models** (2024–2026, 100 s, GAIN 2750, 1×1), one per aperture:
  `models/depth_model_ul5_5.json` (10″, default) and `models/depth_model_ul5_4.json` (5″);
  each band carries its own `zp_e_s` and aperture.
- `ExposureCalculator(aperture=5)` for the 5″ model.
- `calculate_limit`: limiting magnitude or limiting surface brightness.
- Plots: `plot_snr_vs_exptime`, `plot_exptime_vs_snr`, `plot_limit_vs_exptime`, `plot_noise_budget`.
- `observe_and_etc`: target visibility on a night + exposure per filter at the lowest-airmass time.
- Airmass is computed from `ra`, `dec`, `obsdate` (or given directly).
- Tests (`tests/`).

### Changed
- **Single exposure vs 100 s stack**: `calculate_exptime` reports both `T_single` (`exptime_single`, one exposure)
  and `T_stack` (`exptime` = `n_frame` × 100 s, SNR = √N × one-frame SNR); the text report gets a `T_single` column.
- Auxiliary plots: `plot_noise_budget` and `plot_limit_vs_exptime` show a single exposure (1–1000 s);
  `plot_snr_vs_exptime` / `plot_exptime_vs_snr` draw the single exposure as lines and stacked 100 s frames as points
  (100–10000 s). The √t reference line was removed.
- Depth model features: `moon_phase, moon_sep, airmass, seeing` (was `moon_phase, moon_sep, hour_since_sunset, seeing`).
- `filter='all'` returns g, r, i and m400–m875 (25 nm steps). z, m438, m483, m512, m534, m561, m586, m615, m640, m661,
  m769w, m832w remain in the model files and can be requested by name (a warning is printed);
  `filter='available'` lists every band.
- Medium-band plot colors use a fixed rainbow scale over 400–875 nm.
- Packaging moved from `setup.py` to `pyproject.toml`.

### API changes from the development snapshots (0.2.x)
| 0.2.x | 1.0.0 |
|---|---|
| `calculate_exptime(filter, magnitude, snr, ...)` | `calculate_exptime(filter, value, snr, mode='point', area=1.0, ...)` |
| `calculate_snr(filter, magnitude, exptime, ...)` | `calculate_snr(filter, value, exptime, mode='point', area=1.0, ...)` |
| result key `t_exp` | `exptime` (+ `n_frame`) |
| `hour_since_sunset` condition | removed; `airmass` added |
| `calc.load_model()` / `update_depth_model()` / `fit_depth_model()` | removed; models are packaged (refit with the analysis scripts) |
| `observe_and_etc(..., magnitude=...)` | `observe_and_etc(..., value=..., mode=..., area=...)` |

### Removed
- `sevendt_etc/core.py` (previous calculator) and `depth_model_20260329052547.json`.

### Known limitations
- GAIN 2750, 1×1 binning only.
- Most reliable for typical conditions (airmass 1.2–1.7, average Moon); very bright Moon and airmass > 1.7 are less accurate.
- Extended-source limits over large areas assume independent pixel noise and are optimistic.

## Development snapshots (0.x, not released)
- Point-source calculator based on the 100 s 5σ depth model with sqrt(t) scaling.
