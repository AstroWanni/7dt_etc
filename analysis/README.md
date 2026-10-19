# analysis — how the packaged models were built

These scripts rebuild the detector noise constants and the depth models shipped in `sevendt_etc/models/`.
They are **not** part of the installed package: they read the 7DT archive and database directly, so paths
at the top of each script point to the GWUniverse servers (`/lyman`, `/data3`). Edit them for your environment.

Scope: GAIN 2750 (high gain), 1×1 binning, `processing_version = 'v2'`.

| Step | Script | Input | Output | Result used in the package |
|---|---|---|---|---|
| 0 | `00_query_frames.sql` | `gwu` database | `noise_v2_frames.csv` (one row per processed frame) | input of steps 2–3 |
| 1 | `01_detector_noise_raw.py` | raw BIAS / DARK / FLAT frames (`/lyman/data1/obsdata`) | `detector_2026.csv`, `flat_ptc_2026.csv` | `READ_NOISE = 1.50 e⁻`, `DARK_CURRENT = 0.0005 e⁻/s/pix`, `EGAIN = 0.256 e⁻/ADU` (`constants.py`) |
| 2 | `02_skysig_vs_model.py` | `noise_v2_frames.csv` | `skysig_vs_model.png`, per-exposure summary | validation: predicted pixel noise sqrt(sky·t + dark·t + RN²) vs header SKYSIG, 3–120 s |
| 3 | `03_fit_depth_model.py` | `noise_v2_frames.csv` | `depth_model_ul5_4_<stamp>.json`, `depth_model_ul5_5_<stamp>.json`, validation CSV | `sevendt_etc/models/depth_model_ul5_4.json`, `depth_model_ul5_5.json` |
| 4 | `04_validate_snr_T02524.py` | processed frames + catalogs of one tile (`/lyman/data2/processed`) | empty-aperture and repeatability tables | validation: aperture noise and repeated-measurement SNR vs the noise model |

## Updating the models

```bash
psql -U gwuteam -p 5433 gwu --csv -f 00_query_frames.sql -o /data3/cwkim/exptimecal/noise_v2_frames.csv
python 03_fit_depth_model.py
cp <OUT_DIR>/depth_model_ul5_5_<stamp>.json ../sevendt_etc/models/depth_model_ul5_5.json
cp <OUT_DIR>/depth_model_ul5_4_<stamp>.json ../sevendt_etc/models/depth_model_ul5_4.json
python -m pytest -q ../tests
```

Record the model stamp and the training period in `CHANGELOG.md` when the models are replaced.

## Notes

- Step 1 reads raw pixels from a mounted remote disk; it copies only the `.head` files locally and limits
  parallel reads. Run it on a few nights per unit when re-measuring the detector constants.
- Step 3 splits nights 80 / 20 into train / test and fits with iterative 3σ clipping; filters with fewer than
  200 training frames are skipped.
- Frames on a few nights with abnormal background noise (e.g. 2025-06-10, 2024-08-09/10) show SKYSIG well above
  the model; they are removed by the clipping in step 3.
