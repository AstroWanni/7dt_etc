#%%
"""
7DT detector / sky noise characterization (2026 raw data)
  0) Copy .head files from <ROOT>/7DTxx/2026-*/ (mounted remote server) to <HEAD_DIR>/7DTxx/<night>/
     (only new or changed files; low parallelism to limit remote IO)
  1) Index the local .head copies
  2) Per night: read noise (BIAS), dark current + hot pixels (DARK), PTC gain (FLAT pairs)
  3) Per science frame: sky level [e-/s/pix] and empty-sky noise check
     sigma_meas (tile mad_std) vs sigma_pred = sqrt(sky + D*t + RN^2 + calibration noise)
Outputs (OUT_DIR): index_2026.csv, detector_2026.csv, flat_ptc_2026.csv, sky_2026.csv
All pixel statistics use a central CROP x CROP region.
"""
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import mad_std
from tqdm import tqdm

ROOT = Path('/lyman/data1/obsdata')          # remote (mounted) raw data
HEAD_DIR = Path('/data3/cwkim/exptimecal/exposurecalculator/exposurecalculator/update_261006/header')  # local .head copies
UNITS = 'all'                 # 'all' or e.g. ['7DT01', '7DT02']
YEAR = '2026'
OUT_DIR = Path('./noise_2026')
CROP = 2048                   # central crop size [pix]
TILE = 64                     # tile size [pix] for local sky statistics
MAX_SCI_PER_FILTER = 10       # science frames per (night, filter), evenly spaced in time; None = all
HOT_THRESH = 0.5              # [e-/s/pix] dark current above this = hot pixel
PTC_RANGE = (1000, 40000)     # [ADU, bias-subtracted] usable flat level for PTC
N_PROC = 16
N_COPY = 4                    # parallel threads for copying .head files from the remote server
CAL_TYPES = ('BIAS', 'DARK', 'FLAT')
KEYS = ['IMAGETYP', 'OBJECT', 'FILTER', 'EXPTIME', 'GAIN', 'EGAIN', 'CCD-TEMP', 'XBINNING', 'TELESCOP',
        'DATE-OBS', 'MJD', 'AIRMASS', 'MOONPHAS', 'MOONSEP', 'CLUDFRAC', 'RA', 'DEC']
NUM_KEYS = ['EXPTIME', 'GAIN', 'EGAIN', 'CCD-TEMP', 'XBINNING', 'MJD', 'AIRMASS', 'MOONPHAS', 'MOONSEP', 'CLUDFRAC', 'RA', 'DEC']

#%% 0) Copy .head files to local disk
def night_dirs():
    units = sorted(p for p in ROOT.glob('7DT[0-9][0-9]') if p.is_dir())
    if UNITS != 'all':
        units = [p for p in units if p.name in UNITS]
    return sorted(d for u in units for d in u.glob(f'{YEAR}-*') if d.is_dir())

def sync_heads(night_dir):
    """Copy new/changed .head files of one remote night folder; return the local folder."""
    dest = HEAD_DIR / night_dir.parent.name / night_dir.name
    dest.mkdir(parents=True, exist_ok=True)
    for hp in night_dir.glob('*.head'):
        lp = dest / hp.name
        st = hp.stat()
        if lp.exists() and lp.stat().st_size == st.st_size and lp.stat().st_mtime >= st.st_mtime:
            continue
        shutil.copy2(hp, lp)
    return dest

#%% 1) Index (local .head copies)
def read_head(path, fits_path):
    try:
        return fits.Header.fromtextfile(str(path))
    except Exception:
        return fits.getheader(fits_path)      # fallback: read header from the remote FITS

def index_night(local_dir):
    unit, night = local_dir.parent.name, local_dir.name
    m = re.search(r'gain(\d+)', night)
    rows = []
    for hp in sorted(local_dir.glob('*.head')):
        fits_path = str(ROOT / unit / night / hp.with_suffix('.fits').name)
        try:
            hdr = read_head(hp, fits_path)
        except Exception:
            continue
        row = {'unit': unit, 'night': night, 'gain_dir': int(m.group(1)) if m else np.nan, 'path': fits_path}
        for k in KEYS:
            v = hdr.get(k, np.nan)
            row[k] = np.nan if isinstance(v, fits.card.Undefined) else v
        rows.append(row)
    return rows

#%% 2-3) Per-night analysis
def load_crop(path, size=CROP):
    with fits.open(path) as hdul:          # .section reads only the crop (works with BZERO)
        hdu = hdul[0]
        ny, nx = hdu.shape
        y0, x0 = (ny - size) // 2, (nx - size) // 2
        return np.asarray(hdu.section[y0:y0 + size, x0:x0 + size], dtype=np.float32)

def tile_stats(img, tile=TILE):
    n = img.shape[0] // tile
    arr = img[:n * tile, :n * tile].reshape(n, tile, n, tile).swapaxes(1, 2).reshape(n, n, -1)
    return np.median(arr, axis=-1), mad_std(arr, axis=-1)

def analyze_night(records):
    df = pd.DataFrame(records)
    typ = df['IMAGETYP'].astype(str).str.strip().str.upper()
    bias, dark, flat = df[typ == 'BIAS'], df[typ == 'DARK'], df[typ == 'FLAT']
    sci = df[~typ.isin(CAL_TYPES) & (df['EXPTIME'] > 0)]
    g = float(np.nanmedian(df['EGAIN']))
    info = {'unit': df['unit'].iloc[0], 'night': df['night'].iloc[0], 'gain_set': np.nanmedian(df['GAIN']),
            'egain': g, 'ccd_temp': np.nanmedian(df['CCD-TEMP']),
            'n_bias': len(bias), 'n_dark': len(dark), 'n_flat': len(flat), 'n_sci': len(sci)}
    if len(bias) < 2:
        return info, [], []
    # --- BIAS: read noise
    B = np.stack([load_crop(p) for p in bias['path']])
    lev = np.median(B, axis=(1, 2))
    mbias = np.median(B, axis=0)
    Bc = B - lev[:, None, None]
    sig_pix = np.std(Bc, axis=0, ddof=1)
    rn_pair = mad_std(Bc[0] - Bc[1]) / np.sqrt(2) * g
    info.update(bias_level_adu=float(np.median(lev)), bias_level_std_adu=float(np.std(lev)),
                rn_rms_e=float(np.sqrt(np.mean(sig_pix ** 2)) * g), rn_med_e=float(np.median(sig_pix) * g),
                rn_pair_e=float(rn_pair))
    del B, Bc, sig_pix
    nb = len(bias)
    # --- DARK: dark current, hot pixels
    dark_rate_adu = np.zeros_like(mbias)
    t_d, nd, D_med = np.nan, 0, 0.0
    if len(dark) >= 2:
        t_d = float(dark['EXPTIME'].mode().iloc[0])
        dd = dark[dark['EXPTIME'] == t_d]
        nd = len(dd)
        if nd >= 2:
            Dk = np.stack([load_crop(p) for p in dd['path']])
            dark_rate_adu = (np.median(Dk, axis=0) - mbias) / t_d
            dmap_e = dark_rate_adu * g
            D_med = float(np.median(dmap_e))
            v_pair = (mad_std(Dk[0] - Dk[1]) * g) ** 2 / 2
            info.update(dark_exptime=t_d, n_dark_used=nd, dark_med_e_s=D_med,
                        dark_mean_e_s=float(np.mean(dmap_e)), dark_var_e_s=float((v_pair - rn_pair ** 2) / t_d),
                        hot_frac=float(np.mean(dmap_e > HOT_THRESH)))
            del Dk, dmap_e
    hot = dark_rate_adu * g > HOT_THRESH
    # --- FLAT: photon transfer (consecutive pairs per filter)
    rn_adu = rn_pair / g
    ptc_rows = []
    for filt, fg in flat.sort_values('MJD').groupby('FILTER'):
        fg = fg.reset_index(drop=True)
        for i in range(len(fg) - 1):
            t1, t2 = fg.loc[i, 'EXPTIME'], fg.loc[i + 1, 'EXPTIME']
            f1 = load_crop(fg.loc[i, 'path']) - mbias - dark_rate_adu * t1
            f2 = load_crop(fg.loc[i + 1, 'path']) - mbias - dark_rate_adu * t2
            m1, m2 = float(np.median(f1[~hot])), float(np.median(f2[~hot]))
            if not (PTC_RANGE[0] < m1 < PTC_RANGE[1] and PTC_RANGE[0] < m2 < PTC_RANGE[1]):
                continue
            r = m1 / m2
            vd = float(mad_std((f1 - r * f2)[~hot]) ** 2)
            x = m1 + r ** 2 * m2                  # signal term [ADU]
            y = vd - (1 + r ** 2) * rn_adu ** 2    # shot-noise variance [ADU^2] = x / gain
            ptc_rows.append({'unit': info['unit'], 'night': info['night'], 'FILTER': filt, 't1': t1, 't2': t2,
                             'm1_adu': m1, 'm2_adu': m2, 'ratio': r, 'var_diff_adu2': vd, 'x': x, 'y': y,
                             'gain_pair': x / y if y > 0 else np.nan})
    if ptc_rows:
        x = np.array([p['x'] for p in ptc_rows]); y = np.array([p['y'] for p in ptc_rows])
        info.update(gain_ptc=float(np.sum(x ** 2) / np.sum(x * y)),
                    gain_ptc_med=float(np.nanmedian([p['gain_pair'] for p in ptc_rows])), n_ptc_pairs=len(ptc_rows))
    # --- SCIENCE: sky level and empty-sky noise
    sky_rows = []
    sel = []
    for _, sg in sci.sort_values('MJD').groupby('FILTER'):
        if MAX_SCI_PER_FILTER is not None and len(sg) > MAX_SCI_PER_FILTER:
            sg = sg.iloc[np.unique(np.linspace(0, len(sg) - 1, MAX_SCI_PER_FILTER).astype(int))]
        sel.append(sg)
    sel = pd.concat(sel) if sel else sci.iloc[:0]
    for _, row in sel.iterrows():
        t = float(row['EXPTIME'])
        img = load_crop(row['path']) - mbias - dark_rate_adu * t
        img[hot] = np.nan
        tmed, tsig = tile_stats(img)
        good = np.isfinite(tmed) & np.isfinite(tsig)
        good &= tsig < 1.5 * np.nanmedian(tsig)            # reject tiles with bright stars / artifacts
        sky_adu, sig_adu = float(np.median(tmed[good])), float(np.median(tsig[good]))
        sky_e = sky_adu * g
        # calibration noise from master bias (nb frames) and scaled master dark (nd frames)
        var_cal = rn_pair ** 2 / nb + ((t / t_d) ** 2 * (rn_pair ** 2 + D_med * t_d) / nd if nd >= 2 else 0.0)
        sig_pred = np.sqrt(max(sky_e, 0) + D_med * t + rn_pair ** 2 + var_cal)
        sky_rows.append({'unit': info['unit'], 'night': info['night'], 'file': Path(row['path']).name,
                         'OBJECT': row['OBJECT'], 'FILTER': row['FILTER'], 'EXPTIME': t, 'MJD': row['MJD'],
                         'AIRMASS': row['AIRMASS'], 'MOONPHAS': row['MOONPHAS'], 'MOONSEP': row['MOONSEP'],
                         'CLUDFRAC': row['CLUDFRAC'], 'sky_adu': sky_adu, 'sky_e_s': sky_e / t,
                         'sig_meas_e': sig_adu * g, 'sig_pred_e': sig_pred, 'noise_ratio': sig_adu * g / sig_pred,
                         't_cross_s': rn_pair ** 2 / (sky_e / t + D_med) if sky_e > 0 else np.nan})
    return info, ptc_rows, sky_rows

def analyze_night_safe(records):
    try:
        return analyze_night(records)
    except Exception as e:
        return {'unit': records[0]['unit'], 'night': records[0]['night'], 'error': repr(e)}, [], []

#%% Main
if __name__ == '__main__':
    OUT_DIR.mkdir(exist_ok=True)
    dirs = night_dirs()
    print(f'{len(dirs)} night folders found')
    with ThreadPoolExecutor(max_workers=N_COPY) as pool:
        local_dirs = list(tqdm(pool.map(sync_heads, dirs), total=len(dirs), desc='Copying .head'))
    with ProcessPoolExecutor(max_workers=N_PROC) as pool:
        idx = pd.DataFrame([r for rows in tqdm(pool.map(index_night, local_dirs), total=len(local_dirs), desc='Indexing')
                            for r in rows])
    for k in NUM_KEYS:
        idx[k] = pd.to_numeric(idx[k], errors='coerce')
    idx.to_csv(OUT_DIR / 'index_2026.csv', index=False)
    typ = idx['IMAGETYP'].astype(str).str.strip().str.upper()
    print('\n[Frame types]\n', typ.value_counts().to_string())
    print('\n[Gain settings per unit (nights)]\n', idx.groupby(['unit', 'GAIN'])['night'].nunique().to_string())
    sci = idx[~typ.isin(CAL_TYPES)]
    print('\n[Science EXPTIME distribution]\n', sci['EXPTIME'].value_counts().sort_index().to_string())
    groups = [g.to_dict('records') for _, g in idx.groupby(['unit', 'night'])]
    with ProcessPoolExecutor(max_workers=N_PROC) as pool:
        results = list(tqdm(pool.map(analyze_night_safe, groups), total=len(groups), desc='Analyzing'))
    det = pd.DataFrame([r[0] for r in results])
    ptc = pd.DataFrame([p for r in results for p in r[1]])
    sky = pd.DataFrame([s for r in results for s in r[2]])
    det.to_csv(OUT_DIR / 'detector_2026.csv', index=False)
    ptc.to_csv(OUT_DIR / 'flat_ptc_2026.csv', index=False)
    sky.to_csv(OUT_DIR / 'sky_2026.csv', index=False)
    cols = [c for c in ['egain', 'gain_ptc', 'rn_pair_e', 'rn_med_e', 'rn_rms_e', 'dark_med_e_s', 'hot_frac'] if c in det]
    print('\n[Detector (median over nights)]\n', det.groupby(['unit', 'gain_set'])[cols].median().to_string())
    if len(sky):
        print('\n[Sky / noise check (median per filter)]\n',
              sky.groupby('FILTER')[['sky_e_s', 'sig_meas_e', 'sig_pred_e', 'noise_ratio', 't_cross_s']].median().to_string())
    if 'error' in det:
        print(f"\n{det['error'].notna().sum()} nights failed (see 'error' column in detector_2026.csv)")
