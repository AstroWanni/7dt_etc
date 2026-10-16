#%%
"""
Independent noise validation on T02524 (up to 20 epochs per filter, GAIN 2750, 100 s singles)
  0) Copy <PROC_ROOT>/2026-*/T02524/<filter>/singles/*_100s.fits (+ _cat.fits, .header) to DEST/<filter>/
  1) Empty-aperture test: random apertures on source-masked, background-subtracted images
     -> measured aperture noise vs header (SKYSIG*sqrt(Npix)) and noise model sqrt(Npix*(sky_e+RN^2))/g
  2) Repeatability: pair differences of consecutive frames (same visit) -> observed scatter vs predicted SNR
"""
import re
import shutil
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import ndimage

PROC_ROOT = Path('/lyman/data2/processed')
BASE = Path('/data3/cwkim/exptimecal/exposurecalculator/exposurecalculator/update_261006')
DEST = BASE / 'T02524'
TILE, YEAR, MAX_EPOCH = 'T02524', '2026', 20
RN_E = 1.58                       # noise-model read noise [e-]
RADII = [2, 4, 6, 8, 9.901]       # aperture radii [pix]; 9.901 = APER_5/2
N_APER, BOX, EDGE = 3000, 64, 100
N_PROC = 5                        # max cores for analysis
MAX_PAIR_DT = 30 / 1440           # consecutive frames within 30 min form a pair [day]
FILTERS = [f'm{w}' for w in range(400, 876, 25)]   # m400, m425, ..., m875 only
RE_FN = re.compile(rf'^{TILE}_(?P<filt>[^_]+)_(?P<unit>7DT\d\d)_(?P<dt>\d{{8}}_\d{{6}})_100s\.fits$')

#%% 0) Copy
def list_remote():
    rows = []
    for fp in PROC_ROOT.glob(f'{YEAR}-*/{TILE}/*/singles/*_100s.fits'):
        m = RE_FN.match(fp.name)
        if m and m['filt'] in FILTERS and fp.with_name(fp.stem + '_cat.fits').exists():
            rows.append({'path': fp, 'filter': m['filt'], 'unit': m['unit'], 'dt': m['dt']})
    return pd.DataFrame(rows).sort_values('dt')

def select_epochs(df):
    """Up to MAX_EPOCH frames per filter, newest nights first, keeping whole nights together."""
    out = []
    for filt, g in df.groupby('filter'):
        g = g.assign(night=g['dt'].str[:8])
        keep = []
        for _, n in sorted(g.groupby('night'), key=lambda x: x[0], reverse=True):
            if len(keep) + len(n) > MAX_EPOCH:
                continue
            keep.append(n)
        if keep:
            out.append(pd.concat(keep))
    return pd.concat(out).sort_values(['filter', 'dt'])

def copy_frames(sel):
    """Copy via a temporary file and rename; skip files already copied completely (same size)."""
    n_skip = n_copy = 0
    for _, r in sel.iterrows():
        d = DEST / r['filter']
        d.mkdir(parents=True, exist_ok=True)
        for src in [r['path'], r['path'].with_name(r['path'].stem + '_cat.fits'), r['path'].with_suffix('.header')]:
            if not src.exists():
                continue
            dst = d / src.name
            if dst.exists() and dst.stat().st_size == src.stat().st_size:
                n_skip += 1
                continue
            tmp = dst.with_name(dst.name + '.part')
            shutil.copy2(src, tmp)
            tmp.rename(dst)
            n_copy += 1
    print(f'copy: {n_copy} copied, {n_skip} skipped (already complete)')

#%% 1) Empty apertures
def background(img, mask):
    """Coarse median background (BOX x BOX blocks of unmasked pixels), bilinear upsampled."""
    ny, nx = img.shape
    a = np.where(mask, np.nan, img)[: ny // BOX * BOX, : nx // BOX * BOX].reshape(ny // BOX, BOX, nx // BOX, BOX)
    blk = np.nanmedian(a, axis=(1, 3))
    blk = np.where(np.isfinite(blk), blk, np.nanmedian(blk))
    bkg = ndimage.zoom(blk, BOX, order=1)
    return np.pad(bkg, ((0, ny - bkg.shape[0]), (0, nx - bkg.shape[1])), mode='edge')

def source_mask(img, cat, hdr, r_cat=15):
    sig = hdr['SKYSIG']
    bkg0 = background(img, np.zeros_like(img, bool))
    mask = ndimage.binary_dilation((img - bkg0) > 2 * sig, iterations=4) | ~np.isfinite(img)
    yy, xx = np.ogrid[-r_cat:r_cat + 1, -r_cat:r_cat + 1]
    disk = xx ** 2 + yy ** 2 <= r_cat ** 2
    pts = np.zeros_like(mask)
    x = np.clip(np.round(np.asarray(cat['X_IMAGE']) - 1).astype(int), 0, img.shape[1] - 1)
    y = np.clip(np.round(np.asarray(cat['Y_IMAGE']) - 1).astype(int), 0, img.shape[0] - 1)
    pts[y, x] = True
    return mask | ndimage.binary_dilation(pts, structure=disk)

def empty_apertures(img, mask, rng):
    ny, nx = img.shape
    rmax = int(np.ceil(max(RADII)))
    yy, xx = np.mgrid[-rmax:rmax + 1, -rmax:rmax + 1]
    rr = np.hypot(xx, yy)
    sums = {r: [] for r in RADII}
    tries = 0
    while len(sums[RADII[0]]) < N_APER and tries < 50 * N_APER:
        tries += 1
        cx, cy = rng.integers(EDGE, nx - EDGE), rng.integers(EDGE, ny - EDGE)
        sl = np.s_[cy - rmax: cy + rmax + 1, cx - rmax: cx + rmax + 1]
        if mask[sl][rr <= max(RADII)].any():
            continue
        cut = img[sl]
        for r in RADII:
            sums[r].append(cut[rr <= r].sum())
    npix = {r: int((rr <= r).sum()) for r in RADII}
    return {r: np.asarray(v) for r, v in sums.items()}, npix

def frame_empty(fp):
    rng = np.random.default_rng(zlib.crc32(fp.name.encode()))   # reproducible per frame
    with fits.open(fp) as h:
        img, hdr = h[0].data.astype(np.float32), h[0].header
    cat = Table.read(fp.with_name(fp.stem + '_cat.fits'), hdu=2)
    mask = source_mask(img, cat, hdr)
    img = img - background(img, mask)
    sums, npix = empty_apertures(img, mask, rng)
    g, sky_e = hdr['EGAIN'], hdr['SKYVAL'] * hdr['EGAIN']
    rows = []
    for r in RADII:
        meas = sigma_clipped_stats(sums[r], sigma=4)[2]
        rows.append({'file': fp.name, 'filter': hdr['FILTER'], 'unit': hdr['TELESCOP'], 'radius': r, 'npix': npix[r],
                     'n_aper': len(sums[r]), 'sig_meas': meas, 'sig_hdr': hdr['SKYSIG'] * np.sqrt(npix[r]),
                     'sig_model': np.sqrt(npix[r] * (sky_e + RN_E ** 2)) / g, 'mask_frac': mask.mean()})
    return rows

#%% 2) Repeatability (pair differences)
def load_cat(fp):
    with fits.open(fp) as h:
        hdr = h[0].header
    cat = Table.read(fp.with_name(fp.stem + '_cat.fits'), hdu=2)
    g, npix = hdr['EGAIN'], np.pi * (hdr['APER_5'] / 2) ** 2
    f = np.asarray(cat['FLUX_APER_5'], float)
    ok = (f > 0) & (np.asarray(cat['FLAGS']) == 0) & (np.asarray(cat['FLUX_MAX']) < 0.5 * hdr['SATURATE'])
    s = f[ok] * g
    snr_pred = s / np.sqrt(s + npix * (hdr['SKYVAL'] * g + RN_E ** 2))
    return pd.DataFrame({'ra': np.asarray(cat['ALPHA_J2000'])[ok], 'dec': np.asarray(cat['DELTA_J2000'])[ok],
                         'mag': np.asarray(cat[f"MAG_APER_5_{hdr['FILTER']}"])[ok], 'flux': f[ok],
                         'snr_pred': snr_pred, 'snr_cat': f[ok] / np.asarray(cat['FLUXERR_APER_5'], float)[ok]}), hdr['MJD']

def pair_diffs(fps, match_arcsec=1.0):
    cats = [load_cat(fp) for fp in fps]
    rows = []
    for (a, ta), (b, tb), fa, fb in zip(cats[:-1], cats[1:], fps[:-1], fps[1:]):
        if tb - ta > MAX_PAIR_DT:
            continue
        idx, d2d, _ = SkyCoord(a['ra'].values * u.deg, a['dec'].values * u.deg).match_to_catalog_sky(
            SkyCoord(b['ra'].values * u.deg, b['dec'].values * u.deg))
        m = d2d.arcsec < match_arcsec
        A, B = a[m].reset_index(drop=True), b.iloc[idx[m]].reset_index(drop=True)
        dm = -2.5 * np.log10(A['flux'] / B['flux'])
        dm -= np.median(dm[(A['snr_pred'] > 50) & (A['snr_pred'] < 500)])          # per-pair ZP offset
        sig_pred = 1.0857 * np.sqrt(1 / A['snr_pred'] ** 2 + 1 / B['snr_pred'] ** 2)
        sig_cat = 1.0857 * np.sqrt(1 / A['snr_cat'] ** 2 + 1 / B['snr_cat'] ** 2)
        rows.append(pd.DataFrame({'pair': f'{fa.stem}|{fb.stem}', 'mag': A['mag'], 'dm': dm, 'snr_pred': np.sqrt(A['snr_pred'] * B['snr_pred']),
                                  'z_pred': dm / sig_pred, 'z_cat': dm / sig_cat}))
    return pd.concat(rows) if rows else pd.DataFrame()

def bin_repeat(df, edges=np.arange(10, 22.01, 0.5)):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        d = df[(df['mag'] >= lo) & (df['mag'] < hi)]
        if len(d) < 30:
            continue
        sig_single = mad_std(d['dm']) / np.sqrt(2)
        out.append({'mag': 0.5 * (lo + hi), 'n': len(d), 'snr_obs': 1.0857 / sig_single, 'snr_pred': d['snr_pred'].median(),
                    'z_pred_std': mad_std(d['z_pred']), 'z_cat_std': mad_std(d['z_cat'])})
    return pd.DataFrame(out)

#%% Main
if __name__ == '__main__':
    sel = select_epochs(list_remote())
    print(sel.groupby(['filter', 'unit']).size().to_string())
    copy_frames(sel)
    emp, rep = [], []
    for filt in [f for f in FILTERS if f in set(sel['filter'])]:
        fps = sorted(p for p in (DEST / filt).glob(f'{TILE}_*_100s.fits'))
        fps = [p for p in fps if fits.getval(p, 'GAIN') == 2750]
        with ProcessPoolExecutor(max_workers=N_PROC) as pool:
            for fp, rows in zip(fps, pool.map(frame_empty, fps)):
                emp += rows
                print(f'{filt} {fp.name}: empty apertures done')
        pdf = pair_diffs(fps)
        if len(pdf):
            rep.append(bin_repeat(pdf).assign(filter=filt))
            pdf.assign(filter=filt).to_csv(DEST / f'pairdiff_{filt}.csv', index=False)
    emp = pd.DataFrame(emp); emp['ratio_hdr'] = emp['sig_meas'] / emp['sig_hdr']; emp['ratio_model'] = emp['sig_meas'] / emp['sig_model']
    emp.to_csv(DEST / 'empty_aperture.csv', index=False)
    rep = pd.concat(rep) if rep else pd.DataFrame()
    rep.to_csv(DEST / 'repeatability_binned.csv', index=False)
    print('\n[Empty aperture] median sig_meas / prediction')
    print(emp.groupby(['filter', 'radius'])[['ratio_hdr', 'ratio_model']].median().round(3).to_string())
    print('\n[Repeatability] observed vs predicted SNR, z = dm/sigma_pred (std should be 1)')
    print(rep.round(3).to_string(index=False))