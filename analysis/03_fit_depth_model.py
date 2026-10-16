#%%
"""
Per-band linear depth models for 5" (UL5_4) and 10" (UL5_5) aperture photometry, 100 s frames
    UL5 = b0 + b1*moon_phase + b2*moon_sep + b3*airmass + b4*seeing
Data: noise_v2_frames.csv (DB query), GAIN 2750, 1x1 binning, EGAIN ~ 0.256, sanity flag
Validation: nights split into train (80%) / test (20%); iterative 3-sigma clipped least squares on train,
            residual scatter on train and test, predicted vs observed and residual vs feature plots.
Outputs (OUT_DIR): depth_model_ul5_4_<stamp>.json, depth_model_ul5_5_<stamp>.json,
                   depth_model_validation_<stamp>.csv
JSON per filter: intercept, coefficients, features, medians, scatter, n_used, n_total, exptime, aperture, zp_e_s
"""
import json
import re
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import mad_std

CSV = Path('/data3/cwkim/exptimecal/noise_v2_frames.csv')
OUT_DIR = Path('/data3/cwkim/exptimecal/exposurecalculator/exposurecalculator/update_261006')
GAIN_SET, BINNING, EGAIN, EXPTIME = 2750, 1, 0.256, 100.0
FEATURES = ['moon_phase', 'moon_sep', 'airmass', 'seeing']
APERTURES = {'ul5_4': {'zp': 'zp_4', 'diameter': 5.0}, 'ul5_5': {'zp': 'zp_5', 'diameter': 10.0}}
TEST_FRAC, SEED = 0.2, 0
CLIP_SIGMA, CLIP_ITERS, MIN_FRAMES = 3.0, 10, 200
FILT_RE = re.compile(r'_(m\d{3}w?|[ugriz])_7DT\d{2}_')
STAMP = datetime.now().strftime('%Y%m%d%H%M%S')

#%% Load and select
cols = ['processed_filename', 'mjd', 'exptime', 'gain', 'egain', 'binning_x', 'sanity_flag', 'skysig',
        'airmass', 'moon_phase', 'moon_sep', 'seeing', 'zp_4', 'zp_5', 'ul5_4', 'ul5_5']
df = pd.read_csv(CSV, usecols=cols)
df['filter'] = df['processed_filename'].str.extract(FILT_RE, expand=False)
sel = (df['gain'] == GAIN_SET) & (df['binning_x'] == BINNING) & (np.abs(df['egain'] - EGAIN) < 0.01)
sel &= (df['exptime'] == EXPTIME) & (df['sanity_flag'] == 't') & (df['skysig'] > 0)
sel &= (df['airmass'] >= 1) & (df['seeing'] > 0) & df['filter'].notna()
df = df[sel].dropna(subset=FEATURES).copy()
df['night'] = np.floor(df['mjd'] - 16 / 24).astype(int)   # local (Chile) night: boundary at ~local noon (16 UT)
rng = np.random.default_rng(SEED)
nights = df['night'].unique()
test_nights = set(rng.choice(nights, int(len(nights) * TEST_FRAC), replace=False))
df['test'] = df['night'].isin(test_nights)
print(f"{len(df)} frames, {df['filter'].nunique()} filters, {len(nights)} nights ({len(test_nights)} test)")

#%% Fit
def clipped_lstsq(X, y):
    m = np.isfinite(y)
    for _ in range(CLIP_ITERS):
        A = np.column_stack([np.ones(m.sum()), X[m]])
        beta = np.linalg.lstsq(A, y[m], rcond=None)[0]
        res = y - (beta[0] + X @ beta[1:])
        new = np.isfinite(res) & (np.abs(res - np.median(res[m])) < CLIP_SIGMA * mad_std(res[m]))
        if np.array_equal(new, m):
            break
        m = new
    return beta, m

def predict(model, X):
    return model['intercept'] + X @ np.asarray(model['coefficients'])

models, rows = {k: {} for k in APERTURES}, []
for ul5, ap in APERTURES.items():
    for filt, g in df.dropna(subset=[ul5, ap['zp']]).groupby('filter'):
        tr, te = g[~g['test']], g[g['test']]
        if len(tr) < MIN_FRAMES:
            continue
        X, y = tr[FEATURES].to_numpy(float), tr[ul5].to_numpy(float)
        beta, m = clipped_lstsq(X, y)
        zp_e_s = float(np.median(tr[ap['zp']] - 2.5 * np.log10(tr['exptime'] / tr['egain'])))
        res_tr = y - (beta[0] + X @ beta[1:])
        models[ul5][filt] = {'intercept': float(beta[0]), 'coefficients': beta[1:].tolist(), 'features': FEATURES,
                             'medians': np.median(X[m], axis=0).tolist(), 'scatter': float(mad_std(res_tr[m])),
                             'n_used': int(m.sum()), 'n_total': len(tr), 'exptime': EXPTIME, 'filter': filt,
                             'aperture': ap['diameter'], 'zp_e_s': zp_e_s}
        row = {'ul5': ul5, 'filter': filt, 'n_train': len(tr), 'n_test': len(te), 'scatter_train': mad_std(res_tr[m]),
               'bias_train': np.median(res_tr[m]), 'zp_e_s': zp_e_s, **dict(zip(['b0'] + FEATURES, beta))}
        if len(te):
            res_te = te[ul5].to_numpy(float) - predict(models[ul5][filt], te[FEATURES].to_numpy(float))
            row.update(scatter_test=mad_std(res_te), bias_test=np.median(res_te), outlier_frac_test=np.mean(np.abs(res_te) > 1.0))
        rows.append(row)
    out = OUT_DIR / f'depth_model_{ul5}_{STAMP}.json'
    out.write_text(json.dumps(models[ul5], indent=2))
    print(f'saved {out}')
val = pd.DataFrame(rows)
val.to_csv(OUT_DIR / f'depth_model_validation_{STAMP}.csv', index=False)
pd.set_option('display.width', 220)
print(val.round(4).to_string(index=False))

#%% Plot 1: predicted vs observed (test nights)
def band_key(f):
    m = re.match(r'm(\d+)', f)
    return (0, 'griz'.find(f)) if f in 'griz' else (1, int(m.group(1)) if m else 9999)

for ul5 in APERTURES:
    filters = sorted(models[ul5], key=band_key)
    nc = 6
    nr = int(np.ceil(len(filters) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(3 * nc, 3 * nr), squeeze=False)
    for ax, f in zip(axes.flat, filters):
        te = df[(df['filter'] == f) & df['test']].dropna(subset=[ul5])
        if not len(te):
            ax.set_visible(False)
            continue
        p = predict(models[ul5][f], te[FEATURES].to_numpy(float))
        ax.hexbin(p, te[ul5], gridsize=40, bins='log', mincnt=1)
        lo, hi = np.percentile(np.r_[p, te[ul5]], [0.5, 99.5])
        ax.plot([lo, hi], [lo, hi], 'r-', lw=0.8)
        ax.set(xlim=(lo, hi), ylim=(lo, hi))
        ax.set_title(f"{f}  sc={mad_std(te[ul5] - p):.2f}", fontsize=9)
    for ax in axes.flat[len(filters):]:
        ax.set_visible(False)
    fig.supxlabel(f'predicted {ul5.upper()}')
    fig.supylabel(f'observed {ul5.upper()}')
    fig.suptitle(f'{ul5.upper()} model, test nights')
    fig.tight_layout()
    plt.show()

#%% Plot 2: residual (test) vs each feature, all bands
for ul5 in APERTURES:
    te = df[df['test']].dropna(subset=[ul5])
    te = te[te['filter'].isin(models[ul5])]
    res = te[ul5].to_numpy(float) - np.array([predict(models[ul5][f], x[None, :])[0]
                                              for f, x in zip(te['filter'], te[FEATURES].to_numpy(float))])
    fig, axes = plt.subplots(1, len(FEATURES), figsize=(4.5 * len(FEATURES), 3.5), sharey=True)
    for ax, feat in zip(axes, FEATURES):
        ax.hexbin(te[feat], res, gridsize=50, bins='log', mincnt=1)
        bins = np.quantile(te[feat], np.linspace(0, 1, 16))
        mid = 0.5 * (bins[:-1] + bins[1:])
        med = [np.median(res[(te[feat] >= a) & (te[feat] < b)]) if np.any((te[feat] >= a) & (te[feat] < b)) else np.nan
               for a, b in zip(bins[:-1], bins[1:])]
        ax.plot(mid, med, 'r.-')
        ax.axhline(0, color='k', lw=0.8)
        ax.set(xlabel=feat, ylim=(-1.5, 1.5))
    axes[0].set_ylabel(f'{ul5.upper()} obs - pred')
    fig.suptitle(f'{ul5.upper()} residuals vs features (test nights, all bands)')
    fig.tight_layout()
    plt.show()