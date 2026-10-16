#%%
"""
SKYSIG (observed) vs noise-model prediction, gain 2750 / 1x1 binning, all exposure times
  predicted  sigma_pred^2 [e-^2] = SKYVAL*g + DC*t + RN^2      (SKYSIG is NOT used)
  observed   sigma_obs    [e-]   = SKYSIG*g
Outputs (OUT_DIR): skysig_model_frames.csv, skysig_model_by_exptime.csv, skysig_model_by_exptime_filter.csv,
                   skysig_vs_model.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSV = Path('/data3/cwkim/exptimecal/noise_v2_frames.csv')
OUT_DIR = Path('/data3/cwkim/exptimecal')
GAIN_SET, BINNING = 2750, 1
G = 0.256          # [e-/ADU] PTC gain
RN = 1.50          # [e-] read noise (bias, RMS)
DC = 0.0005        # [e-/s/pix] dark current
SKY_BINS = np.logspace(-1, 4, 26)   # [e-] sky per pixel per frame
COLS = ['processed_id', 'processed_filename', 'unit_id', 'filter_id', 'exptime', 'gain', 'egain', 'binning_x',
        'sanity_flag', 'reduction_status', 'skyval', 'skysig', 'airmass', 'moon_phase', 'obstime']

#%% Load and select
df = pd.read_csv(CSV, usecols=COLS)
sel = (df['gain'] == GAIN_SET) & (df['binning_x'] == BINNING) & (np.abs(df['egain'] - G) < 0.01)
sel &= (df['sanity_flag'] == 't') & (df['skyval'] > 0) & (df['skysig'] > 0)
df = df[sel].copy()
print(f'{len(df)} frames (gain {GAIN_SET}, {BINNING}x{BINNING}, sanity, skyval/skysig > 0)')

#%% Prediction vs observation
t = df['exptime']
df['sky_e'] = df['skyval'] * G
df['sky_rate'] = df['sky_e'] / t
df['sig_obs_e'] = df['skysig'] * G
df['sig_pred_e'] = np.sqrt(df['sky_e'] + DC * t + RN ** 2)
df['ratio'] = df['sig_obs_e'] / df['sig_pred_e']
df['excess_var'] = df['sig_obs_e'] ** 2 - df['sig_pred_e'] ** 2     # [e-^2]
df['rn_frac'] = RN ** 2 / df['sig_pred_e'] ** 2                      # RN share of predicted variance
df.to_csv(OUT_DIR / 'skysig_model_frames.csv', index=False)

def summary(g):
    return pd.Series({'n': len(g), 'sky_e_med': g['sky_e'].median(), 'rn_frac_med': g['rn_frac'].median(),
                      'ratio_med': g['ratio'].median(), 'ratio_p16': g['ratio'].quantile(0.16),
                      'ratio_p84': g['ratio'].quantile(0.84), 'excess_var_med': g['excess_var'].median()})

by_t = df.groupby('exptime').apply(summary)
by_t.to_csv(OUT_DIR / 'skysig_model_by_exptime.csv')
by_tf = df.groupby(['exptime', 'filter_id']).apply(summary)
by_tf.to_csv(OUT_DIR / 'skysig_model_by_exptime_filter.csv')
pd.set_option('display.width', 200)
print('\n[sigma_obs / sigma_pred per exposure time]')
print(by_t.round(3).to_string())

#%% Plots
fig, axes = plt.subplots(1, 3, figsize=(17, 5))
ax = axes[0]
hb = ax.hexbin(df['sig_pred_e'], df['sig_obs_e'], xscale='log', yscale='log', bins='log', gridsize=80, mincnt=1)
lim = [df['sig_pred_e'].min() * 0.8, df['sig_obs_e'].max() * 1.2]
ax.plot(lim, lim, 'r-', lw=1)
ax.set_xlabel(r'$\sigma_{pred}$ [e$^-$]  (SKYVAL, RN, DC)')
ax.set_ylabel(r'$\sigma_{obs}$ = SKYSIG$\cdot g$ [e$^-$]')
fig.colorbar(hb, ax=ax, label='N frames')
ax = axes[1]
df['sky_bin'] = pd.cut(df['sky_e'], SKY_BINS)
mid = np.sqrt(SKY_BINS[:-1] * SKY_BINS[1:])
for et, g in df.groupby('exptime'):
    if len(g) < 50:
        continue
    med = g.groupby('sky_bin', observed=False)['ratio'].median().to_numpy()
    cnt = g.groupby('sky_bin', observed=False).size().to_numpy()
    ok = cnt >= 20
    ax.plot(mid[ok], med[ok], 'o-', ms=3, label=f'{et:g} s')
ax.axhline(1, color='k', lw=0.8)
ax.axvline(RN ** 2, color='gray', ls='--', lw=0.8)
ax.text(RN ** 2, 1.0, ' sky = RN$^2$', color='gray', fontsize=8, va='bottom')
ax.set_xscale('log')
ax.set_xlabel(r'sky per pixel per frame [e$^-$]')
ax.set_ylabel(r'median $\sigma_{obs}/\sigma_{pred}$')
ax.legend(fontsize=7, ncol=2)
ax = axes[2]
ok = by_t['n'] >= 50
ax.errorbar(by_t.index[ok], by_t.loc[ok, 'ratio_med'],
            yerr=[by_t.loc[ok, 'ratio_med'] - by_t.loc[ok, 'ratio_p16'], by_t.loc[ok, 'ratio_p84'] - by_t.loc[ok, 'ratio_med']],
            fmt='o', capsize=3)
ax.axhline(1, color='k', lw=0.8)
ax.set_xscale('log')
ax.set_xlabel('exposure time [s]')
ax.set_ylabel(r'$\sigma_{obs}/\sigma_{pred}$ (median, 16-84%)')
fig.suptitle(f'SKYSIG vs noise model (gain {GAIN_SET}, {BINNING}x{BINNING}; g={G}, RN={RN}, DC={DC})')
fig.tight_layout()
fig.savefig(OUT_DIR / 'skysig_vs_model.png', dpi=130)
