#%%
"""
WRITER: Hyeonho Choi (hhchoi1022@gmail.com)
DATE: 2026-09-14
DESCRIPTION: Standalone version of bridge.utils.exposurecalculator.ExposureCalculator.
             Empirical 7DT exposure time calculator (ETC) based on depth (5-sigma limiting
             magnitude) statistics of archived 7DT frames.

             Model (per filter, fitted at a reference exposure time, 100 s):
                 UL5_ref = b0 + b1*moon_phase + b2*moon_sep + b3*hour_since_sunset + b4*seeing
                 t_exp   = t_ref * 10**((mag - UL5_ref) / 1.25) * (snr / 5)**2
                 SNR(t)  = 5 * 10**(-0.4 * (mag - UL5_ref - 1.25*log10(t / t_ref)))

             Dependencies: numpy, astropy, astroplan, matplotlib (plots), tqdm.
             scikit-learn is only needed for refitting the model (fit_depth_model), and the
             BRIDGE package (bridge.connector.GWPortalConnector) only when update_depth_model
             is asked to query the 7DT GW portal itself. Neither is imported at module level.

             Usage:
                 calc = ExposureCalculator()
                 calc.load_model()                      # latest depth_model_*.json next to this file
                 calc.calculate_exptime(filter='all', magnitude=19.5, snr=10, ra=30.97, dec=4.25, obsdate=Time.now())
                 calc.calculate_snr(filter=['g', 'r'], magnitude=19.5, exptime=300, seeing=1.5)
"""
#%%
import glob
import json
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import astropy.units as u
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from astroplan import Observer
from astropy.coordinates import EarthLocation, SkyCoord, get_body
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm

#%%
# 7DT site (same values as bridge/configuration/observer.config)
OBSERVER_NAME = '7DT'
OBSERVER_LATITUDE = -30.4704
OBSERVER_LONGITUDE = -70.7804
OBSERVER_ELEVATION = 1580
OBSERVER_TIMEZONE = 'America/Santiago'


def get_observer() -> Observer:
    location = EarthLocation.from_geodetic(lat=OBSERVER_LATITUDE, lon=OBSERVER_LONGITUDE, height=OBSERVER_ELEVATION)
    return Observer(location=location, name=OBSERVER_NAME, timezone=OBSERVER_TIMEZONE)


def calculate_moonsep_worker(ra, dec, mjd):
    """Moon-target separation [deg] at the given MJD (geocentric moon position)."""
    target = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.deg, u.deg))
    moon = get_body('moon', Time(mjd, format='mjd'))
    return float(np.round(moon.separation(target).deg, 2))


def calculate_moonphase_worker(mjd):
    """Fraction of the Moon illuminated (0 = new, 1 = full) at the given MJD."""
    return float(get_observer().moon_illumination(Time(mjd, format='mjd')))


def calculate_hour_since_sunset_worker(mjd):
    """
    Hours elapsed since astronomical sunset (sun altitude -18 deg) at the given MJD.

    The sunset is located by first finding the *next* sunrise (sun altitude +10 deg) and then
    the sunset preceding it, so the result stays on the same night through morning twilight.
    Do not change this definition: the depth model was fitted with hour_since_sunset computed
    exactly this way.
    """
    observer = get_observer()
    sunrise = observer.sun_rise_time(Time(mjd, format='mjd'), which='next', horizon=10 * u.deg)
    sunset = observer.sun_set_time(sunrise, which='previous', horizon=-18 * u.deg)
    hour_since_sunset = (mjd - sunset.mjd) * 24
    if hour_since_sunset < 0:
        hour_since_sunset += 24
    return float(hour_since_sunset)


#%%
class ExposureCalculator:
    """
    Empirical 7DT exposure time calculator.

    Observing conditions used by calculate_exptime / calculate_snr / plot_*:
      - If ``ra`` and ``dec`` are given, ``obsdate`` (default: now) is used to compute
        moon_phase, moon_separation and hour_since_sunset, overriding any given values.
      - If only ``obsdate`` is given, moon_phase and hour_since_sunset are computed;
        moon_separation cannot be computed without a target position.
      - ``seeing`` is never computed; it is always taken from the argument.
      - Any condition still ``None`` after that falls back to ``DEFAULT_CONDITIONS``.
    """

    # Typical observing conditions used when a condition is not given:
    # [moon_phase, moon_separation (deg), hour_since_sunset (h), seeing (arcsec)]
    DEFAULT_CONDITIONS = (0.5, 90, 5, 2.0)

    BROAD_BAND_COLORS = {
        "g": "green",
        "r": "red",
        "i": "gold",
        "z": "k"
    }

    BROAD_ORDER = ["g", "r", "i", "z"]

    def __init__(self):
        self.db_tbl = Table()
        self.model = None

    # ------------------------------------------------------------
    # Filter helpers
    # ------------------------------------------------------------
    @staticmethod
    def _filter_wavelength(filt):
        """
        Extract wavelength from filter name.
        Examples:
          m400 -> 400
          m875 -> 875
          g/r/i -> None
        """
        m = re.match(r"m(\d+)", filt)
        return int(m.group(1)) if m else None

    def get_filter_colors(self, filters):
        """
        Return dict: {filter_name: color}
        """
        filters = list(filters)

        colors = {}

        # Broad bands
        for f in filters:
            if f in self.BROAD_BAND_COLORS:
                colors[f] = self.BROAD_BAND_COLORS[f]

        # Medium bands
        mfilters = [f for f in filters if f.startswith("m")]

        if len(mfilters) > 0:
            waves = np.array([self._filter_wavelength(f) for f in mfilters])

            norm = mcolors.Normalize(vmin=waves.min(), vmax=waves.max())
            cmap = plt.get_cmap("jet")

            for f, w in zip(mfilters, waves):
                colors[f] = cmap(norm(w))

        return colors

    def sort_filters_broad_then_medium(self, filters):
        broad_rank = {f: i for i, f in enumerate(self.BROAD_ORDER)}

        def key(f):
            f = str(f)
            # 0) broad
            if f in broad_rank:
                return (0, broad_rank[f], -1, f)

            # 1) medium
            w = self._filter_wavelength(f)
            if w is not None:
                return (1, 0, w, f)

            # 2) others
            return (2, 0, 10**9, f)

        return sorted(list(filters), key=key)

    def _filter_list(self, filter) -> list:
        """'all' -> every filter in the model; str -> [str]; iterable -> list of str."""
        if isinstance(filter, str):
            filter = list(self.model.keys()) if filter == 'all' else [filter]
        return [str(f) for f in np.atleast_1d(filter)]

    # ------------------------------------------------------------
    # Model fitting
    # ------------------------------------------------------------
    def update_depth_model(self,
                           depth_table=None,
                           obs_start_date: str = '2025-01-01',
                           obs_end_date: str = '2026-01-01',
                           model_save_path: str = None,
                           depth_table_save_path: str = None,
                           n_processes: int = 32
                           ):
        """
        Fit the depth model for every filter in ``depth_table`` and save it as JSON.

        depth_table : astropy.table.Table or path to a CSV with columns
            filter, exptime, ul5, seeing, ra_center, dec_center, mjd
            (moon_sep, moon_phase, hour_since_sunset are computed if missing).
            None -> query the 7DT GW portal for [obs_start_date, obs_end_date] through
            bridge.connector.GWPortalConnector (requires the BRIDGE package and its API key).
        """
        date_end = Time.now()
        current_dir = Path(__file__).resolve().parent
        if model_save_path is None:
            model_save_path = f'{current_dir}/depth_model_{date_end.datetime.strftime("%Y%m%d%H%M%S")}.json'
        if depth_table_save_path is None:
            depth_table_save_path = f'{current_dir}/depth_data_{date_end.datetime.strftime("%Y%m%d%H%M%S")}.csv'

        if depth_table is None:
            from bridge.connector import GWPortalConnector  # only needed for querying the portal
            self.db_tbl = GWPortalConnector().query(obs_start_date=obs_start_date, obs_end_date=obs_end_date)
        elif isinstance(depth_table, (str, Path)):
            self.db_tbl = Table.read(depth_table)
        else:
            self.db_tbl = depth_table

        if any(col not in self.db_tbl.colnames for col in ('moon_sep', 'moon_phase', 'hour_since_sunset')):
            self._add_condition_columns(n_processes=n_processes)

        filter_names = sorted(list(set(self.db_tbl['filter'])))
        model_dict = {}
        for filter_name in filter_names:
            fit_result = self.fit_depth_model(self.db_tbl, filter_name)
            model_dict[str(filter_name)] = fit_result

        with open(model_save_path, 'w') as f:
            json.dump(model_dict, f, indent=4)
        self.model = model_dict

        for col in ('poly', 'poly_galactic'):
            if col in self.db_tbl.colnames:
                self.db_tbl.remove_column(col)
        self.db_tbl.write(depth_table_save_path, format='csv', overwrite=True)

    def _add_condition_columns(self, n_processes: int = 32):
        """Compute moon_sep, moon_phase, hour_since_sunset for every row of self.db_tbl."""
        with ProcessPoolExecutor(max_workers=n_processes) as pool:
            all_moonsep = list(tqdm(pool.map(calculate_moonsep_worker, self.db_tbl['ra_center'], self.db_tbl['dec_center'], self.db_tbl['mjd']), desc='Calculating moon separation...', total=len(self.db_tbl)))
        with ProcessPoolExecutor(max_workers=n_processes) as pool:
            all_moonphase = list(tqdm(pool.map(calculate_moonphase_worker, self.db_tbl['mjd']), desc='Calculating moon phase...', total=len(self.db_tbl)))
        with ProcessPoolExecutor(max_workers=n_processes) as pool:
            all_hour_since_sunset = list(tqdm(pool.map(calculate_hour_since_sunset_worker, self.db_tbl['mjd']), desc='Calculating hour since sunset...', total=len(self.db_tbl)))

        self.db_tbl['moon_sep'] = all_moonsep
        self.db_tbl['moon_phase'] = all_moonphase
        self.db_tbl['hour_since_sunset'] = all_hour_since_sunset

    def fit_depth_model(
        self,
        tbl,
        filter_name,
        exptime=100,
        sigma=3.0,
        maxiters=3,
        verbose=True,
        visualize=True
        ):
        """
        Fit limiting magnitude (UL5) as a linear function of observing conditions
        at a fixed exposure time, with sigma clipping on residuals.

        Model:
            UL5 = b0
                + b1 * moon_phase
                + b2 * moon_sep
                + b3 * hour_since_sunset
                + b4 * seeing
        """
        from sklearn.linear_model import LinearRegression  # only needed for fitting

        def _to_float_array(col):
            """Safely convert Astropy Column / MaskedColumn to float ndarray."""
            if hasattr(col, 'filled'):          # MaskedColumn
                col = col.filled(np.nan)
            return np.asarray(col, dtype=float)

        # ------------------------------------------------------------
        # 1. Select data
        # ------------------------------------------------------------
        m = (tbl['filter'] == filter_name) & (tbl['exptime'] == exptime)
        tbl_group = tbl[m]

        if len(tbl_group) < 10:
            raise RuntimeError("Not enough data points for fitting")

        moon_phase = _to_float_array(tbl_group['moon_phase'])
        moon_sep   = _to_float_array(tbl_group['moon_sep'])
        hour       = _to_float_array(tbl_group['hour_since_sunset'])
        seeing     = _to_float_array(tbl_group['seeing'])
        ul5        = _to_float_array(tbl_group['ul5'])

        finite = (
            np.isfinite(moon_phase) &
            np.isfinite(moon_sep) &
            np.isfinite(hour) &
            np.isfinite(seeing) &
            np.isfinite(ul5)
        )

        tbl_group = tbl_group[finite]

        X = np.column_stack([
            tbl_group['moon_phase'],
            tbl_group['moon_sep'],
            tbl_group['hour_since_sunset'],
            tbl_group['seeing'],
        ]).astype(float)

        y = np.asarray(tbl_group['ul5'], dtype=float)

        # ------------------------------------------------------------
        # 2. Initial fit
        # ------------------------------------------------------------
        model_init = LinearRegression()
        model_init.fit(X, y)
        y_pred_init = model_init.predict(X)

        # ------------------------------------------------------------
        # 3. Sigma clipping on residuals
        # ------------------------------------------------------------
        resid = y - y_pred_init
        finite = np.isfinite(resid)

        clipped = sigma_clip(
            resid[finite],
            sigma=sigma,
            maxiters=maxiters,
            cenfunc='median',
            stdfunc='mad_std',
        )

        mask_good = np.zeros_like(y, dtype=bool)
        mask_good[finite] = ~clipped.mask

        if verbose:
            print(
                f"[{filter_name}, {exptime}s] "
                f"Rejected {np.sum(~mask_good)}/{len(y)} points "
                f"({100*np.sum(~mask_good)/len(y):.1f}%)"
            )

        # ------------------------------------------------------------
        # 4. Final fit
        # ------------------------------------------------------------
        model = LinearRegression()
        model.fit(X[mask_good], y[mask_good])

        scatter = np.std(y[mask_good] - model.predict(X[mask_good]))
        median = [np.nanmedian(moon_phase), np.nanmedian(moon_sep), np.nanmedian(hour), np.nanmedian(seeing)]

        result = {
            "intercept": model.intercept_,
            "coefficients": list(model.coef_),
            "medians": median,
            "scatter": scatter,
            "n_used": int(np.sum(mask_good)),
            "n_total": int(len(y)),
            "exptime": exptime,
            "filter": str(filter_name)}

        # ------------------------------------------------------------
        # 5. Visualization (PARTIAL DEPENDENCE — 2x2)
        # ------------------------------------------------------------
        if visualize:
            features = [
                ("moon_phase", tbl_group['moon_phase']),
                ("moon_sep", tbl_group['moon_sep']),
                ("hour_since_sunset", tbl_group['hour_since_sunset']),
                ("seeing", tbl_group['seeing']),
            ]

            medians = {
                "moon_phase": np.median(tbl_group['moon_phase']),
                "moon_sep": np.median(tbl_group['moon_sep']),
                "hour_since_sunset": np.median(tbl_group['hour_since_sunset']),
                "seeing": np.median(tbl_group['seeing']),
            }

            fig, axes = plt.subplots(2, 2, figsize=(10, 8))
            axes = axes.flatten()

            for i, (name, values) in enumerate(features):
                ax = axes[i]

                x = np.asarray(values, dtype=float)
                x_grid = np.linspace(np.nanmin(x), np.nanmax(x), 200)

                # Partial-dependence design matrix
                X_grid = np.zeros((len(x_grid), X.shape[1]))
                X_grid[:, 0] = medians["moon_phase"]
                X_grid[:, 1] = medians["moon_sep"]
                X_grid[:, 2] = medians["hour_since_sunset"]
                X_grid[:, 3] = medians["seeing"]
                X_grid[:, i] = x_grid

                y_grid = model.predict(X_grid)

                # Used points
                ax.scatter(
                    x[mask_good],
                    y[mask_good],
                    s=20,
                    facecolor = 'none',
                    edgecolor = 'black',
                    alpha=0.6,
                    label="Used in fit",
                )

                # Rejected points
                ax.scatter(
                    x[~mask_good],
                    y[~mask_good],
                    s=30,
                    facecolor = 'none',
                    edgecolor = 'red',
                    alpha=0.2,
                    label="Sigma-clipped",
                )

                # Linear model
                ax.plot(
                    x_grid,
                    y_grid,
                    "r--",
                    lw=3,
                    label="Linear model",
                )

                ax.set_xlabel(name.replace("_", " ").title(), fontsize = 16)
                ax.set_ylabel("UL5", fontsize =16)
                ax.grid(alpha=0.3)
                ax.set_title(name.replace("_", " ").title(), fontsize = 16)

            # One legend for all panels
            handles, labels = axes[0].get_legend_handles_labels()

            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=3,
                bbox_to_anchor=(0.5, 0.92),
                frameon=True,
                fontsize = 16
            )

            fig.suptitle(f"{filter_name} band", fontsize=25, y=0.96)

            fig.tight_layout(rect=[0, 0, 1, 0.88])
            plt.show()

            # ------------------------------------------------------------
            # 6. Residual diagnostic (UL5 - model prediction, all components)
            # ------------------------------------------------------------
            y_pred_all = model.predict(X)
            residual_all = y - y_pred_all

            resid_used = residual_all[mask_good]
            resid_rej  = residual_all[~mask_good]

            mean_used = np.mean(resid_used)
            std_used  = np.std(resid_used)

            fig_r, axes_r = plt.subplots(1, 2, figsize=(12, 5))

            # --- Left: residual vs predicted UL5 ---
            ax0 = axes_r[0]
            ax0.scatter(
                y_pred_all[mask_good],
                resid_used,
                s=20,
                facecolor='none',
                edgecolor='black',
                alpha=0.6,
                label="Used in fit",
            )
            ax0.scatter(
                y_pred_all[~mask_good],
                resid_rej,
                s=30,
                facecolor='none',
                edgecolor='red',
                alpha=0.3,
                label="Sigma-clipped",
            )
            ax0.axhline(0.0, color='red', ls='--', lw=2)
            ax0.axhline(mean_used + std_used, color='gray', ls=':', lw=1)
            ax0.axhline(mean_used - std_used, color='gray', ls=':', lw=1)
            ax0.set_xlabel("Predicted UL5", fontsize=16)
            ax0.set_ylabel("Residual (UL5$_{obs}$ - UL5$_{model}$)", fontsize=16)
            ax0.set_title("Residual vs Predicted UL5", fontsize=16)
            ax0.grid(alpha=0.3)
            ax0.legend(fontsize=12, loc='best')

            # --- Right: residual histogram ---
            ax1 = axes_r[1]
            bins = np.linspace(
                np.nanmin(residual_all[np.isfinite(residual_all)]),
                np.nanmax(residual_all[np.isfinite(residual_all)]),
                40,
            )
            ax1.hist(
                resid_used, bins=bins,
                color='black', alpha=0.5, label="Used in fit",
            )
            ax1.hist(
                resid_rej, bins=bins,
                color='red', alpha=0.4, label="Sigma-clipped",
            )
            ax1.axvline(0.0, color='red', ls='--', lw=2)
            ax1.axvline(mean_used, color='blue', ls='-', lw=1.5,
                        label=f"mean = {mean_used:+.3f}")
            ax1.axvline(mean_used + std_used, color='gray', ls=':', lw=1)
            ax1.axvline(mean_used - std_used, color='gray', ls=':', lw=1,
                        label=f"std = {std_used:.3f}")
            ax1.set_xlabel("Residual (UL5$_{obs}$ - UL5$_{model}$)", fontsize=16)
            ax1.set_ylabel("Count", fontsize=16)
            ax1.set_title("Residual distribution", fontsize=16)
            ax1.grid(alpha=0.3)
            ax1.legend(fontsize=11, loc='best')

            fig_r.suptitle(
                f"{filter_name} band — full-model residual "
                f"(N$_{{used}}$ = {int(np.sum(mask_good))}, scatter = {scatter:.3f})",
                fontsize=18, y=0.99,
            )
            fig_r.tight_layout(rect=[0, 0, 1, 0.93])
            plt.show()

        return result

    # ------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------
    def load_model(self, model_path: str = None):
        """Load a depth model JSON. Default: the latest depth_model_*.json next to this file."""
        if model_path is None:
            model_list = sorted(glob.glob(f'{Path(__file__).parent}/depth_model_*.json'))
            if len(model_list) == 0:
                raise FileNotFoundError(f"No depth_model_*.json found in {Path(__file__).parent}")
            model_path = model_list[-1]

        with open(model_path, 'r') as f:
            self.model = json.load(f)
        print(f"Loading model from {model_path}")

    def load_depth_table(self, depth_table_path: str = None):
        """Load a depth table CSV. Default: the latest depth_data_*.csv next to this file."""
        if depth_table_path is None:
            table_list = sorted(glob.glob(f'{Path(__file__).parent}/depth_data_*.csv'))
            if len(table_list) == 0:
                raise FileNotFoundError(f"No depth_data_*.csv found in {Path(__file__).parent}")
            depth_table_path = table_list[-1]
        self.db_tbl = Table.read(depth_table_path)
        print(f"Loading depth table from {depth_table_path}")

    # ------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------
    def _check_model(self):
        if self.model is None:
            raise ValueError("Model is not loaded. Please load the model first with self.load_model(model_path)")

    def _resolve_conditions(self,
                            ra, dec, obsdate,
                            moon_phase, moon_separation, hour_since_sunset, seeing,
                            verbose: bool) -> dict:
        """
        Determine the observing conditions to feed the depth model (see class docstring).
        Returns dict with moon_phase, moon_separation, hour_since_sunset, seeing (floats)
        and obsdate (astropy Time or None).
        """
        has_target = (ra is not None) and (dec is not None)
        if has_target and obsdate is None:
            obsdate = Time.now()
        if obsdate is not None:
            obsdate = Time(obsdate)
            moon_phase = calculate_moonphase_worker(obsdate.mjd)
            hour_since_sunset = calculate_hour_since_sunset_worker(obsdate.mjd)
            if has_target:
                moon_separation = calculate_moonsep_worker(ra, dec, obsdate.mjd)
            if verbose:
                print(f"Calculated from the target information (RA = {ra}, Dec = {dec}, Obsdate = {obsdate.iso})")
                print(f"Moon phase = {moon_phase} \nMoon separation = {moon_separation} \nHour since sunset = {hour_since_sunset}")

        defaults = self.DEFAULT_CONDITIONS
        return {
            'moon_phase': float(defaults[0] if moon_phase is None else moon_phase),
            'moon_separation': float(defaults[1] if moon_separation is None else moon_separation),
            'hour_since_sunset': float(defaults[2] if hour_since_sunset is None else hour_since_sunset),
            'seeing': float(defaults[3] if seeing is None else seeing),
            'obsdate': obsdate,
        }

    def _broadcast(self, filter, *others) -> list:
        """Broadcast filter and the other inputs to a common length (length-1 inputs are repeated)."""
        arrays = [np.atleast_1d(a) for a in (self._filter_list(filter), *others)]
        n = max(len(a) for a in arrays)
        for i, a in enumerate(arrays):
            if len(a) == 1:
                arrays[i] = np.repeat(a, n)
            elif len(a) != n:
                raise ValueError("filter, magnitude, and snr/exptime must have the same length (or length 1)")
        return arrays

    def _predict_ul5_ref(self, filt: str, conditions: dict):
        """UL5 (5-sigma limiting magnitude) at the model's reference exposure time. Returns (ul5_ref, t_ref)."""
        if filt not in self.model:
            raise KeyError(f"Filter '{filt}' is not in the depth model. Available: {list(self.model.keys())}")
        model = self.model[filt]
        beta0 = float(model["intercept"])
        beta = np.asarray(model["coefficients"], dtype=float)  # [b_phase, b_sep, b_hour, b_seeing]
        t_ref = float(model["exptime"])  # e.g. 100

        ul5_ref = (
            beta0
            + beta[0] * conditions['moon_phase']
            + beta[1] * conditions['moon_separation']
            + beta[2] * conditions['hour_since_sunset']
            + beta[3] * conditions['seeing']
        )
        return float(ul5_ref), t_ref

    def calculate_exptime(self,
                          filter: str or list[str],
                          magnitude: float or list[float],
                          snr: float or list[float] = 5.0,
                          ra: float = None,
                          dec: float = None,
                          obsdate: Time = None,
                          moon_phase: float = None,
                          moon_separation: float = None,
                          hour_since_sunset: float = None,
                          seeing: float = None,
                          verbose: bool = True) -> dict:
        """
        Exposure time [s] needed to detect ``magnitude`` at ``snr`` in each filter.

        filter    : filter name, list of names, or 'all'
        magnitude : AB magnitude (scalar or list matching filter)
        snr       : target signal-to-noise ratio (scalar or list matching filter)
        ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing : see class docstring

        Returns {filter: {t_exp, magnitude, snr, ul5_ref, moon_phase, moon_separation, hour_since_sunset, seeing}}
        """
        self._check_model()
        conditions = self._resolve_conditions(ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing, verbose)
        filter, magnitude, snr = self._broadcast(filter, magnitude, snr)

        result = dict()
        for filt, mag, snr_i in zip(filter, magnitude, snr):
            filt = str(filt)
            mag = float(mag)
            snr_i = float(snr_i)
            ul5_ref, t_ref = self._predict_ul5_ref(filt, conditions)

            # UL5 is defined at SNR=5; scale by (snr/5)^2 for the requested SNR.
            t_exp = t_ref * 10 ** ((mag - ul5_ref) / 1.25) * (snr_i / 5.0) ** 2
            result[filt] = {
                't_exp': float(t_exp),
                'magnitude': mag,
                'snr': snr_i,
                'ul5_ref': ul5_ref,
                'moon_phase': conditions['moon_phase'],
                'moon_separation': conditions['moon_separation'],
                'hour_since_sunset': conditions['hour_since_sunset'],
                'seeing': conditions['seeing'],
            }
        return result

    def calculate_snr(self,
                      filter: str or list[str],
                      magnitude: float or list[float],
                      exptime: float or list[float],
                      ra: float = None,
                      dec: float = None,
                      obsdate: Time = None,
                      moon_phase: float = None,
                      moon_separation: float = None,
                      hour_since_sunset: float = None,
                      seeing: float = None,
                      verbose: bool = True) -> dict:
        """
        Expected SNR of a ``magnitude`` source in an ``exptime`` [s] exposure in each filter.

        filter    : filter name, list of names, or 'all'
        magnitude : AB magnitude (scalar or list matching filter)
        exptime   : exposure time in seconds (scalar or list matching filter)
        ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing : see class docstring

        Returns {filter: {snr, magnitude, exptime, ul5_exptime, filter, moon_phase, moon_separation, hour_since_sunset, seeing}}
        """
        self._check_model()
        conditions = self._resolve_conditions(ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing, verbose)
        filter, magnitude, exptime = self._broadcast(filter, magnitude, exptime)

        result = {}
        for filt, mag, t in zip(filter, magnitude, exptime):
            filt = str(filt)
            mag = float(mag)
            t = float(t)
            ul5_ref, t_ref = self._predict_ul5_ref(filt, conditions)

            # Scale UL5 to the requested exposure time, then convert to SNR
            ul5_t = ul5_ref + 1.25 * np.log10(t / t_ref)
            snr_val = 5.0 * 10 ** (-0.4 * (mag - ul5_t))

            result[filt] = {
                "snr": float(snr_val),
                "magnitude": mag,
                "exptime": t,
                "ul5_exptime": float(ul5_t),
                "filter": filt,
                "moon_phase": conditions['moon_phase'],
                "moon_separation": conditions['moon_separation'],
                "hour_since_sunset": conditions['hour_since_sunset'],
                "seeing": conditions['seeing'],
            }
        return result

    # ------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------
    def _conditions_text(self, ra, dec, conditions: dict) -> str:
        lines = []
        if ra is not None:
            lines.append(f"RA = {ra:.2f}")
        if dec is not None:
            lines.append(f"Dec = {dec:.2f}")
        if conditions['obsdate'] is not None:
            lines.append(f"Obsdate = {conditions['obsdate'].iso}")
        lines += [
            f"Moon phase = {conditions['moon_phase']:.2f}",
            f"Moon separation = {conditions['moon_separation']:.2f}",
            f"Hour since sunset = {conditions['hour_since_sunset']:.2f}",
            f"Seeing = {conditions['seeing']:.2f}",
        ]
        return "\n".join(lines)

    def plot_snr_vs_exptime(self,
                            filter: str or list[str],
                            magnitude=20.0,
                            exptime_range=(1, 1000),
                            exptime_ref = 300,
                            n_points=100,
                            ra=None,
                            dec=None,
                            obsdate=None,
                            moon_phase=None,
                            moon_separation=None,
                            hour_since_sunset=None,
                            seeing=None,
                            save_path: str = None):
        self._check_model()
        conditions = self._resolve_conditions(ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing, verbose=False)
        filters = self._filter_list(filter)

        exptimes = np.logspace(
            np.log10(exptime_range[0]),
            np.log10(exptime_range[1]),
            n_points,
        )
        exptimes = np.append(exptimes, exptime_ref)

        tbl = Table()
        tbl["exptime"] = exptimes
        for filt in filters:
            tbl[filt] = np.nan
        for i, t in enumerate(exptimes):
            result = self.calculate_snr(
                filter=filters,
                magnitude=magnitude,
                exptime=t,
                moon_phase=conditions['moon_phase'],
                moon_separation=conditions['moon_separation'],
                hour_since_sunset=conditions['hour_since_sunset'],
                seeing=conditions['seeing'],
                verbose=False,
            )
            for filt in filters:
                tbl[filt][i] = result[filt]["snr"]

        plt.figure(figsize=(8, 12))
        xmax = exptimes[-2]
        plt.xlim(exptimes[0], xmax * 2)
        colors = self.get_filter_colors(filters)
        for filt in filters:
            plt.plot(exptimes, tbl[filt], label=f'{filt} (SNR(t_exp={exptime_ref}s) = {tbl[filt][-1]:.1f})', c = colors[filt])
            plt.text(exptimes[-2] * 1.05, tbl[filt][-2], filt, color=colors[filt], fontsize=12, va="center")
        plt.axvline(x=exptime_ref, color='black', linestyle='--', alpha = 0.5)

        plt.text(
            0.02, 0.98,
            self._conditions_text(ra, dec, conditions),
            transform=plt.gca().transAxes,
            fontsize=12,
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round",
                facecolor="white",
                edgecolor="gray",
                alpha=0.9,
            ),
        )
        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("Exposure time (s)", fontsize=16)
        plt.ylabel("SNR", fontsize=16)
        plt.title(f"SNR vs Exposure Time (mag={magnitude})", fontsize=20)
        plt.grid(alpha=0.3)
        plt.legend(ncol=2, loc = 'lower right')
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.show()

    def plot_exptime_vs_snr(self,
                            filter: str or list[str],
                            magnitude=20.0,
                            snr_range=(1, 50),
                            snr_ref = 5,
                            n_points=100,
                            ra=None,
                            dec=None,
                            obsdate=None,
                            moon_phase=None,
                            moon_separation=None,
                            hour_since_sunset=None,
                            seeing=None,
                            save_path: str = None):
        self._check_model()
        conditions = self._resolve_conditions(ra, dec, obsdate, moon_phase, moon_separation, hour_since_sunset, seeing, verbose=False)
        filters = self._filter_list(filter)

        snrs = np.logspace(
            np.log10(snr_range[0]),
            np.log10(snr_range[1]),
            n_points,
        )
        snrs = np.append(snrs, snr_ref)

        tbl = Table()
        tbl["snr"] = snrs
        for filt in filters:
            tbl[filt] = np.nan
        for i, s in enumerate(snrs):
            result = self.calculate_exptime(
                filter=filters,
                magnitude=magnitude,
                snr=s,
                moon_phase=conditions['moon_phase'],
                moon_separation=conditions['moon_separation'],
                hour_since_sunset=conditions['hour_since_sunset'],
                seeing=conditions['seeing'],
                verbose=False,
            )
            for filt in filters:
                tbl[filt][i] = result[filt]["t_exp"]

        plt.figure(figsize=(8, 12))
        xmax = snrs[-2]
        plt.xlim(snrs[0], xmax * 2)
        colors = self.get_filter_colors(filters)
        for filt in filters:
            plt.plot(snrs, tbl[filt], label=f'{filt} (t_exp(SNR={snr_ref}) = {tbl[filt][-1]:.1f}s)', c = colors[filt])
            plt.text(snrs[-2] * 1.05, tbl[filt][-2], filt, color=colors[filt], fontsize=12, va="center")
        plt.axvline(x=snr_ref, color='black', linestyle='--', alpha = 0.5)

        plt.text(
            0.02, 0.98,
            self._conditions_text(ra, dec, conditions),
            transform=plt.gca().transAxes,
            fontsize=12,
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round",
                facecolor="white",
                edgecolor="gray",
                alpha=0.9,
            ),
        )

        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("SNR", fontsize=16)
        plt.ylabel("Exposure time (s)", fontsize=16)
        plt.title(f"Exposure time vs SNR (mag={magnitude})", fontsize=20)
        plt.grid(alpha=0.3)
        plt.legend(ncol=2, loc = 'lower right')
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.show()

    def plot_depth_model(self, exptime = 100, save_path: str = None):
        """Violin plot of the measured UL5 distribution per filter (needs self.db_tbl, see load_depth_table)."""
        if self.db_tbl is None or len(self.db_tbl) == 0:
            raise RuntimeError("Depth table is empty. Run update_depth_model or load_depth_table first.")

        filters = self.sort_filters_broad_then_medium(set(self.db_tbl["filter"]))
        colors = self.get_filter_colors(filters)

        # ------------------------------
        # 1. Collect data
        # ------------------------------
        data = []
        valid_filters = []

        for f in filters:
            m = (self.db_tbl["filter"] == f) & (self.db_tbl["exptime"] == exptime)
            ul5 = np.asarray(self.db_tbl["ul5"][m], dtype=float)
            ul5 = ul5[np.isfinite(ul5) & (ul5 > 0)]

            if len(ul5) < 5:
                continue

            data.append(ul5)
            valid_filters.append(f)

        positions = np.arange(len(valid_filters))

        # ------------------------------
        # 2. Plot
        # ------------------------------
        fig, ax = plt.subplots(figsize=(9, 6))

        vp = ax.violinplot(
            data,
            positions=positions,
            widths=0.8,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        # ------------------------------
        # 3. Color each violin
        # ------------------------------
        for body, f, ul5 in zip(vp["bodies"], valid_filters, data):
            verts = body.get_paths()[0].vertices
            ymin = np.min(ul5)
            ymax = np.max(ul5)
            verts[:, 1] = np.clip(verts[:, 1], ymin, ymax)
            body.set_facecolor(colors[f])
            body.set_alpha(0.7)

        # Median line style
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(2)

        # ------------------------------
        # 4. Model point at the fitted median conditions + median label
        # ------------------------------
        for i, (f, ul5) in enumerate(zip(valid_filters, data)):
            model = self.model.get(f) if self.model is not None else None
            if model is not None:
                beta0 = model["intercept"]
                beta = np.asarray(model["coefficients"])
                med = model["medians"]

                ul5_model = (
                    beta0
                    + beta[0] * med[0]
                    + beta[1] * med[1]
                    + beta[2] * med[2]
                    + beta[3] * med[3]
                )

                ax.scatter(
                    i,
                    ul5_model,
                    s=70,
                    color="white",
                    edgecolor="black",
                    zorder=5,
                )

            # robust upper edge (avoid tails)
            ul5_top = np.nanpercentile(ul5, 99)
            ul5_median = np.nanpercentile(ul5, 50)

            ax.text(
                i,
                ul5_top + 0.55,
                f"{ul5_median:.2f}",
                ha="center",
                va="bottom",
                fontsize=11,
                rotation=90,
                color=colors[f],
                clip_on=False,
                zorder=6,
            )
        # ------------------------------
        # 5. Formatting
        # ------------------------------
        ax.set_xticks(positions)
        ax.set_xticklabels(valid_filters, fontsize=12, rotation = 45)
        ax.tick_params(axis="y", labelsize=14)
        ax.set_ylabel("Limiting magnitude (5σ)", fontsize=16)
        ax.set_title(f"Limiting magnitude distribution (exptime={exptime}s)", fontsize=16)

        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.show()


#%%
if __name__ == "__main__":
    calc = ExposureCalculator()
    calc.load_model()

    ra, dec = 30.96774, 4.24528
    obsdate = Time.now() + 5 * u.hour

    # Exposure time to reach mag 19.5 at SNR = 10 in g, r, i for this target / time
    result_exptime = calc.calculate_exptime(filter=['g', 'r', 'i'], magnitude=19.5, snr=10, ra=ra, dec=dec, obsdate=obsdate)
    for filt, res in result_exptime.items():
        print(f"{filt}: t_exp = {res['t_exp']:.1f} s (UL5 at 100 s = {res['ul5_ref']:.2f})")

    # Expected SNR of a mag 19.5 source in 300 s under typical conditions with 1.5" seeing
    result_snr = calc.calculate_snr(filter='all', magnitude=19.5, exptime=300, seeing=1.5)
    for filt, res in result_snr.items():
        print(f"{filt}: SNR = {res['snr']:.1f}")

    calc.plot_exptime_vs_snr(filter='all', magnitude=19.5, ra=ra, dec=dec, obsdate=obsdate)
    calc.plot_snr_vs_exptime(filter='all', magnitude=19.5, ra=ra, dec=dec, obsdate=obsdate)
# %%
