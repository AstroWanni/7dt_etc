"""
7DT exposure time calculator with an empirical noise model (GAIN 2750, 1x1 binning).

Conditions -> depth model (UL5 of a 100 s frame, aperture D) -> per-pixel background
    UL5 = b0 + b1*moon_phase + b2*moon_sep + b3*airmass + b4*seeing
    sigma_ap = 10^(-0.4 (UL5 - ZP)) * t_ref / 5,  var_pix = sigma_ap^2 / N_ref,  B = (var_pix - RN^2) / t_ref
Total exposure T taken in N frames (read noise counted once per frame):
    SNR = S*T / sqrt( S*T + Npix * (B*T + N*RN^2) + (f_sys * B_sky*T*Npix)^2 )
    single exposure : N = 1
    stacked frames  : T = N * frame_time (N frames of 100 s) -> SNR = sqrt(N) * SNR(one frame)
    given T only    : N = ceil(T / frame_time)  (calculate_snr / calculate_limit)
    point    : S = 10^(-0.4 (m  - ZP)),     Npix = pixels in the aperture
    extended : S = 10^(-0.4 (mu - ZP)) * A, Npix = A / pixscale^2   (mu [mag/arcsec^2], A [arcsec^2])
"""
import json
import warnings
from importlib import resources
from pathlib import Path
import numpy as np
from . import constants as C
from . import plot as _plot
from .conditions import resolve_conditions


def _version():
    from . import __version__
    return __version__


def load_depth_model(aperture=C.DEFAULT_APERTURE):
    """Packaged depth model for the aperture diameter [arcsec] (10 or 5)."""
    if aperture not in C.APERTURE_MODELS:
        raise ValueError(f'aperture must be one of {list(C.APERTURE_MODELS)}')
    return json.loads(resources.files('sevendt_etc.models').joinpath(C.APERTURE_MODELS[aperture]).read_text())


class ExposureCalculator:
    def __init__(self, aperture=C.DEFAULT_APERTURE, model=None, rn=C.READ_NOISE, dark=C.DARK_CURRENT,
                 pixscale=C.PIXSCALE, frame_time=C.FRAME_TIME, sky_sys_frac=0.0):
        """
        aperture     : 10 (default) or 5 ["]; selects the packaged depth model
        model        : path to a depth model JSON (overrides aperture)
        rn, dark     : read noise [e-], dark current [e-/s/pix]
        pixscale     : ["/pix]
        frame_time   : single-frame exposure [s]
        sky_sys_frac : fractional sky systematic for extended sources (0 = pure random noise)
        """
        self.depth = json.loads(Path(model).read_text()) if model else load_depth_model(aperture)
        aps = {v['aperture'] for v in self.depth.values()}
        self.aperture = float(aps.pop()) if len(aps) == 1 else float(aperture)
        self.rn, self.dark, self.pixscale, self.frame_time, self.sky_sys_frac = rn, dark, pixscale, frame_time, sky_sys_frac
        self.n_ref = np.pi * (self.aperture / 2 / pixscale) ** 2
        self.features = {k: [{'moon_sep': 'moon_separation'}.get(x, x) for x in v['features']] for k, v in self.depth.items()}
        self.available = sorted(self.depth, key=C.band_key)
        self.default_bands = [f for f in C.DEFAULT_BANDS if f in self.depth]

    # ---------------- filters ----------------
    def filter_list(self, filter):
        """'all' -> default bands; 'available' -> every band in the model; str or list -> as given."""
        if isinstance(filter, str):
            if filter == 'all':
                return list(self.default_bands)
            if filter == 'available':
                return list(self.available)
            filter = [filter]
        filters = [str(f) for f in filter]
        bad = [f for f in filters if f not in self.depth]
        if bad:
            raise KeyError(f'{bad} not in the model. Available: {self.available}')
        extra = [f for f in filters if f not in C.DEFAULT_BANDS]
        if extra:
            warnings.warn(f'{extra}: {C.EXTRA_BANDS_NOTE}', stacklevel=3)
        return filters

    def _broadcast(self, filter, *others):
        arrs = [np.atleast_1d(a) for a in (self.filter_list(filter), *others)]
        n = max(len(a) for a in arrs)
        if any(len(a) not in (1, n) for a in arrs):
            raise ValueError('filter and value/exptime/snr must have the same length (or length 1)')
        arrs = [np.repeat(a, n) if len(a) == 1 else a for a in arrs]
        f = [str(x) for x in arrs[0]]
        keys = f if len(set(f)) == n else [f'{x}[{i}]' for i, x in enumerate(f)]
        return (keys, f, *arrs[1:])

    # ---------------- background from the depth model ----------------
    def band_params(self, filt, cond):
        """UL5 (100 s) under the conditions -> sky + dark rate per pixel [e-/s/pix]."""
        d = self.depth[filt]
        ul5 = float(d['intercept'] + np.dot(d['coefficients'], [cond[k] for k in self.features[filt]]))
        t_ref, zp = float(d['exptime']), float(d['zp_e_s'])
        var_pix = (10 ** (-0.4 * (ul5 - zp)) * t_ref / 5) ** 2 / self.n_ref
        bkg = max(var_pix - self.rn ** 2, 0.0) / t_ref
        sky = max(bkg - self.dark, 0.0)
        skysig = float(np.sqrt(bkg * self.frame_time + self.rn ** 2))   # per pixel, one frame (sky + dark + read)
        return {'ul5_ref': ul5, 't_ref': t_ref, 'zp_e_s': zp, 'bkg_rate': bkg, 'sky_rate': sky,
                'mu_sky': zp - 2.5 * np.log10(sky / self.pixscale ** 2) if sky > 0 else np.nan,
                'skysig_e': skysig, 'skysig_adu': skysig / C.EGAIN,
                'rn_frac_ref': min(self.rn ** 2 / var_pix, 1.0) if var_pix > 0 else 1.0,
                'floor_hit': bool(var_pix <= self.rn ** 2 + self.dark * t_ref)}

    def geometry(self, mode, area):
        """(Npix, area factor) for 'point' (aperture) or 'extended' (area [arcsec^2])."""
        if mode == 'point':
            return self.n_ref, 1.0
        if mode == 'extended':
            return area / self.pixscale ** 2, area
        raise ValueError("mode must be 'point' or 'extended'")

    # ---------------- frames ----------------
    def n_frame(self, T):
        """Number of frames for total exposure T: ceil(T / frame_time), at least 1."""
        return np.maximum(1, np.ceil(np.asarray(T, float) / self.frame_time - 1e-9)).astype(int)

    # ---------------- noise model ----------------
    def var_components(self, S, p, npix, T, n_frame=None):
        """Variance terms for total exposure T in n_frame frames (default ceil(T / frame_time); 1 = single exposure)."""
        T = np.asarray(T, float)
        n = self.n_frame(T) if n_frame is None else np.asarray(n_frame)
        return {'source': S * T, 'sky': npix * p['sky_rate'] * T, 'dark': npix * (p['bkg_rate'] - p['sky_rate']) * T,
                'read': npix * n * self.rn ** 2,
                'systematic': (self.sky_sys_frac * p['sky_rate'] * T * npix) ** 2}

    def snr(self, S, p, npix, T, n_frame=None):
        T = np.asarray(T, float)
        return S * T / np.sqrt(sum(self.var_components(S, p, npix, T, n_frame).values()))

    def exptime_single(self, S, p, npix, snr):
        """Single-exposure time [s] reaching snr (read noise once)."""
        if self.sky_sys_frac == 0:
            k = snr ** 2 * (S + npix * p['bkg_rate'])
            return float((k + np.sqrt(k ** 2 + 4 * S ** 2 * snr ** 2 * npix * self.rn ** 2)) / (2 * S ** 2))
        T = np.logspace(-3, 9, 20000)
        ok = self.snr(S, p, npix, T, n_frame=1) >= snr
        return float(T[np.argmax(ok)]) if ok.any() else np.inf

    def n_frame_stack(self, S, p, npix, snr, n_max=10 ** 7):
        """Number of frame_time frames to stack for snr (SNR = sqrt(N) * one-frame SNR); -1 if unreachable."""
        f = lambda n: self.snr(S, p, npix, n * self.frame_time, n_frame=n)
        if f(n_max) < snr:
            return -1
        lo, hi = 0, n_max
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if f(mid) >= snr else (mid, hi)
        return hi

    def limit_flux(self, p, npix, T, snr, n_frame=None):
        """Flux S [e-/s] reaching snr in total exposure T (n_frame frames; default ceil(T / frame_time))."""
        T = np.asarray(T, float)
        V = sum(v for k, v in self.var_components(0.0, p, npix, T, n_frame).items() if k != 'source')
        return (snr ** 2 * T + np.sqrt(snr ** 4 * T ** 2 + 4 * T ** 2 * snr ** 2 * V)) / (2 * T ** 2)

    @staticmethod
    def snr_sqrt_t(mag, p, T):
        """Reference only: UL5 scaled as sqrt(t), no read/source noise (previous ETC)."""
        return 5 * 10 ** (-0.4 * (mag - p['ul5_ref'] - 1.25 * np.log10(np.asarray(T, float) / p['t_ref'])))

    # ---------------- public calculations ----------------
    def calculate_snr(self, filter, value, exptime, mode='point', area=1.0, verbose=True, **cond_kw):
        """
        SNR for total exposure 'exptime' [s].
        value : AB mag (point) or surface brightness [mag/arcsec^2] (extended); scalar or per filter
        area  : measurement area [arcsec^2] (extended)
        cond_kw : ra, dec, obsdate, moon_phase, moon_separation, airmass, seeing
        """
        cond = resolve_conditions(**cond_kw)
        npix, A = self.geometry(mode, area)
        keys, filters, vals, ts = self._broadcast(filter, value, exptime)
        res = {}
        for key, f, v, T in zip(keys, filters, vals, ts):
            p = self.band_params(f, cond)
            S = 10 ** (-0.4 * (v - p['zp_e_s'])) * A
            var = self.var_components(S, p, npix, T)
            tot = sum(var.values())
            res[key] = {'filter': f, 'value': float(v), 'exptime': float(T), 'n_frame': int(self.n_frame(T)),
                        'snr': float(self.snr(S, p, npix, T)),
                        'snr_sqrt_t': float(self.snr_sqrt_t(v, p, T)) if mode == 'point' else np.nan,
                        **{f'frac_{k}': float(x / tot) for k, x in var.items()}, **p}
        if verbose:
            self._print(res, cond, mode, area, 'snr')
        return res

    def calculate_exptime(self, filter, value, snr=5.0, mode='point', area=1.0, verbose=True, **cond_kw):
        """
        Exposure to reach 'snr' for 'value' (mag or mag/arcsec^2).
        exptime_single : one exposure [s]
        n_frame        : number of frame_time (100 s) frames to stack
        exptime        : stacked total = n_frame * frame_time [s]
        """
        cond = resolve_conditions(**cond_kw)
        npix, A = self.geometry(mode, area)
        keys, filters, vals, snrs = self._broadcast(filter, value, snr)
        res = {}
        for key, f, v, s in zip(keys, filters, vals, snrs):
            p = self.band_params(f, cond)
            S = 10 ** (-0.4 * (v - p['zp_e_s'])) * A
            n = self.n_frame_stack(S, p, npix, s)
            res[key] = {'filter': f, 'value': float(v), 'snr': float(s), 'exptime_single': self.exptime_single(S, p, npix, s),
                        'n_frame': n, 'exptime': n * self.frame_time if n > 0 else np.inf, **p}
        if verbose:
            print(self.exptime_report(res, cond, mode, area) + '\n')
        return res

    def calculate_limit(self, filter, exptime, snr=5.0, mode='point', area=1.0, verbose=True, **cond_kw):
        """Limiting magnitude (point) or surface brightness (extended) at 'snr' for total exposure 'exptime'."""
        cond = resolve_conditions(**cond_kw)
        npix, A = self.geometry(mode, area)
        res = {}
        for f in self.filter_list(filter):
            p = self.band_params(f, cond)
            S = self.limit_flux(p, npix, exptime, snr)
            res[f] = {'filter': f, 'limit': float(p['zp_e_s'] - 2.5 * np.log10(S / A)), 'exptime': float(exptime),
                      'snr': float(snr), 'n_frame': int(self.n_frame(exptime)), **p}
        if verbose:
            self._print(res, cond, mode, area, 'limit')
        return res

    # ---------------- text ----------------
    def conditions_text(self, cond, mode, area):
        lines = []
        if cond['ra'] is not None:
            lines.append(f"RA, Dec = {cond['ra']:.4f}, {cond['dec']:.4f}")
        if cond['obsdate'] is not None:
            lines.append(f"Obsdate = {cond['obsdate'].iso[:19]} UTC")
        lines += [f"Moon phase = {cond['moon_phase']:.2f}, Moon sep = {cond['moon_separation']:.1f} deg",
                  f"Airmass = {cond['airmass']:.2f}, Seeing = {cond['seeing']:.2f}\"",
                  f'Point source, {self.aperture:g}" aperture' if mode == 'point'
                  else f'Extended source, area = {area:g} arcsec^2' + (f', sky sys = {self.sky_sys_frac:g}' if self.sky_sys_frac else ''),
                  f'Frame = {self.frame_time:g} s']
        return '\n'.join(lines)

    def detector_text(self):
        """Detector setup common to every band."""
        return '\n'.join([f'Setting        : GAIN {C.GAIN_MODE}, 1x1 binning, pixel scale = {self.pixscale:g} "/pix',
                          f'Gain           : {C.EGAIN:g} e-/ADU',
                          f'Read noise     : {self.rn:.2f} e-/pix (per frame)',
                          f'Dark current   : {self.dark:g} e-/s/pix',
                          f'Frame          : {self.frame_time:g} s per frame (N_frame = number of frames)'])

    def source_text(self, res, mode, area):
        """Source type, brightness, SNR and measurement geometry."""
        vals, snrs = {r['value'] for r in res.values()}, {r['snr'] for r in res.values()}
        unit = 'mag (AB)' if mode == 'point' else 'mag/arcsec^2'
        val = f'{vals.pop():g} {unit}' if len(vals) == 1 else f'per filter [{unit}]'
        snr = f'{snrs.pop():g}' if len(snrs) == 1 else 'per filter'
        geo = (f'aperture diameter = {self.aperture:g}", enclosed pixels = {self.n_ref:.1f}' if mode == 'point' else
               f'area = {area:g} arcsec^2 ({area / self.pixscale ** 2:.1f} pixels)'
               + (f', sky sys = {self.sky_sys_frac:g}' if self.sky_sys_frac else ''))
        return '\n'.join([f'Source         : {mode} source, {val}, target SNR = {snr}', f'Measurement    : {geo}'])

    def cond_text(self, cond):
        """Observing conditions used for the calculation."""
        lines = []
        if cond['ra'] is not None:
            lines.append(f"Target         : RA = {cond['ra']:.4f}, Dec = {cond['dec']:.4f} deg")
        if cond['obsdate'] is not None:
            lines.append(f"Obsdate        : {cond['obsdate'].iso[:19]} UTC")
        lines += [f"Moon           : phase = {cond['moon_phase']:.2f}, separation = {cond['moon_separation']:.1f} deg",
                  f"Airmass        : {cond['airmass']:.2f}", f"Seeing         : {cond['seeing']:.2f} arcsec"]
        return '\n'.join(lines)

    def exptime_report(self, res, cond, mode, area, plan=None, warn_threshold=None):
        """
        Text report of calculate_exptime (detector, source, conditions, [observation plan], per-band table).
        sky    : sky surface brightness [mag/arcsec^2] and sky rate [e-/s/pix]
        skysig : background RMS per pixel in one frame (sky + dark + read noise) [e-] and [ADU]
        T_single : one exposure [s];  T_stack = N_frame x frame_time (stacked frames) [s]
        """
        bar = '=' * 78
        out = [bar, f'7DT Exposure Time Calculator (sevendt_etc {_version()})', bar,
               '[Detector]', self.detector_text(), '', '[Source]', self.source_text(res, mode, area), '',
               '[Conditions]', self.cond_text(cond), '']
        if plan:
            out += ['[Observation plan]', plan, '']
        vary = len({(r['value'], r['snr']) for r in res.values()}) > 1
        vh, vu = (f"{'value':>7} {'SNR':>5} ", ' ' * 14) if vary else ('', '')
        ft = f'{self.frame_time:g}'
        out += [f'[Exposure time per band]  (skysig: per pixel, one {ft} s frame)',
                f'  T_single : one exposure;  T_stack = N_frame x {ft} s frames (SNR = sqrt(N_frame) x one-frame SNR)',
                f"{'filter':6} {vh}{'sky':>9} {'sky':>10} {'skysig':>7} {'skysig':>7} {'T_single':>10} {'T_stack':>10} {'N_frame':>8}",
                f"{'':6} {vu}{'[mag/as2]':>9} {'[e-/s/pix]':>10} {'[e-]':>7} {'[ADU]':>7} {'[s]':>10} {'[s]':>10} {f'[{ft} s]':>8}"]
        fmt = lambda x, w: f'{x:{w}.1f}' if np.isfinite(x) else f"{'-':>{w}}"
        for r in res.values():
            v = f"{r['value']:7.2f} {r['snr']:5.1f} " if vary else ''
            n = f"{r['n_frame']:8d}" if r['n_frame'] > 0 else f"{'-':>8}"
            flag = '  !' if warn_threshold is not None and r['exptime'] > warn_threshold else ''
            out.append(f"{r['filter']:6} {v}{r['mu_sky']:9.2f} {r['sky_rate']:10.4f} {r['skysig_e']:7.2f} "
                       f"{r['skysig_adu']:7.2f} {fmt(r['exptime_single'], 10)} {fmt(r['exptime'], 10)} {n}{flag}")
        if warn_threshold is not None and any(r['exptime'] > warn_threshold for r in res.values()):
            out.append(f'! : T_stack longer than {warn_threshold:g} s')
        bad = [r['filter'] for r in res.values() if r['floor_hit']]
        if bad:
            out.append(f'WARNING {bad}: predicted background below the read-noise floor (conditions outside the training range?)')
        return '\n'.join(out + [bar])

    def _print(self, res, cond, mode, area, kind):
        print(self.conditions_text(cond, mode, area))
        unit = 'mag' if mode == 'point' else 'mag/arcsec^2'
        if kind == 'snr':
            print(f"{'filter':6} {'value':>7} {'T[s]':>9} {'SNR':>8} {'N_frame':>8}")
            for r in res.values():
                print(f"{r['filter']:6} {r['value']:7.2f} {r['exptime']:9.1f} {r['snr']:8.2f} {r['n_frame']:8d}")
        else:
            print(f"{'filter':6} {'T[s]':>9} {'SNR':>5} {'limit':>7} {'N_frame':>8}")
            for r in res.values():
                print(f"{r['filter']:6} {r['exptime']:9.1f} {r['snr']:5.1f} {r['limit']:7.2f} {r['n_frame']:8d}")
        print(f'(value / limit in {unit})')
        bad = [r['filter'] for r in res.values() if r['floor_hit']]
        if bad:
            print(f'WARNING {bad}: predicted background below the read-noise floor (conditions outside the training range?)')
        print()

    # ---------------- plots ----------------
    def plot_snr_vs_exptime(self, filter, value, **kw):
        return _plot.plot_snr_vs_exptime(self, filter, value, **kw)

    def plot_exptime_vs_snr(self, filter, value, **kw):
        return _plot.plot_exptime_vs_snr(self, filter, value, **kw)

    def plot_limit_vs_exptime(self, filter, **kw):
        return _plot.plot_limit_vs_exptime(self, filter, **kw)

    def plot_noise_budget(self, filt, value, **kw):
        return _plot.plot_noise_budget(self, filt, value, **kw)
