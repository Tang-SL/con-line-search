"""Blended features: candidate inventory, shared-profile decomposition, component maps.

A broad ESPADA feature can contain several transitions (for example HCN 3-2 and
HCN v2=1 l=1e, 40 km/s apart). The centroid of such a feature has no meaning, so
instead of matching one catalogue line to it:

1. inventory() lists every catalogue transition whose predicted position falls
   inside the feature (plus a margin);
2. decompose() fits the aperture spectrum with one copy of a shared line shape
   per candidate transition. Lines from one nucleus share its kinematics, so a
   single shape (one Gaussian, two Gaussians, or the measured profile of another
   line) is placed at each candidate's predicted frequency; only the amplitudes
   differ. Amplitudes can be tied (A_a = ratio * A_b), e.g. HCN v2=1 l=1e to
   l=1f in another spectral window. Each line is tested by refitting without it
   (Delta BIC);
3. component_maps() repeats the fit in every pixel with the shape fixed (a linear
   least-squares problem), giving one map per transition.

Velocities are radio-convention offsets from each line's predicted frequency
f_e = rest/(1+z): v = c (f_e - f) / f_e (positive = redshifted), the same
convention as the velocity offsets elsewhere in the app.
"""
import copy
import itertools
import numpy as np
from scipy.optimize import least_squares, lsq_linear
from .core import C
from .measurements import matching_frequency
from .solutions import transition_key

PROFILES = {
    'gauss': 'One Gaussian',
    'gauss2': 'Two Gaussians (a negative second part models a self-absorption dip)',
    'template': 'Measured shape of another line (template)',
}


def label(row):
    return f"{row['species']} {row.get('qns', '')}".strip()+f" ({row['rest_GHz']:.5f} GHz)"


def channel_bounds(freq, lo, hi):
    half = abs(freq[1]-freq[0])/2
    a, b = sorted([float(freq[lo]), float(freq[hi])])
    return a-half, b+half


def inventory(feature, freq, rows, z, margin_kms=100.):
    """Every transition predicted inside the feature's channel range (+/- margin)."""
    low, high = channel_bounds(freq, feature['lo'], feature['hi'])
    pad = (low+high)/2*margin_kms/C
    observed = matching_frequency(feature)
    out, seen = [], set()
    for r in rows:
        expected = r['rest_GHz']/(1+z)
        key = transition_key(r)
        if low-pad <= expected <= high+pad and key not in seen:
            seen.add(key)
            out.append(dict(copy.deepcopy(r), expected_GHz=float(expected), inside=bool(low <= expected <= high),
                            velocity_offset_kms=float(C*(expected-observed)/expected)))
    return sorted(out, key=lambda r: r['expected_GHz'])


# ------------------------------------------------------------------ line shapes
def shape(kind, v, p, template=None):
    v = np.asarray(v, float)
    if kind == 'gauss':
        v0, s = p
        return np.exp(-.5*((v-v0)/s)**2)
    if kind == 'gauss2':
        v1, s1, v2, s2, r = p
        return np.exp(-.5*((v-v1)/s1)**2)+r*np.exp(-.5*((v-v2)/s2)**2)
    if kind == 'template':
        (shift,) = p
        return np.interp(v-shift, template['v'], template['y'], left=0., right=0.)
    raise ValueError(f'Unknown line shape {kind!r}.')


def shape_integral(kind, p, template=None):
    v = np.arange(-5000., 5000., .5)
    return float(np.sum(shape(kind, v, p, template))*.5)


def make_template(freq, spectrum, rest_GHz, z, lo, hi):
    """Normalised profile of a reference line in velocity (peak |y| = 1, positive integral)."""
    f = np.asarray(freq[lo:hi+1], float); y = np.asarray(spectrum[lo:hi+1], float)
    e = rest_GHz/(1+z); v = C*(e-f)/e
    order = np.argsort(v); v, y = v[order], y[order]
    total = np.trapezoid(y, v) if hasattr(np, 'trapezoid') else np.trapz(y, v)
    y = y*np.sign(total or 1.)/max(np.max(np.abs(y)), 1e-30)
    return dict(v=v.tolist(), y=y.tolist())


def _bounds_and_starts(kind, vmax, smin, smax, s0):
    if kind == 'gauss':
        lb, ub = [-vmax, smin], [vmax, smax]
        starts = [[v, s] for v in (0., -vmax/3, vmax/3) for s in (s0, max(smin, s0/2.5), min(smax, s0*2.5))]
    elif kind == 'gauss2':
        lb, ub = [-vmax, smin, -vmax, smin, -1.], [vmax, smax, vmax, smax, 3.]
        starts = [[0., s0, 0., max(smin, s0/3), -.5], [-s0/2, max(smin, s0/2), s0/2, max(smin, s0/2), 1.],
                  [0., s0, 0., min(smax, 2*s0), .3]]
    else:
        lb, ub = [-vmax], [vmax]
        starts = [[0.]]
    clip = lambda x: [min(max(a, l+1e-6), u-1e-6) for a, l, u in zip(x, lb, ub)]
    return (np.array(lb), np.array(ub)), [clip(s) for s in starts]


# ------------------------------------------------------------------ decomposition
def _columns(regions, components, kind, p, template):
    cols = []
    for c in components:
        parts = []
        for i, r in enumerate(regions):
            if i == c['region']:
                e = c['expected_GHz']; v = C*(e-np.asarray(r['freq'], float))/e
                parts.append(shape(kind, v, p, template))
            else:
                parts.append(np.zeros(len(r['freq'])))
        cols.append(np.concatenate(parts))
    return np.array(cols).T if cols else np.zeros((sum(len(r['freq']) for r in regions), 0))


def _reduce(X, components, baseline, regions):
    """Apply ties (A_a = ratio A_b) and add per-region offsets. Returns design, bounds, mapping."""
    index = {c['key']: j for j, c in enumerate(components)}
    free = [j for j, c in enumerate(components) if not c.get('tie')]
    Xr = [X[:, j].copy() for j in free]; lb, ub = [], []
    for j in free:
        sign = 1 if components[j]['sign'] == 'Emission' else -1
        lb.append(0. if sign > 0 else -np.inf); ub.append(np.inf if sign > 0 else 0.)
    for j, c in enumerate(components):
        if c.get('tie'):
            target, ratio = c['tie']
            Xr[free.index(index[target])] += ratio*X[:, j]
    n_comp = len(Xr)
    if baseline:
        start = 0
        for r in regions:
            col = np.zeros(X.shape[0]); col[start:start+len(r['freq'])] = 1.; start += len(r['freq'])
            Xr.append(col); lb.append(-np.inf); ub.append(np.inf)
    Xr = np.array(Xr).T if Xr else np.zeros((X.shape[0], 0))
    return Xr, (np.array(lb), np.array(ub)), free, n_comp


def _solve(Xr, y, w, bounds):
    if Xr.shape[1] == 0:
        return np.zeros(0), float(np.sum((y*w)**2))
    res = lsq_linear(Xr*w[:, None], y*w, bounds=bounds, method='bvls')
    return res.x, float(np.sum(((y-Xr@res.x)*w)**2))


def decompose(regions, components, kind='gauss', z=0., template=None, vmax=300., baseline=True, sigma0=None, assess=True, starts_override=None):
    """Fit the regions' spectra with a shared line shape placed at every component.

    regions: [dict(id, spw, freq, spectrum, noise)] (channel ranges of one or more features)
    components: [dict(key, row, expected_GHz, region, sign, tie=None or (key, ratio))]
    Returns amplitudes, uncertainties, integrated fluxes, Delta BIC per component and models.
    assess=False skips the per-line refits (Delta BIC, verdicts), e.g. for noise Monte Carlo runs;
    starts_override replaces the starting shapes (e.g. the best-fit shape of the observed data).
    """
    if not components:
        raise ValueError('Choose at least one line to fit.')
    for c in components:
        if c.get('tie'):
            target, ratio = c['tie']
            partner = next((d for d in components if d['key'] == target), None)
            if partner is None or partner.get('tie'):
                raise ValueError(f"{label(c['row'])} is tied to a line that is not in the fit or is itself tied.")
            if partner['sign'] != c['sign'] or not ratio > 0:
                raise ValueError('Tied lines must have the same sign and a positive ratio.')
    y = np.concatenate([np.asarray(r['spectrum'], float) for r in regions])
    w = 1/np.concatenate([np.asarray(r['noise'], float) for r in regions])
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(w))):
        raise ValueError('The spectrum contains blank channels in the fitted range.')
    dv = min(C*abs(r['freq'][1]-r['freq'][0])/r['freq'][0] for r in regions)
    width = max(C*abs(r['freq'][-1]-r['freq'][0])/r['freq'][0] for r in regions)
    s0 = sigma0 or width/6
    bounds, starts = _bounds_and_starts(kind, vmax, dv/2, width, s0)
    if kind == 'template' and template is None:
        raise ValueError('Choose a template line.')

    def fit_shape(comps, starts_):
        def residual(p):
            X = _columns(regions, comps, kind, p, template)
            Xr, b, _, _ = _reduce(X, comps, baseline, regions)
            a, _ = _solve(Xr, y, w, b)
            return (y-Xr@a)*w
        best = None
        for p0 in starts_:
            fit = least_squares(residual, p0, bounds=bounds, x_scale='jac')
            if best is None or fit.cost < best.cost:
                best = fit
        return best.x

    p = fit_shape(components, starts_override or starts)
    X = _columns(regions, components, kind, p, template)
    Xr, b, free, n_comp = _reduce(X, components, baseline, regions)
    a, chi2 = _solve(Xr, y, w, b)
    n = len(y); k = len(a)+len(p)
    # Uncertainties from the full model (amplitudes, offsets and shape together), so
    # that lines which can trade flux through the shared velocity/width show it.
    def model(theta):
        aa, pp = theta[:len(a)], theta[len(a):]
        Xp, _, _, _ = _reduce(_columns(regions, components, kind, pp, template), components, baseline, regions)
        return Xp@aa*w
    theta = np.r_[a, p]; base = model(theta); J = np.zeros((n, len(theta)))
    for t in range(len(theta)):
        step = 1e-6*max(abs(theta[t]), 1e-3) if t < len(a) else 1e-4*max(abs(theta[t]), 1.)
        tt = theta.copy(); tt[t] += step; J[:, t] = (model(tt)-base)/step
    cov = np.linalg.pinv(J.T@J)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    corr = cov/np.outer(np.where(err > 0, err, 1), np.where(err > 0, err, 1))
    integral = shape_integral(kind, p, template)
    amp = np.zeros(len(components)); amp_err = np.zeros(len(components))
    for i, j in enumerate(free):
        amp[j], amp_err[j] = a[i], err[i]
    index = {c['key']: j for j, c in enumerate(components)}
    for j, c in enumerate(components):
        if c.get('tie'):
            target, ratio = c['tie']
            amp[j], amp_err[j] = ratio*amp[index[target]], ratio*amp_err[index[target]]
    free_pos = {components[j]['key']: i for i, j in enumerate(free)}
    def pos(c):
        return free_pos[c['tie'][0] if c.get('tie') else c['key']]
    # Is each line needed? Refit everything (shape included) without it.
    results = []
    for j, c in enumerate(components):
        dropped = {c['key']} | ({d['key'] for d in components if d.get('tie') and d['tie'][0] == c['key']} if not c.get('tie') else set())
        rest = [d for d in components if d['key'] not in dropped]
        if not assess:
            chi2_0 = None
        elif rest:
            p0 = fit_shape(rest, [list(np.clip(p, bounds[0]+1e-6, bounds[1]-1e-6))])
            Xr0, b0, _, _ = _reduce(_columns(regions, rest, kind, p0, template), rest, baseline, regions)
            _, chi2_0 = _solve(Xr0, y, w, b0)
        else:
            Xr0, b0, _, _ = _reduce(np.zeros((n, 0)), [], baseline, regions)
            _, chi2_0 = _solve(Xr0, y, w, b0)
        delta_bic = chi2_0-chi2-(0 if c.get('tie') else np.log(n)) if chi2_0 is not None else None
        verdict = ('not tested' if delta_bic is None else 'Needed' if delta_bic > 10 else 'Weak evidence' if delta_bic > 2 else 'Not needed')
        degenerate = [d['key'] for d in components if d['key'] != c['key'] and pos(d) != pos(c)
                      and corr[pos(c), pos(d)] < -.8 and d['region'] == c['region']]
        results.append(dict(key=c['key'], row=c['row'], region=c['region'], sign=c['sign'], tie=c.get('tie'),
                            amplitude=float(amp[j]), amplitude_error=float(amp_err[j]),
                            flux=float(amp[j]*integral), flux_error=float(amp_err[j]*abs(integral)),
                            delta_bic=None if delta_bic is None else float(delta_bic), verdict=verdict, degenerate_with=degenerate))
    # Lines that can trade flux: test the group as a whole.
    def chi2_without(dropped):
        rest = [d for d in components if d['key'] not in dropped and not (d.get('tie') and d['tie'][0] in dropped)]
        if rest:
            p0 = fit_shape(rest, [list(np.clip(p, bounds[0]+1e-6, bounds[1]-1e-6))])
            Xr0, b0, _, _ = _reduce(_columns(regions, rest, kind, p0, template), rest, baseline, regions)
        else:
            Xr0, b0, _, _ = _reduce(np.zeros((n, 0)), [], baseline, regions)
        return _solve(Xr0, y, w, b0)[1]
    for q in results:
        if assess and q['degenerate_with'] and q['verdict'] != 'Needed':
            group = {q['key'], *q['degenerate_with']}
            free_count = sum(1 for d in components if d['key'] in group and not d.get('tie'))
            if chi2_without(group)-chi2-free_count*np.log(n) > 10:
                q['verdict'] = 'Needed, but not separable'
    # A tied line is fixed by its partner: report that instead of a test it cannot pass on its own.
    verdicts = {q['key']: q['verdict'] for q in results}
    for q in results:
        if q['tie'] and assess:
            q['verdict'] = 'Set by tie' if verdicts[q['tie'][0]] != 'Not needed' else 'Not needed (partner not needed)'
    names = {'gauss': ['velocity', 'width'], 'gauss2': ['main velocity', 'main width', 'second velocity', 'second width', 'second height'],
             'template': ['shift']}[kind]
    span = bounds[1]-bounds[0]
    at_bounds = [nm for nm, x, lo_, hi_, sp in zip(names, p, bounds[0], bounds[1], span) if min(x-lo_, hi_-x) < 1e-3*sp]
    for i, r in enumerate(regions):
        total = sum(abs(q['flux']) for q in results if q['region'] == i) or 1.
        for q in results:
            if q['region'] == i:
                q['fraction'] = abs(q['flux'])/total
    # Models per region for plotting.
    models, start = [], 0
    for i, r in enumerate(regions):
        sl = slice(start, start+len(r['freq'])); start += len(r['freq'])
        parts = {q['key']: (X[sl, j]*amp[j]).tolist() for j, q in enumerate(results) if q['region'] == i}
        offset = float(a[n_comp+i]) if baseline else 0.
        total = (np.sum([np.array(v) for v in parts.values()], axis=0) if parts else np.zeros(len(r['freq'])))+offset
        models.append(dict(components=parts, offset=offset, total=total.tolist(),
                           residual=(np.asarray(r['spectrum'])-total).tolist()))
    return dict(kind=kind, params=[float(x) for x in p], template=template, z=float(z), baseline=baseline, at_bounds=at_bounds, vmax=float(vmax),
                chi2=chi2, dof=int(n-k), reduced_chi2=float(chi2/max(n-k, 1)), n=n, bic=float(chi2+k*np.log(n)), n_parameters=int(k),
                components=results, models=models, integral_kms=integral,
                regions=[dict(id=r['id'], spw=r['spw'], lo=int(r['lo']), hi=int(r['hi'])) for r in regions])


def describe_shape(fit):
    p = fit['params']
    fwhm = 2*np.sqrt(2*np.log(2))
    if fit['kind'] == 'gauss':
        return f"shared velocity {p[0]:+.1f} km/s, FWHM {p[1]*fwhm:.1f} km/s"
    if fit['kind'] == 'gauss2':
        return (f"main part {p[0]:+.1f} km/s, FWHM {p[1]*fwhm:.1f} km/s; second part {p[2]:+.1f} km/s, "
                f"FWHM {p[3]*fwhm:.1f} km/s, relative height {p[4]:+.2f}"+(' (a dip)' if p[4] < 0 else ''))
    return f"template shifted by {p[0]:+.1f} km/s"


def assignments_from_fit(fit, accept=('Needed', 'Needed, but not separable', 'Set by tie')):
    """Per fitted feature: main line = strongest accepted component, blend = the others."""
    out = {}
    for i, r in enumerate(fit['regions']):
        kept = sorted([q for q in fit['components'] if q['region'] == i and q['verdict'] in accept],
                      key=lambda q: -abs(q['flux']))
        rows = [q['row'] for q in kept]
        out[r['id']] = dict(transition=copy.deepcopy(rows[0]) if rows else None,
                            blend=copy.deepcopy(rows[1:]), reviewed=True, origin='blend fit')
    return out


# ------------------------------------------------------------------ maps
def _resample(image, wcs_from, wcs_to, shape_to):
    """Nearest-pixel sample of image (on wcs_from) at the pixel centres of wcs_to."""
    yy, xx = np.indices(shape_to)
    world = wcs_to.pixel_to_world(xx.ravel(), yy.ravel())
    x, y = wcs_from.world_to_pixel(world)
    x, y = np.rint(x).astype(int), np.rint(y).astype(int)
    ok = (x >= 0) & (y >= 0) & (x < image.shape[1]) & (y < image.shape[0])
    out = np.full(len(x), np.nan)
    out[ok] = image[y[ok], x[ok]]
    return out.reshape(shape_to)


def component_maps(fit, cubes, noises, local=True, snr_cut=5., rows_per_chunk=32):
    """One integrated-intensity map per fitted line, from a per-pixel fit.

    cubes: region index -> Cube (the region's spectral window); noises: region index -> per-channel noise.
    In every pixel the lines are fitted together (strengths linear, plus an optional
    offset), exactly as for the aperture spectrum, in one of two ways:

    local=True (default): in pixels with signal (|integrated S/N| > snr_cut over the
      fitted channels) all lines share one shape whose velocity follows the pixel:
      either a Gaussian with its own velocity and width, or the aperture shape shifted
      to the pixel's velocity (keeps double-peaked, self-absorbed profiles at the
      nucleus), whichever fits better. With a fixed shape, a line shifted by the local
      velocity field would leak into its blended neighbours.
      The velocity is searched within half the smallest separation between the lines
      whose strengths are free (within that window two lines cannot swap places).
      With several regions, the gas has one velocity field: it is first measured in the
      least ambiguous region (widest window, ignoring tied lines) and every region is
      then fitted near that velocity, pixel by pixel.
      Other pixels use the aperture shape.
    local=False: every pixel uses the aperture shape, unshifted.

    Lines tied to a line in another region take ratio x that line's map (resampled).
    Returns ({component key: (image, region index)},
             {region index: dict(rms, velocity, sigma, note, limited, window_kms)}), where
    velocity is the per-pixel velocity and sigma the Gaussian width (NaN where the shifted
    aperture shape was used, outside the signal mask, or if local=False).
    Maps are signed and unclipped; units are the cube's unit x km/s.
    """
    from astropy.wcs import WCS
    comps = fit['components']; kind, p, template = fit['kind'], fit['params'], fit['template']
    integral = fit['integral_kms']
    by_key = {q['key']: q for q in comps}
    targets = {q['tie'][0] for q in comps if q.get('tie')}
    gauss_area = np.sqrt(2*np.pi)
    vv_ = np.arange(-5000., 5000., 1.); prof = shape(kind, vv_, p, template)
    v_ap = float(np.sum(vv_*prof)/np.sum(prof)) if np.sum(prof) else 0.
    flux = {}

    def onto(image, src_cube, dst_cube):
        ny, nx = dst_cube.data.shape[1:]
        if src_cube.data.shape[1:] == (ny, nx) and WCS(src_cube.header).celestial.to_header() == WCS(dst_cube.header).celestial.to_header():
            return image
        return _resample(image, WCS(src_cube.header).celestial, WCS(dst_cube.header).celestial, (ny, nx))

    def window(i, ignore_ties):
        here = [q for q in comps if q['region'] == i and not (q.get('tie') and (ignore_ties or by_key[q['tie'][0]]['region'] == i))]
        free_e = [q['row']['rest_GHz']/(1+fit['z']) for q in here if not q.get('tie')]
        seps = [abs(C*(a_-b_)/a_) for n_, a_ in enumerate(free_e) for b_ in free_e[n_+1:]]
        return (min(.5*min(seps), fit.get('vmax', 300.)) if seps else fit.get('vmax', 300.)), bool(seps)

    def region_pass(i, reference=None, ignore_external=False):
        r = fit['regions'][i]; cube = cubes[i]
        lo, hi = r['lo'], r['hi']
        freq = np.asarray(cube.frequency[lo:hi+1], float)
        ny, nx = cube.data.shape[1:]
        here = [q for q in comps if q['region'] == i]
        expected = {q['key']: q['row']['rest_GHz']/(1+fit['z']) for q in here}
        vel = {k: C*(e-freq)/e for k, e in expected.items()}
        free_keys = [q['key'] for q in here if not q.get('tie')]
        same = [q for q in here if q.get('tie') and by_key[q['tie'][0]]['region'] == i]
        other = [] if ignore_external else [q for q in here if q.get('tie') and by_key[q['tie'][0]]['region'] != i]
        w = 1/np.asarray(noises[i][lo:hi+1], float)
        partner = {q['key']: q['tie'][1]*onto(flux[q['tie'][0]], cubes[by_key[q['tie'][0]]['region']], cube) for q in other}

        def design(profile_of, norm):
            """Columns for one shape; ties in this region merged. profile_of(key) -> column."""
            cols = [profile_of(k) for k in free_keys]
            for q in same:
                j = free_keys.index(q['tie'][0]); cols[j] = cols[j]+q['tie'][1]*profile_of(q['key'])
            X = np.array(cols+([np.ones(len(freq))] if fit['baseline'] else [])).T
            fixed_cols = [(q['key'], profile_of(q['key'])/norm) for q in other]   # partner maps are fluxes
            return X, fixed_cols

        Xa, fa = design(lambda k: shape(kind, vel[k], p, template), integral)
        dv = C*abs(freq[1]-freq[0])/freq[0]; span = C*abs(freq[-1]-freq[0])/freq[0]
        sgrid = np.geomspace(max(1.5*dv, 5.), max(span/4, 3*dv), 24)
        half, has_seps = window(i, False)
        step = max(dv, 8.)
        offsets = np.arange(-half, half+1e-9, step) if half > step else np.array([0.])
        near = min(half, max(3*step, 30.))   # search half-width around a reference velocity
        near_offsets = np.arange(-near, near+1e-9, step) if near > step else np.array([0.])
        maps = {k: np.full((ny, nx), np.nan) for k in free_keys}
        rms_map = np.full((ny, nx), np.nan); vmap = np.full((ny, nx), np.nan); smap = np.full((ny, nx), np.nan)
        # Signal cut: per-pixel noise = robust spread of the aperture-shape fit residual in that
        # pixel; a sum over correlated channels is wider than sqrt(N) x that, so the S/N is
        # divided by a factor measured on the same number of line-free channels (if any).
        nchan = cube.data.shape[0]
        outside = [c for c in range(max(0, lo-40), lo)]+[c for c in range(hi+1, min(nchan, hi+41))]
        correlation = 1.
        if len(outside) >= 10:
            side = np.asarray(cube.data[outside], dtype=float)
            side_noise = 1.4826*np.nanmedian(np.abs(side-np.nanmedian(side, axis=0)), axis=0)
            with np.errstate(divide='ignore', invalid='ignore'):
                z_ = np.nansum(side, axis=0)/(side_noise*np.sqrt(len(outside)))
            spread = 1.4826*np.nanmedian(np.abs(z_-np.nanmedian(z_)))
            if np.isfinite(spread) and spread > 1:
                correlation = float(spread)
        collected = []
        for y0 in range(0, ny, rows_per_chunk):
            y1 = min(ny, y0+rows_per_chunk)
            D = np.asarray(cube.data[lo:hi+1, y0:y1, :], dtype=float).reshape(len(freq), -1)
            blank = ~np.all(np.isfinite(D), axis=0); D = np.where(np.isfinite(D), D, 0.)
            fixed = np.zeros_like(D)
            for k, col in fa:
                fixed += np.outer(col, np.nan_to_num(partner[k][y0:y1].ravel()))
            A = np.linalg.pinv(Xa*w[:, None])@((D-fixed)*w[:, None])
            resid = D-fixed-Xa@A
            rms = np.sqrt(np.mean(resid**2, axis=0)); rms[blank] = np.nan
            rms_map[y0:y1] = rms.reshape(y1-y0, nx)
            for j, k in enumerate(free_keys):
                a = A[j]*integral; a[blank] = np.nan
                maps[k][y0:y1] = a.reshape(y1-y0, nx)
            if local:
                pix_noise = 1.4826*np.median(np.abs(resid-np.median(resid, axis=0)), axis=0)
                with np.errstate(divide='ignore', invalid='ignore'):
                    snr = np.abs(np.sum(D, axis=0))/(pix_noise*np.sqrt(len(freq))*correlation)
                sig = (snr > snr_cut) & ~blank
                if sig.any():
                    collected.append((np.flatnonzero(sig)+y0*nx, D[:, sig]))
        limited = False
        if collected:
            pix = np.concatenate([c[0] for c in collected]); Ds = np.concatenate([c[1] for c in collected], axis=1)
            m = Ds.shape[1]; nf = len(free_keys); nb = Xa.shape[1]-nf
            signs = np.array([1 if by_key[k]['sign'] == 'Emission' else -1 for k in free_keys])
            subsets = [list(c) for n_ in range(nf+1) for c in itertools.combinations(range(nf), n_)] if nf <= 6 else [list(range(nf))]
            fixed_amp = {k: np.nan_to_num(partner[k].ravel()[pix]) for k in partner}
            bA = np.zeros((nf, m)); brms = np.zeros(m); bv = np.zeros(m); bs = np.zeros(m)
            # Pixels are grouped by the velocity they are searched around.
            ref = reference.ravel()[pix] if reference is not None else np.full(m, np.nan)
            centre = np.where(np.isfinite(ref), np.round(ref/step)*step, np.nan)
            groups = {}
            for n_, c_ in enumerate(centre):
                groups.setdefault(None if not np.isfinite(c_) else float(c_), []).append(n_)
            for c_, members in groups.items():
                members = np.array(members)
                vgrid = v_ap+offsets if c_ is None else c_+near_offsets
                Dg = Ds[:, members]; best = np.full(len(members), np.inf)
                for v0 in vgrid:
                    candidates = [(s0, (lambda k, v0=v0, s0=s0: shape('gauss', vel[k], (v0, s0))), s0*gauss_area) for s0 in sgrid]
                    candidates.append((np.nan, (lambda k, d=v0-v_ap: shape(kind, vel[k]-d, p, template)), integral))
                    for s0, profile_of, norm in candidates:
                        X, fc = design(profile_of, norm)
                        Dr = Dg.copy()
                        for k, col in fc:
                            Dr -= np.outer(col, fixed_amp[k][members])
                        Drw = Dr*w[:, None]
                        for subset in subsets:
                            cols = subset+list(range(nf, nf+nb))
                            Xc = X[:, cols]
                            Ag = np.linalg.pinv(Xc*w[:, None])@Drw
                            rg = Dr-Xc@Ag
                            cg = np.sum((rg*w[:, None])**2, axis=0)
                            if subset:   # strengths keep the line's sign, as in the aperture fit
                                cg = np.where(np.all(Ag[:len(subset)]*signs[subset][:, None] >= 0, axis=0), cg, np.inf)
                            better = cg < best
                            if not better.any():
                                continue
                            best = np.where(better, cg, best)
                            idx = members[better]
                            full = np.zeros((nf, len(idx)))
                            if subset:
                                full[subset] = Ag[:len(subset), better]*norm
                            bA[:, idx] = full
                            brms[idx] = np.sqrt(np.mean(rg[:, better]**2, axis=0)); bv[idx] = v0; bs[idx] = s0
            for j, k in enumerate(free_keys):
                maps[k].ravel()[pix] = bA[j]
            rms_map.ravel()[pix] = brms; vmap.ravel()[pix] = bv; smap.ravel()[pix] = bs
            # Evidence that the velocity field reaches beyond the window: many pixels at its edge
            # (only meaningful without a reference velocity).
            free_search = ~np.isfinite(centre)
            if has_seps and free_search.sum() >= 10:
                limited = bool((np.abs(bv[free_search]-v_ap) >= half-step/2).mean() > .1)
        if reference is not None:
            note = f'per-pixel velocity searched within ±{near:.0f} km/s of the velocity field measured in the least ambiguous window'
        elif has_seps:
            note = f'per-pixel velocity searched within ±{half:.0f} km/s of the aperture velocity (half the smallest separation between free lines)'
        else:
            note = 'per-pixel velocity searched freely (one free line)'
        out = {k: maps[k] for k in free_keys}
        for q in same:
            out[q['key']] = q['tie'][1]*maps[q['tie'][0]]
        for q in other:
            out[q['key']] = np.where(np.isnan(rms_map), np.nan, partner[q['key']])
        return out, dict(rms=rms_map, velocity=vmap, sigma=smap, note=note if local else 'aperture shape in every pixel',
                         limited=bool(local and limited), window_kms=float(near if reference is not None else half))

    n = len(fit['regions'])
    reference, ref_region = None, None
    if local and n > 1:
        ref_region = max(range(n), key=lambda i: (window(i, True)[0], sum(abs(q['flux']) for q in comps if q['region'] == i and not q.get('tie'))))
        _, d0 = region_pass(ref_region, ignore_external=True)
        reference = d0['velocity']
    order = sorted(range(n), key=lambda i: 0 if any(q['region'] == i and q['key'] in targets for q in comps) else 1)
    diag = {}
    for i in order:
        ref_i = onto(reference, cubes[ref_region], cubes[i]) if reference is not None else None
        out, diag[i] = region_pass(i, ref_i)
        flux.update(out)
    return {q['key']: (flux[q['key']], q['region']) for q in comps}, diag
