"""Static figures (PNG) for the weblog: every step gets plots that can be checked by eye."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ..core import C

TEAL, ORANGE, PURPLE, GREY, RED, BLUE = '#167f83', '#c2410c', '#7356ad', '#6b7280', '#b91c1c', '#1d4ed8'
COLOURS = ['#d95f02', '#7570b3', '#1b9e77', '#e7298a', '#1f78b4', '#66a61e', '#a6761d', '#e6ab02']    # distinct categorical colours
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'figure.dpi': 110, 'axes.spines.top': False, 'axes.spines.right': False})


def _save(fig, path):
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


FULL_VIEW_PIXELS = 64     # images this small (cutouts) are always shown whole


def _window(centre, zoom, nx, ny, box=None):
    """Zoom window around centre, kept inside the data: (cx, cy, half_x, half_y) in pixels.

    box = (x0, x1, y0, y1), the pixel range that holds finite data (default: the whole image).
    Small data areas (cutouts, <= FULL_VIEW_PIXELS on both axes) are shown whole. Larger ones get
    a square window that is never larger than the data and is shifted, not cut, at the edges, so a
    panel never shows blank space beyond the data."""
    x0, x1, y0, y1 = box if box else (0, nx-1, 0, ny-1)
    w, h = x1-x0+1, y1-y0+1
    if max(w, h) <= FULL_VIEW_PIXELS or not zoom or centre is None:
        return (x0+x1)/2, (y0+y1)/2, w/2, h/2
    half = min(max(12., min(nx, ny)/(2*zoom)), w/2, h/2)
    cx = float(np.clip(centre[0], x0+half-.5, x1+.5-half))
    cy = float(np.clip(centre[1], y0+half-.5, y1+.5-half))
    return cx, cy, half, half


def _finite_box(data):
    rows, cols = np.flatnonzero(np.isfinite(data).any(axis=1)), np.flatnonzero(np.isfinite(data).any(axis=0))
    return (cols[0], cols[-1], rows[0], rows[-1]) if rows.size else None


def _image(ax, image, title, centre=None, zoom=None, cmap='viridis', symmetric=False, label=None, marks=()):
    data = np.asarray(image, float)
    ny, nx = data.shape
    box = _finite_box(data)
    cx, cy, hx, hy = _window(centre, zoom, nx, ny, box)
    view = data[max(0, int(np.floor(cy-hy+.5))):int(np.ceil(cy+hy+.5)), max(0, int(np.floor(cx-hx+.5))):int(np.ceil(cx+hx+.5))]
    finite = view[np.isfinite(view)]
    if finite.size:
        if symmetric:            # colour scale from what is shown, so a compact source is not saturated
            v = np.nanpercentile(np.abs(finite), 99.8) or 1
            vmin, vmax = -v, v
        else:
            vmin, vmax = np.nanpercentile(finite, [1, 99.7])
    else:
        vmin, vmax = 0, 1
    h = ax.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_xlim(cx-hx, cx+hx); ax.set_ylim(cy-hy, cy+hy)
    for x, y, style in marks:
        ax.plot(x, y, style, ms=9, mfc='none', mew=1.6)
    ax.set_title(title, fontsize=9); ax.set_xlabel('x (pixel)'); ax.set_ylabel('y (pixel)')
    cb = plt.colorbar(h, ax=ax, fraction=.046, pad=.03)
    if label:
        cb.set_label(label)
    return h


def inspection(path, cube, info, median_map, peak_snr, z=None):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(width_ratios=[1, 1, 1.6]), layout='constrained')
    lx, ly = info['line_peak_pixel']; cx, cy = info['continuum_pixel']
    real = (info.get('continuum_snr') or 0) >= 5          # otherwise the 'peak' is noise: centre on the line source
    _image(ax[0], median_map, 'Continuum estimate (lower of the two band-edge means)' + ('' if real else f"\nno significant continuum (peak S/N {info.get('continuum_snr') or 0:.0f})"),
           (cx, cy) if real else (lx, ly), 3, label=info['unit'], marks=[(cx, cy, 'w+')] if real else [(lx, ly, 'c+')])
    _image(ax[1], peak_snr, 'Peak |S/N| over channels (line emission)', (lx, ly), 3, cmap='magma', marks=[(lx, ly, 'c+')])
    f = cube.frequency
    for (x, y), colour, name in [((lx, ly), TEAL, f'brightest line pixel {lx},{ly}'), ((cx, cy), GREY, f'continuum peak {cx},{cy}')]:
        ax[2].step(f, cube.data[:, y, x], where='mid', color=colour, lw=.8, label=name)
    ax[2].fill_between(f, -3*info['noise_per_channel'], 3*info['noise_per_channel'], color='#cbd5e1', alpha=.6, step='mid', label='±3σ per channel')
    ax[2].set_xlabel(f'Observed frequency (GHz, {info["frame"]})'); ax[2].set_ylabel(info['unit']); ax[2].legend(fontsize=7)
    ax[2].set_title(f"{info['key']}: single-pixel spectra")
    if z is not None:
        top = ax[2].secondary_xaxis('top', functions=(lambda x: x*(1+z), lambda x: x/(1+z)))
        top.set_xlabel(f'Rest frequency at z = {z:.5f} (GHz)')
    return _save(fig, path)


def coverage(path, rows, z=None):
    fig, ax = plt.subplots(figsize=(10, 1.4+.35*len(rows)))
    for i, r in enumerate(rows):
        lo, hi = (r['freq_min_GHz'], r['freq_max_GHz'])
        if z is not None:
            lo, hi = lo*(1+z), hi*(1+z)
        ax.barh(i, hi-lo, left=lo, color=TEAL, alpha=.7)
        ax.text(hi, i, '  '+r['key'], va='center', fontsize=8)
    ax.set_yticks([]); ax.set_xlabel('Rest frequency (GHz)' if z is not None else 'Observed frequency (GHz)')
    ax.set_title('Spectral coverage' + (f' (rest frame, z = {z:.5f})' if z is not None else ''))
    return _save(fig, path)


def components(path, image, components_, key, centre=None):
    fig, ax = plt.subplots(figsize=(6.2, 5.4), layout='constrained')
    _image(ax, image, f'{key}: ESPADA components (red = emission, blue = absorption)', centre, 2 if centre else None)
    for c in components_:
        colour = RED if c['sign'] == 'Emission' else BLUE
        ax.plot(c['x'], c['y'], 'o', mfc='none', mec=colour, ms=8, mew=1.3)
        ax.annotate(c['id'].rsplit('-', 1)[-1], (c['x'], c['y']), color=colour, fontsize=6, xytext=(3, 3), textcoords='offset points')
    return _save(fig, path)


def spectrum(path, freq, spec, noise, features=(), title='', unit='Jy/beam', z=None, marks=(), xlim=None, models=()):
    fig, ax = plt.subplots(figsize=(12, 3.6))
    ax.step(freq, spec, where='mid', color=TEAL, lw=.9, label='aperture spectrum')
    ax.fill_between(freq, -3*noise, 3*noise, color='#cbd5e1', alpha=.6, step='mid', label='±3σ (spatial controls)')
    top = np.nanmax(spec) if np.isfinite(np.nanmax(spec)) else 1
    for q in features:
        a, b = sorted([freq[q['lo']], freq[q['hi']]])
        good = q.get('significant') and q.get('reliable', True)
        ax.axvspan(a, b, color='#16a34a' if good else '#d97706', alpha=.18 if good else .1, lw=0)
        ax.text((a+b)/2, top*1.02, q['id'].rsplit('-', 1)[-1], ha='center', fontsize=7, color='#166534' if good else '#92400e', rotation=90, va='bottom')
    for x, label, colour in marks:
        ax.axvline(x, color=colour, ls='--', lw=.8)
        ax.text(x, ax.get_ylim()[0], ' '+label, rotation=90, fontsize=6, color=colour, va='bottom')
    for f_, y_, label, colour, style in models:
        ax.plot(f_, y_, style, color=colour, lw=1.2, label=label)
    ax.axhline(0, color='#94a3b8', lw=.5)
    ax.set_xlabel('Observed frequency (GHz)'); ax.set_ylabel(unit); ax.set_title(title)
    if xlim:
        ax.set_xlim(*xlim)
    if z is not None:
        sec = ax.secondary_xaxis('top', functions=(lambda x: x*(1+z), lambda x: x/(1+z)))
        sec.set_xlabel(f'Rest frequency (GHz), z = {z:.5f}')
    ax.legend(fontsize=7, loc='upper right')
    return _save(fig, path)


GREEN = '#15803d'
STATUS_STYLE = {'consistent': ('seen', GREEN), 'missing': ('MISSING', RED), 'detected': ('emission here, but predicted below noise', GREY),
                'not detected': ('below the noise: no test', GREY)}
VERDICT_COLOUR = {'rejected': RED, 'partial': '#b45309', 'supported': GREEN, 'inconclusive': '#b45309', 'untested': GREY}


def _template(feature, spectra):
    """The feature's own line shape: velocity axis (km/s, around its centroid) and a profile with unit
    integral, so a companion predicted at P Jy/beam km/s is drawn as P x profile."""
    sp = spectra[feature['spw']]; f = np.asarray(sp['frequency']); y = np.asarray(sp['spectrum'])
    sel = (f >= feature['freq_lo_GHz']) & (f <= feature['freq_hi_GHz'])
    if sel.sum() < 2:
        return None
    sign = 1 if feature['sign'] == 'Emission' else -1
    v = C*(feature['freq_GHz']-f[sel])/feature['freq_GHz']
    prof = np.clip(sign*y[sel], 0, None)
    order = np.argsort(v); v, prof = v[order], prof[order]
    area = float(np.sum(.5*(prof[1:]+prof[:-1])*np.diff(v)))
    return (v, prof/area) if area > 0 else None


def companions(path, feature, cand, spectra, max_panels=8):
    """Does the candidate's other lines appear where, and as strongly as, LTE predicts?

    Panel 1: the feature itself. Further panels: the observed spectrum around each of the
    candidate's other lines in the observed windows (velocity relative to the line, at the
    feature's velocity), with the line it predicts drawn on top. The predicted line has the
    feature's own shape, scaled by LTE (temperature in the title) as if the candidate made the
    whole feature (filled); dashed is the largest share of the feature it can make without
    predicting a missing line (f_max)."""
    sign = 1 if feature['sign'] == 'Emission' else -1
    rank = {'missing': 0, 'consistent': 1, 'detected': 3, 'not detected': 4}
    rows = sorted(cand['companions'], key=lambda r: (rank.get(r['status'], 2), -(r['predicted'] or 0)))
    shown = rows[:max_panels-1]
    n = 1+max(1, len(shown)); cols = min(4, n); nrows = int(np.ceil(n/cols))
    fig, axes = plt.subplots(nrows, cols, figsize=(3.6*cols, 2.9*nrows+.9), squeeze=False, layout='constrained')
    tmpl = _template(feature, spectra)
    width = feature['width_kms']; span = max(1.5*width, 300.)
    fmax = cand.get('fmax')
    dashed = fmax is not None and .05 <= fmax < .98

    def panel(ax, key, centre, title, pred=None, status=None, meas=None, err=None, here=False):
        sp = spectra[key]; f = np.asarray(sp['frequency']); y = np.asarray(sp['spectrum']); nz = np.asarray(sp['noise'])
        v = C*(centre-f)/centre
        sel = np.abs(v) <= span
        if sel.sum() < 2:
            ax.text(.5, .5, 'outside the data', transform=ax.transAxes, ha='center', color=GREY); ax.set_title(title, fontsize=8); return
        o = np.argsort(v[sel]); vv, yy, nn = v[sel][o], y[sel][o], nz[sel][o]
        ax.fill_between(vv, -3*nn, 3*nn, color='#e2e8f0', step='mid', lw=0)
        ax.axvspan(-width/2, width/2, color='#f1f5f9', zorder=0)
        ax.step(vv, yy, where='mid', color='#0f172a', lw=.8)
        ax.axhline(0, color='#94a3b8', lw=.5)
        model = None
        if pred is not None and tmpl is not None:
            model = sign*pred*np.interp(vv, tmpl[0], tmpl[1], left=0, right=0)
            ax.fill_between(vv, 0, model, color=PURPLE, alpha=.28, lw=0)
            ax.plot(vv, model, color=PURPLE, lw=1)
            if dashed:
                ax.plot(vv, fmax*model, color=PURPLE, lw=1.1, ls='--')
        if here:
            for t in cand['transitions']:
                if t.get('inside'):
                    ax.axvline(C*(centre-t['expected_GHz']*(1-cand['offset_kms']/C))/centre, color=PURPLE, lw=1, ls=':')
        ax.set_xlim(-span, span)
        noise3 = 3*float(np.nanmedian(nn))
        d_hi, d_lo = max(float(np.nanmax(yy)), noise3), min(float(np.nanmin(yy)), -noise3)
        peak = float(np.nanmax(np.abs(model))) if model is not None else 0.
        scale = max(abs(d_hi), abs(d_lo))
        if peak > 3*scale:          # a huge predicted line would flatten the data: keep the data readable
            top = 3*scale if sign > 0 else -d_lo*1.15
            ax.set_ylim(d_lo*1.15 if sign > 0 else -3*scale, top*1.3 if sign > 0 else d_hi*1.3)
            ax.text(.5, .80 if sign > 0 else .08, f"{'▲' if sign > 0 else '▼'} predicted peak {sign*peak:.3g} Jy/beam (off scale)",
                    transform=ax.transAxes, ha='center', fontsize=7, color=PURPLE, fontweight='bold')
        else:
            lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi+.25*(hi-lo))
        ax.set_title(title, fontsize=8)
        if status:
            label, colour = STATUS_STYLE.get(status, ('blended: another feature', '#b45309'))
            if status == 'missing' and meas is not None and err and meas > 3*err:
                label = 'TOO WEAK'                 # present, but well below the prediction
            if status.startswith('blended'):
                label = status.replace('blended with', 'on feature').replace(feature['spw']+'-', '')
            ax.text(.03, .95, label, transform=ax.transAxes, va='top', fontsize=8, color='white', fontweight='bold',
                    bbox=dict(facecolor=colour, edgecolor='none', pad=2.2))
            ax.text(.97, .95, f'Jy/beam km/s\npredicted {pred:.3g}\nobserved {meas:.2g} ± {err:.1g}', transform=ax.transAxes, va='top', ha='right', fontsize=6.5, color='#334155')
        ax.tick_params(labelsize=7)

    main = next((t for t in cand['transitions'] if t.get('inside')), cand['transitions'][0])
    panel(axes[0, 0], feature['spw'], feature['freq_GHz'], f"the feature: {feature['id']}\n{(main.get('qns') or '')[:30]}", here=True)
    for ax, r in zip(axes.flat[1:], shown):
        panel(ax, r['window'], r['centre_GHz'], f"{r['window']}: {(r.get('qns') or '')[:26]}\n{r['rest_GHz']:.4f} GHz, E_u {r['Eu_K'] or 0:.0f} K",
              r['predicted'], r['status'], r['measured'], r['error'])
    if not shown:
        ax = axes.flat[1]; ax.axis('off')
        ax.text(.5, .5, f"No other line of {cand['species']}\nfalls in the observed windows:\nnothing to test it against.",
                transform=ax.transAxes, ha='center', va='center', fontsize=9, color=GREY,
                bbox=dict(facecolor='#f8fafc', edgecolor='#cbd5e1', boxstyle='round,pad=.8'))
    for ax in list(axes.flat)[n:]:
        ax.axis('off')
    for ax in axes[-1]:
        ax.set_xlabel('velocity relative to the line (km/s)', fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel('Jy/beam', fontsize=7)
    v, why = cand.get('verdict') or 'untested', cand.get('verdict_reason') or ''
    T = cand.get('best_T'); Tr = cand.get('T_range')
    temp = '' if T is None else (f'LTE, T = {T:g} K' + (f' (same result for {Tr[0]:g}–{Tr[1]:g} K)' if Tr and Tr[0] != Tr[1] else ''))
    left = len(rows)-len(shown)
    fig.suptitle(f"Could {feature['id']} be {cand['species']}?  If so, its other lines should appear as drawn in purple ({temp})", fontsize=10)
    fig.text(.01, -.01, f"{v.upper()}: {why}." + (f"  It can make at most {fmax:.0%} of the feature (f_max)." if fmax is not None and v not in ('untested', 'partial') else '')
             + (f"  {left} more weaker lines not shown." if left > 0 else '')
             + "\nGrey band: ±3σ noise; light column: integration window; purple: predicted line (feature's shape) if "
               "the candidate made the whole feature" + (f"; dashed: the largest share it can make (f_max = {fmax:.2f})." if dashed else '.'),
             fontsize=7.5, color=VERDICT_COLOUR.get(v, GREY), va='top')
    return _save(fig, path)


def strategy(path, freq, spec, noise, features, lines, title, z=None, xlim=None):
    """Aperture spectrum with the strong lines of the species to be tested.

    lines: list of (observed GHz, species, transition label, colour). One legend entry per
    species; transition labels sit at the top of the panel."""
    fig, ax = plt.subplots(figsize=(12, 4.2), layout='constrained')
    ax.fill_between(freq, -3*noise, 3*noise, color='#e2e8f0', alpha=.8, step='mid', label='±3σ')
    for q in features:
        a, b = sorted([freq[q['lo']], freq[q['hi']]])
        ax.axvspan(a, b, color='#16a34a' if q.get('significant') and q.get('reliable', True) else '#d97706', alpha=.08, lw=0)
    ax.step(freq, spec, where='mid', color='#0f172a', lw=.8, label='aperture spectrum')
    ax.axhline(0, color='#94a3b8', lw=.5)
    lo, hi = np.nanmin(spec), np.nanmax(spec)
    ax.set_ylim(min(lo, -3*np.nanmedian(noise))*1.1, hi+.35*(hi-lo))
    seen = set()
    for f, species, label, colour in sorted(lines):
        ax.axvline(f, color=colour, lw=1.1, ls='--', label=None if species in seen else species)
        seen.add(species)
        ax.text(f, ax.get_ylim()[1], label, rotation=90, fontsize=6, color=colour, ha='center', va='top',
                bbox=dict(facecolor='white', edgecolor='none', alpha=.75, pad=.5))
    ax.set_xlim(*(xlim or (freq.min(), freq.max())))
    ax.set_xlabel('Observed frequency (GHz)'); ax.set_ylabel('Jy/beam'); ax.set_title(title, fontsize=9)
    if z is not None:
        sec = ax.secondary_xaxis('top', functions=(lambda x: x*(1+z), lambda x: x/(1+z)))
        sec.set_xlabel(f'Rest frequency (GHz), z = {z:.5f}', fontsize=8)
    ax.legend(fontsize=7, loc='center left', bbox_to_anchor=(1.01, .5), frameon=False)
    return _save(fig, path)


def blend(path, fit, regions, names, unit):
    n = len(regions)
    fig, axes = plt.subplots(2, n, figsize=(6.5*n, 5), gridspec_kw=dict(height_ratios=[3, 1]), squeeze=False)
    for i, reg in enumerate(regions):
        ax, axr = axes[0, i], axes[1, i]
        m = fit['models'][i]; f = reg['freq']
        ax.step(f, reg['spectrum'], where='mid', color=TEAL, lw=.9, label='aperture spectrum')
        ax.plot(f, m['total'], color='k', lw=1.4, label='total fit')
        base = np.broadcast_to(np.asarray(m['offset'], float), np.shape(f))
        for j, (k, v) in enumerate(m['components'].items()):      # each line on its own (baseline not added)
            ax.plot(f, np.array(v), '--', color=COLOURS[j % len(COLOURS)], lw=1.1, label=names.get(k, k))
        if np.any(base != 0):
            ax.plot(f, base, ':', color=GREY, lw=1, label='baseline')
        ax.axhline(0, color='#94a3b8', lw=.5)
        ax.set_title(reg['id']); ax.set_ylabel(unit); ax.legend(fontsize=7)
        axr.step(f, m['residual'], where='mid', color=GREY, lw=.8)
        axr.fill_between(f, -3*np.asarray(reg['noise']), 3*np.asarray(reg['noise']), color='#cbd5e1', alpha=.6, step='mid')
        axr.axhline(0, color='#94a3b8', lw=.5); axr.set_xlabel('Observed frequency (GHz)'); axr.set_ylabel('residual')
    return _save(fig, path)


def robustness(path, tallies, title):
    labels = list(tallies)
    fig, ax = plt.subplots(figsize=(7, .5+.35*len(labels)))
    ax.barh(range(len(labels)), [tallies[k] for k in labels], color=TEAL)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(0, 1); ax.set_xlabel('Fraction of noise realisations'); ax.set_title(title)
    for i, k in enumerate(labels):
        ax.text(tallies[k], i, f' {tallies[k]:.2f}', va='center', fontsize=7)
    return _save(fig, path)


def maps(path, images, titles, centre=None, zoom=3, unit='Jy/beam km/s', units=None):
    """Moment / blend-fit maps; panels whose unit is km/s are drawn as velocity fields."""
    n = len(images)
    units = list(units) if units else [unit]*n
    cols = min(4, n); rows = int(np.ceil(n/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6*cols, 4.0*rows), squeeze=False, layout='constrained')
    for ax, image, title, u in zip(axes.flat, images, titles, units):
        velocity = u == 'km/s'
        marks = [(centre[0], centre[1], 'k+')] if centre is not None else []
        _image(ax, image, title, centre, zoom, cmap='coolwarm' if velocity else 'RdBu_r', symmetric=True, label=u or unit, marks=marks)
    for ax in list(axes.flat)[n:]:
        ax.axis('off')
    return _save(fig, path)


def sed(path, photometry, title):
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    pts = [(p['wavelength_um'], p['flux_Jy']) for p in photometry if p.get('wavelength_um') and p.get('flux_Jy') and p['flux_Jy'] > 0]
    if pts:
        w, fl = zip(*pts)
        ax.loglog(w, fl, 'o', color=TEAL, ms=4)
    ax.set_xlabel('Wavelength (µm)'); ax.set_ylabel('Flux density (Jy)'); ax.set_title(title)
    return _save(fig, path)
