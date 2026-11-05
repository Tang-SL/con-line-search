"""Static figures (PNG) for the weblog: every step gets plots that can be checked by eye."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ..core import C

TEAL, ORANGE, PURPLE, GREY, RED, BLUE = '#167f83', '#c2410c', '#7356ad', '#6b7280', '#b91c1c', '#1d4ed8'
COLOURS = [ORANGE, PURPLE, '#0f766e', '#b45309', BLUE, '#be185d', '#4d7c0f', GREY]
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'figure.dpi': 110, 'axes.spines.top': False, 'axes.spines.right': False})


def _save(fig, path):
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def _image(ax, image, title, centre=None, zoom=None, cmap='viridis', symmetric=False, label=None, marks=()):
    data = np.asarray(image, float)
    ny, nx = data.shape
    finite = data[np.isfinite(data)]
    if finite.size:
        if symmetric:
            v = np.nanpercentile(np.abs(finite), 99.5) or 1
            vmin, vmax = -v, v
        else:
            vmin, vmax = np.nanpercentile(finite, [1, 99.7])
    else:
        vmin, vmax = 0, 1
    h = ax.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    if centre is not None and zoom:
        cx, cy = centre; half = max(8, min(nx, ny)/(2*zoom))
        ax.set_xlim(cx-half, cx+half); ax.set_ylim(cy-half, cy+half)
    for x, y, style in marks:
        ax.plot(x, y, style, ms=9, mfc='none', mew=1.6)
    ax.set_title(title); ax.set_xlabel('x (pixel)'); ax.set_ylabel('y (pixel)')
    cb = plt.colorbar(h, ax=ax, fraction=.046, pad=.03)
    if label:
        cb.set_label(label)
    return h


def inspection(path, cube, info, median_map, peak_snr, z=None):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.6]))
    lx, ly = info['line_peak_pixel']; cx, cy = info['continuum_pixel']
    _image(ax[0], median_map, 'Continuum estimate (25th percentile over channels)', (cx, cy), 3, label=info['unit'], marks=[(cx, cy, 'w+')])
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
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
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


def companions(path, rows, title):
    """Predicted vs measured integrated intensity for every transition of one candidate species."""
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    x = np.arange(len(rows))
    pred = np.array([r['predicted'] for r in rows]); meas = np.array([r['measured'] for r in rows]); err = np.array([r['error'] for r in rows])
    ax.bar(x-.2, pred, .4, color=PURPLE, alpha=.7, label='LTE prediction (scaled to the feature)')
    ax.errorbar(x+.2, meas, yerr=3*err, fmt='o', color=TEAL, label='measured ±3σ')
    for i, r in enumerate(rows):
        if r.get('status') == 'missing':
            ax.text(i, max(pred[i], meas[i]), ' missing', color=RED, fontsize=7, rotation=90, va='bottom')
        elif str(r.get('status', '')).startswith('blended'):
            ax.text(i, max(pred[i], meas[i]), ' blended', color=GREY, fontsize=7, rotation=90, va='bottom')
    ax.set_xticks(x); ax.set_xticklabels([r['label'] for r in rows], rotation=60, ha='right', fontsize=6)
    ax.set_ylabel('Integrated intensity (Jy/beam km/s)'); ax.axhline(0, color='#94a3b8', lw=.5)
    ax.set_title(title); ax.legend(fontsize=7)
    return _save(fig, path)


def blend(path, fit, regions, names, unit):
    n = len(regions)
    fig, axes = plt.subplots(2, n, figsize=(6.5*n, 5), gridspec_kw=dict(height_ratios=[3, 1]), squeeze=False)
    for i, reg in enumerate(regions):
        ax, axr = axes[0, i], axes[1, i]
        m = fit['models'][i]; f = reg['freq']
        ax.step(f, reg['spectrum'], where='mid', color=TEAL, lw=.9, label='aperture spectrum')
        ax.plot(f, m['total'], color='k', lw=1.4, label='total fit')
        for j, (k, v) in enumerate(m['components'].items()):
            ax.plot(f, np.array(v)+m['offset'], '--', color=COLOURS[j % len(COLOURS)], lw=1.1, label=names.get(k, k))
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


def maps(path, images, titles, centre=None, zoom=3, unit='Jy/beam km/s'):
    n = len(images)
    cols = min(4, n); rows = int(np.ceil(n/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.4*cols, 3.8*rows), squeeze=False)
    for ax, image, title in zip(axes.flat, images, titles):
        _image(ax, image, title, centre, zoom, cmap='RdBu_r', symmetric=True, label=unit)
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
