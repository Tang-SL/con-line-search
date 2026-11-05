"""Step: moment maps of the adopted lines, and per-line maps of blended features.

For every line in solution.json: a signed moment-0 map over the feature's channels (baseline
fitted per pixel beside the line, nothing clipped), an intensity-weighted velocity map
(moment 1, pixels with a 3σ moment-0 only) and, for features with a blend fit, one map per
line from the per-pixel blend fit (conlines.blends.component_maps), plus its velocity field.
"""
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from ..core import read_cube, moment, C
from .. import blends
from . import figs


def _header(cube, unit, history=()):
    h = WCS(cube.header).celestial.to_header()
    for k in ('BMAJ', 'BMIN', 'BPA', 'OBJECT', 'SPECSYS'):
        if k in cube.header:
            h[k] = cube.header[k]
    h['BUNIT'] = unit
    for line in history:
        h.add_history(str(line).encode('ascii', 'replace').decode()[:70])
    return h


def moment1(cube, lo, hi, rest, z, mom0):
    f = cube.frequency[lo:hi+1]; e = rest/(1+z)
    v = C*(e-f)/e
    data = cube.data[lo:hi+1]
    noise = np.nanmedian(np.abs(data-np.nanmedian(data, axis=(1, 2), keepdims=True)), axis=(1, 2))*1.4826
    dv = abs(v[1]-v[0])
    err0 = np.sqrt(np.sum(noise**2))*dv
    w = np.clip(data, 0, None) if np.nansum(mom0) >= 0 else np.clip(-data, 0, None)
    with np.errstate(invalid='ignore', divide='ignore'):
        m1 = np.nansum(w*v[:, None, None], axis=0)/np.nansum(w, axis=0)
    m1[~(np.abs(mom0) > 3*err0)] = np.nan
    return m1, err0


def run(run, baseline='linear', zoom=4., local=True):
    sol = run.load('solution'); z = sol['z']
    feats = {f['id']: f for f in run.load('features')['features']}
    pos = run.state.get('position')
    cubes = {}
    def cube_of(key):
        if key not in cubes:
            cubes[key] = read_cube(run.cubes()[key]['path'])
        return cubes[key]
    products, panels = [], []
    done_fits = set()
    for e in sol['lines']:
        if e['status'] == 'unidentified' and not e['lines']:
            continue
        f = feats[e['feature']]; cube = cube_of(f['spw'])
        main = e['lines'][0] if e['lines'] else None
        rest = main['rest_GHz'] if main else f['freq_GHz']*(1+z)
        im, h = moment(cube, f['lo'], f['hi'], rest, z, baseline)
        name = (main['species'] if main else 'unidentified').replace(' ', '').replace(',', '_')[:30]
        h['LINEID'] = (main['species'] if main else 'unidentified').encode('ascii', 'replace').decode()[:60]
        h['IDSTAT'] = e['status'].upper()
        if len(e['lines']) > 1:
            h.add_history('Blended feature: this map integrates every line in it.')
        path = run.product(f"{e['feature']}__{name}__mom0.fits"); fits.writeto(path, im.astype('float32'), h, overwrite=True); products.append(path.name)
        m1, err0 = moment1(cube, f['lo'], f['hi'], rest, z, im)
        path = run.product(f"{e['feature']}__{name}__mom1.fits")
        fits.writeto(path, m1.astype('float32'), _header(cube, 'km/s', ['Intensity-weighted velocity relative to the redshift', 'Pixels with |mom0| > 3 sigma only']), overwrite=True)
        products.append(path.name)
        x, y = (float(v) for v in WCS(cube.header).celestial.world_to_pixel_values(pos['ra'], pos['dec'])) if pos else (None, None)
        panels.append(dict(images=[im, m1], titles=[f"{e['feature']} {h['LINEID']} [{e['status']}]: moment 0", 'moment 1 (km/s)'],
                           centre=(x, y), units=[h['BUNIT'], 'km/s']))
        # per-line maps from a blend fit
        if e.get('blend_fit') and e['blend_fit'] not in done_fits:
            slug = ''.join(ch if ch.isalnum() else '_' for ch in e['blend_fit'])[:60]
            fit = run.load(f'blend_{slug}')
            done_fits.add(e['blend_fit'])
            cubes_r = {i: cube_of(r['spw']) for i, r in enumerate(fit['regions'])}
            noises = {}
            for i, r in enumerate(fit['regions']):
                d = np.load(run.path/'cache'/f"spectrum_{r['spw']}.npz"); noises[i] = d['noise']
            maps, diag = blends.component_maps(fit, cubes_r, noises, local=local)
            images, titles = [], []
            for key, (img, i) in maps.items():
                cube_i = cubes_r[i]; nm = fit['names'].get(key, key)
                hh = _header(cube_i, cube_i.unit+' km/s', [f"Per-pixel blend fit {e['blend_fit']}: {nm}", diag[i]['note']])
                hh['LINEID'] = nm.encode('ascii', 'replace').decode()[:60]
                safe = ''.join(ch if ch.isalnum() else '_' for ch in nm)[:40]
                p = run.product(f"{fit['regions'][i]['feature']}__{safe}__blendfit.fits"); fits.writeto(p, img.astype('float32'), hh, overwrite=True); products.append(p.name)
                images.append(img); titles.append(f'{nm} (blend fit)')
            for i, d in diag.items():
                if np.isfinite(d['velocity']).any():
                    p = run.product(f"{fit['regions'][i]['feature']}__blendfit_velocity.fits")
                    fits.writeto(p, d['velocity'].astype('float32'), _header(cubes_r[i], 'km/s', [d['note']]), overwrite=True); products.append(p.name)
                    images.append(d['velocity']); titles.append(f"{fit['regions'][i]['feature']}: per-pixel velocity (km/s)")
            panels.append(dict(images=images, titles=titles, centre=(x, y), units=[None]*len(images), blend=e['blend_fit'],
                               warnings=[f"{fit['regions'][i]['feature']}: velocity field reaches the edge of the search window" for i, d in diag.items() if d.get('limited')]))
    plots = []
    for n, p in enumerate(panels):
        path = run.plot(f'maps_{n+1}')
        figs.maps(path, p['images'], p['titles'], p['centre'] if p['centre'][0] is not None else None, zoom)
        plots.append(dict(plot=path.name, blend=p.get('blend'), warnings=p.get('warnings', [])))
    result = dict(products=products, plots=plots, baseline=baseline, local=local)
    run.save('maps', result)
    run.log('maps', f'{len(products)} FITS products', baseline=baseline)
    return result


def summary(result):
    return '\n'.join([f"{len(result['products'])} FITS maps in products/:"]+['  '+p for p in result['products']]
                     + ['  ! '+w for p in result['plots'] for w in p['warnings']])
