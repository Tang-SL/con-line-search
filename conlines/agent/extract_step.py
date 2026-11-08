"""Step: aperture spectra at the source, feature significance and detection reliability.

Significance (per feature, in the aperture): integrated S/N under two baselines, channel
support, baseline sensitivity, edges (espada_backend.aperture_features). Reliability (per
feature, catalogue-free): how often equally strong and wide components occur away from the
source (reliability.detection_reliability), i.e. an empirical artefact rate, in the spirit of
the fidelity of interferometric line searches (e.g. ASPECS: 1 - N_neg/N_pos).
"""
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from ..core import read_cube, extract, C
from ..espada_backend import aperture_features
from ..measurements import matching_frequency
from ..reliability import channel_noise, component_properties, detection_reliability, neighbour_counts
from . import figs


def choose_position(run, how='auto', ra=None, dec=None):
    inv = run.load('inventory')['cubes']
    if how == 'radec' or (ra is not None and dec is not None):
        return SkyCoord(ra, dec, unit='deg'), f'given RA/Dec {ra:.6f}, {dec:.6f}'
    if how in ('auto', 'line'):
        best = max(inv, key=lambda r: r['line_peak_snr'])
        if how == 'line' or best['line_peak_snr'] > 8:
            return SkyCoord(best['line_peak_ra'], best['line_peak_dec'], unit='deg'), f"brightest line emission ({best['key']}, S/N {best['line_peak_snr']:.0f})"
    if how in ('auto', 'continuum'):
        best = max(inv, key=lambda r: r['continuum_snr'] or 0)
        if (best['continuum_snr'] or 0) > 5:
            cube = read_cube(best['path'])
            x, y = best['continuum_pixel']
            w = WCS(cube.header).celestial.pixel_to_world(x, y)
            return w, f"continuum peak ({best['key']}, S/N {best['continuum_snr']:.0f})"
    r = inv[0]
    return SkyCoord(r['centre_ra'], r['centre_dec'], unit='deg'), 'image centre'


def run(run, how='auto', ra=None, dec=None, radius_arcsec=None, threshold=8., region_beams=5., min_reliability=.5):
    position, reason = choose_position(run, how, ra, dec)
    comps = run.load('components')
    inv = {r['key']: r for r in run.load('inventory')['cubes']}
    radius = radius_arcsec or max(r['beam_major_arcsec'] for r in inv.values())
    region = region_beams*max(r['beam_major_arcsec'] for r in inv.values())
    features, spectra, scores, positions, ids = [], {}, {}, [], []
    for key, record in run.cubes().items():
        cube = read_cube(record['path'])
        w = WCS(cube.header).celestial
        x, y = (float(v) for v in w.world_to_pixel(position))
        ex = extract(cube, x, y, 'Circle', radius)
        np.savez_compressed(run.path/'cache'/f'spectrum_{key}.npz', frequency=cube.frequency, spectrum=ex['spectrum'],
                            noise=ex['noise'], controls=ex['controls'], mask=ex['mask'])
        spectra[key] = dict(x=x, y=y, pixels=ex['pixels'])
        run_k = comps.get(key, {'components': []})
        found = aperture_features(cube, ex, run_k, threshold)
        # Reliability of every component of this cube (artefact rate away from the source).
        cs = run_k['components']
        if cs:
            props = component_properties(cube, cs, channel_noise(cube))
            world = w.pixel_to_world([c['x'] for c in cs], [c['y'] for c in cs])
            dra, ddec = position.spherical_offsets_to(world)
            offsets = np.c_[dra.to_value('arcsec'), ddec.to_value('arcsec')]
            dist = np.hypot(offsets[:, 0], offsets[:, 1])
            yy, xx = np.indices(cube.data.shape[1:]); sx, sy = cube.pixel_arcsec
            rr = np.hypot((xx-x)*sx, (yy-y)*sy); data_ok = np.isfinite(cube.data[len(cube.frequency)//2])
            areas = dict(inside=int(((rr <= region) & data_ok).sum()), near=int(((rr > region) & (rr <= 3*region) & data_ok).sum()),
                         far=int(((rr > 3*region) & data_ok).sum()))
            for c, p, d, rel, xy in zip(cs, props, dist, detection_reliability(props, dist, region, areas), offsets):
                scores[c['id']] = dict(p, **rel, distance_arcsec=float(d))
                positions.append(xy); ids.append(c['id'])
        df = abs(cube.frequency[1]-cube.frequency[0])
        for q in found:
            m = q.get('aperture_measurement', {})
            s = scores.get(q['id'], {})
            fq = matching_frequency(q)
            lo_f, hi_f = sorted([cube.frequency[q['lo']], cube.frequency[q['hi']]])
            sign = 1 if q['sign'] == 'Emission' else -1
            v = C*df/fq
            flux = float(sign*np.sum(ex['spectrum'][q['lo']:q['hi']+1])*v)          # Jy/beam km/s (signed: >0 means the stated sign)
            flux_err = float(np.sqrt(np.sum(ex['noise'][q['lo']:q['hi']+1]**2))*v)
            features.append(dict(id=q['id'], spw=key, sign=q['sign'], lo=q['lo'], hi=q['hi'], channels=q['hi']-q['lo']+1,
                                 freq_GHz=fq, freq_lo_GHz=float(lo_f-df/2), freq_hi_GHz=float(hi_f+df/2),
                                 width_kms=float(C*(hi_f-lo_f+df)/fq), snr=q['snr'], significant=q['significant'], flags=q['flags'],
                                 centroid_error_kms=m.get('centroid_error_kms') if m.get('reliable') else None,
                                 sigma_kms=m.get('sigma_kms') if m.get('reliable') else None,
                                 flux_Jykms=flux, flux_error_Jykms=flux_err,
                                 reliability=s.get('reliability'), reliability_best=s.get('reliability_best'),
                                 expected_false=s.get('expected_false'), peak_snr=s.get('peak_snr'),
                                 distance_arcsec=s.get('distance_arcsec')))
        figs.spectrum(run.plot(f'spectrum_{key}'), cube.frequency, ex['spectrum'], ex['noise'],
                      [dict(f, reliable=(f['reliability'] or 0) >= min_reliability) for f in features if f['spw'] == key],
                      title=f'{key}: aperture spectrum (r = {radius:.2f}″) at {reason}', unit=cube.unit,
                      z=run.state.get('z'))
        del cube
    for ident, n in zip(ids, neighbour_counts(positions, max(r['beam_major_arcsec'] for r in inv.values()))):
        scores[ident]['neighbours'] = n
    for f in features:
        f['neighbours'] = scores.get(f['id'], {}).get('neighbours')
        f['reliable'] = f['reliability'] is not None and f['reliability'] >= min_reliability
        if f['reliability'] is None and 'reliability untestable' not in f['flags']:
            f['flags'] = list(f['flags'])+['reliability untestable']      # no control zone as large as the source region
        f['verdict'] = ('significant and reliable' if f['significant'] and f['reliable'] else
                        'reliable but not significant in the aperture' if f['reliable'] else
                        'significant; reliability untestable (field too small for controls)' if f['significant'] and f['reliability'] is None else
                        'significant in the aperture but possibly an artefact' if f['significant'] else
                        'weak; reliability untestable' if f['reliability'] is None else 'weak / likely artefact')
    # Trials (look-elsewhere): independent spatial-spectral elements searched.
    trials = {}
    for key, r in inv.items():
        beam_area = 1.1331*r['beam_major_arcsec']*r['beam_minor_arcsec']
        field = r['shape'][1]*r['shape'][2]*r['finite_fraction']*r['pixel_arcsec']**2
        trials[key] = dict(beams=field/beam_area, channels=r['nchan'], independent_elements=field/beam_area*r['nchan']/3)
    result = dict(position=dict(ra=position.ra.deg, dec=position.dec.deg, reason=reason), radius_arcsec=radius,
                  region_radius_arcsec=region, threshold=threshold, min_reliability=min_reliability,
                  spectra=spectra, features=features, control_scores=scores, trials=trials)
    run.save('features', result)
    run.update(position=result['position'], radius_arcsec=radius)
    run.log('extract', f'{len(features)} features at {reason}', radius=radius, threshold=threshold)
    return result


def summary(result):
    lines = [f"Aperture r = {result['radius_arcsec']:.3f}\" at {result['position']['reason']} "
             f"(RA {result['position']['ra']:.6f}, Dec {result['position']['dec']:.6f}); source region {result['region_radius_arcsec']:.2f}\""]
    for f in sorted(result['features'], key=lambda f: f['freq_GHz']):
        lines.append(f"  {f['id']}: {f['sign'][:3]} {f['freq_GHz']:.5f} GHz, {f['channels']} ch ({f['width_kms']:.0f} km/s), "
                     f"S/N {f['snr'] if f['snr'] is None else round(f['snr'], 1)}, reliability {f['reliability'] if f['reliability'] is None else round(f['reliability'], 2)}, "
                     f"flux {f['flux_Jykms']:.3g}±{f['flux_error_Jykms']:.2g} Jy/b km/s -> {f['verdict']}" + (f" [{', '.join(f['flags'])}]" if f['flags'] else ''))
    return '\n'.join(lines)
