"""Steps: blend decomposition of features, and robustness of competing solutions under the noise.

A *hypothesis* (JSON) names the features to fit together, the lines in them, the shared shape
and any ties:

    {"name": "HCN+vib, e=f", "shape": "gauss2",
     "features": ["spw27-F006", "spw25-F003"],
     "lines": [{"feature": "spw27-F006", "species": "HCN, v = 0", "rest_GHz": 265.8864},
               {"feature": "spw27-F006", "species": "HCN, v2 = 1", "rest_GHz": 265.8527, "tie": 2, "ratio": 1},
               {"feature": "spw25-F003", "species": "HCN, v2 = 1", "rest_GHz": 267.1993},
               {"feature": "spw25-F003", "species": "HCO+, v = 0", "rest_GHz": 267.5576}]}

("tie": index of the line, 0-based, whose strength this one follows, times "ratio").
Robustness: every hypothesis is refitted to the observed spectrum plus many independent noise
realisations taken from the spatial control apertures (real noise, with its channel
correlations). The fraction of realisations in which each hypothesis has the lowest BIC, and
the spread of every line's flux, show which solution survives the noise.
"""
import json
import numpy as np
from .. import blends
from ..core import C
from . import figs


def spectra(run):
    out = {}
    for p in (run.path/'cache').glob('spectrum_*.npz'):
        d = np.load(p)
        out[p.stem[len('spectrum_'):]] = dict(frequency=d['frequency'], spectrum=d['spectrum'], noise=d['noise'], controls=d['controls'])
    return out


def _catalogue_row(run, line):
    """Attach catalogue metadata (E_up, A, S mu^2) to a line given by species and rest frequency."""
    try:
        cat = run.load('catalogue')
    except FileNotFoundError:
        return dict(species=line['species'], rest_GHz=float(line['rest_GHz']), qns=line.get('qns', ''))
    best = None
    for w in cat['windows'].values():
        for r in w['rows']:
            if abs(r['rest_GHz']-float(line['rest_GHz'])) < 5e-4 and (line['species'].split(',')[0].strip() in r['species']):
                if best is None or abs(r['rest_GHz']-line['rest_GHz']) < abs(best['rest_GHz']-line['rest_GHz']):
                    best = r
    row = dict(best) if best else {}
    row.update(species=line['species'], rest_GHz=float(line['rest_GHz']))
    row.setdefault('qns', line.get('qns', ''))
    return row


def build(run, hyp, z=None, spec=None, margin=10):
    feats = {f['id']: f for f in run.load('features')['features']}
    z = run.state['z'] if z is None else z
    spec = spec or spectra(run)
    regions = []
    for fid in hyp['features']:
        f = feats[fid]; sp = spec[f['spw']]; n = len(sp['frequency'])
        m = int(hyp.get('margin_channels', margin))
        lo, hi = max(0, f['lo']-m), min(n-1, f['hi']+m)
        regions.append(dict(id=fid, spw=f['spw'], lo=lo, hi=hi, freq=np.asarray(sp['frequency'][lo:hi+1], float),
                            spectrum=np.asarray(sp['spectrum'][lo:hi+1], float), noise=np.asarray(sp['noise'][lo:hi+1], float)))
    comps = []
    for i, line in enumerate(hyp['lines']):
        row = _catalogue_row(run, line)
        region = hyp['features'].index(line['feature'])
        key = f"{line['feature']}|{row['species']}|{row['rest_GHz']:.6f}"
        comps.append(dict(key=key, row=row, expected_GHz=row['rest_GHz']/(1+z), region=region, sign=feats[line['feature']]['sign'], tie=None))
    for i, line in enumerate(hyp['lines']):
        if line.get('tie') is not None:
            comps[i]['tie'] = (comps[int(line['tie'])]['key'], float(line.get('ratio', 1.)))
    template = None
    if hyp.get('shape') == 'template':
        t = hyp['template']; f = feats[t['feature']]; sp = spec[f['spw']]
        template = blends.make_template(sp['frequency'], sp['spectrum'], float(t['rest_GHz']), z, f['lo'], f['hi'])
    return regions, comps, template, z


def fit(run, hyp, save=True):
    regions, comps, template, z = build(run, hyp)
    result = blends.decompose(regions, comps, hyp.get('shape', 'gauss'), z, template, float(hyp.get('vmax', 300.)), hyp.get('baseline', True))
    names = {c['key']: f"{c['row']['species']} {c['row'].get('qns', '')}".strip()[:40] for c in comps}
    for r, reg in zip(result['regions'], regions):
        r['feature'] = reg['id']
    result.update(name=hyp['name'], hypothesis=hyp, names=names, shape_text=blends.describe_shape(result))
    if save:
        slug = ''.join(ch if ch.isalnum() else '_' for ch in hyp['name'])[:60]
        run.save(f'blend_{slug}', result)
        unit = 'Jy/beam'
        figs.blend(run.plot(f'blend_{slug}'), result, regions, names, unit)
        run.log('blend', f"{hyp['name']}: chi2_r {result['reduced_chi2']:.2f}", hypothesis=hyp)
    return result


def robustness(run, hyps, n=100, seed=1, name='robustness'):
    """Refit every hypothesis to observed + control-noise realisations; tally BIC winners."""
    rng = np.random.default_rng(seed)
    spec = spectra(run)
    base = [fit(run, h, save=False) for h in hyps]
    built = [build(run, h, spec=spec) for h in hyps]
    wins = np.zeros(len(hyps)); fluxes = [{q['key']: [] for q in b['components']} for b in base]
    for i in range(n):
        noisy = {}
        for key, sp in spec.items():
            ctrl = sp['controls']-np.median(sp['controls'], axis=0)
            noisy[key] = dict(sp, spectrum=sp['spectrum']+ctrl[rng.integers(len(ctrl))])
        bics = []
        for h, b, (regions, comps, template, z) in zip(hyps, base, built):
            regs = [dict(r, spectrum=np.asarray(noisy[r['spw']]['spectrum'][r['lo']:r['hi']+1], float)) for r in regions]
            try:
                res = blends.decompose(regs, comps, h.get('shape', 'gauss'), z, template, float(h.get('vmax', 300.)),
                                       h.get('baseline', True), assess=False, starts_override=[b['params']])
            except Exception:
                bics.append(np.inf); continue
            bics.append(res['bic'])
            for q in res['components']:
                fluxes[hyps.index(h)][q['key']].append(q['flux'])
        wins[int(np.argmin(bics))] += 1
    tallies = {h['name']: float(w/n) for h, w in zip(hyps, wins)}
    lines = []
    for h, b, fl in zip(hyps, base, fluxes):
        for q in b['components']:
            arr = np.array(fl[q['key']])
            lines.append(dict(hypothesis=h['name'], line=b['names'][q['key']], flux=q['flux'], flux_error=q['flux_error'],
                              median=float(np.median(arr)) if arr.size else None,
                              p16=float(np.percentile(arr, 16)) if arr.size else None, p84=float(np.percentile(arr, 84)) if arr.size else None,
                              fraction_positive=float(np.mean(arr*np.sign(q['flux'] or 1) > 0)) if arr.size else None))
    result = dict(name=name, n=n, note='observed spectrum + one control-aperture noise realisation each (conservative: doubles the noise variance)',
                  hypotheses=[dict(name=h['name'], bic=b['bic'], reduced_chi2=b['reduced_chi2'], n_parameters=b['n_parameters']) for h, b in zip(hyps, base)],
                  win_fraction=tallies, lines=lines)
    slug = ''.join(ch if ch.isalnum() else '_' for ch in name)[:60]
    run.save(f'robust_{slug}', result)
    figs.robustness(run.plot(f'robust_{slug}'), tallies, f'{name}: lowest BIC in {n} noise realisations')
    run.log('robust', f'{name}: {tallies}', n=n)
    return result


def load_hypotheses(path_or_json):
    text = path_or_json
    try:
        from pathlib import Path
        if Path(path_or_json).exists():
            text = Path(path_or_json).read_text()
    except OSError:
        pass
    data = json.loads(text)
    return data if isinstance(data, list) else [data]


def summary_fit(r):
    lines = [f"{r['name']}: shape {r['shape_text']}; reduced chi2 {r['reduced_chi2']:.2f}, BIC {r['bic']:.1f}" + (f"; AT BOUNDS: {r['at_bounds']}" if r['at_bounds'] else '')]
    for q in r['components']:
        lines.append(f"    {r['names'][q['key']]:<40} flux {q['flux']:.3g} ± {q['flux_error']:.2g} Jy/b km/s  share {q.get('fraction', 0):.2f}  "
                     f"dBIC {q['delta_bic'] if q['delta_bic'] is None else round(q['delta_bic'], 1)}  {q['verdict']}"
                     + (f"  trades flux with {[r['names'][d] for d in q['degenerate_with']]}" if q['degenerate_with'] else ''))
    return '\n'.join(lines)


def summary_robust(r):
    lines = [f"{r['name']} ({r['n']} realisations): " + ', '.join(f'{k}: {v:.2f}' for k, v in r['win_fraction'].items())]
    for q in r['lines']:
        lines.append(f"    [{q['hypothesis']}] {q['line']:<40} flux {q['flux']:.3g}; noise realisations: median {q['median']:.3g} (16-84%: {q['p16']:.3g} to {q['p84']:.3g}), same sign in {q['fraction_positive']:.0%}")
    return '\n'.join(lines)
