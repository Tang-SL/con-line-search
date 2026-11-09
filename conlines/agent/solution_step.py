"""Step: record the adopted line solution and derive line properties.

The agent writes its decision as JSON (one entry per feature):

    {"lines": [
      {"feature": "spw27-F006", "status": "secure",          # secure | tentative | unidentified
       "main":  {"species": "HCN, v = 0", "rest_GHz": 265.886434, "qns": "J=3-2"},
       "blend": [{"species": "HCN, v2 = 1", "rest_GHz": 265.852709, "qns": "J=3-2 l=1e",
                  "status": "tentative"}],                   # optional per-line status (default: the feature's)
       "blend_fit": "HCN+vib, e=f",                          # optional: name of a blend fit
       "evidence": "why this is right (significance, companions, literature)",
       "alternatives": [{"solution": "CH2NH ...", "why_not": "..."}]}]}

Derived for every line: observed and rest frequency, velocity relative to the redshift, width,
integrated flux (from the blend fit when given), aperture flux (Jy km/s), and line luminosities
L' (K km/s pc^2) and L (L_sun) from Solomon & Vanden Bout (2005) when the distance is known.
"""
import json
from pathlib import Path
import numpy as np
from ..core import C

STATUSES = ('secure', 'tentative', 'unidentified')


def _aperture_factor(run):
    """Jy/beam (aperture mean) x factor = Jy summed over the aperture."""
    feats = run.load('features'); inv = {r['key']: r for r in run.load('inventory')['cubes']}
    out = {}
    for key, s in feats['spectra'].items():
        r = inv[key]
        beam_area = 1.1331*r['beam_major_arcsec']*r['beam_minor_arcsec']
        out[key] = s['pixels']*r['pixel_arcsec']**2/beam_area
    return out


def luminosities(flux_Jykms, nu_rest_GHz, z, dl_Mpc):
    if not dl_Mpc or flux_Jykms is None:
        return None, None
    nu_obs = nu_rest_GHz/(1+z)
    lprime = 3.25e7*flux_Jykms*nu_obs**-2*dl_Mpc**2*(1+z)**-3
    lsun = 1.04e-3*flux_Jykms*nu_rest_GHz*dl_Mpc**2
    return float(lprime), float(lsun)


def run(run, solution):
    if isinstance(solution, (str, Path)):
        p = Path(solution)
        solution = json.loads(p.read_text() if p.exists() else str(solution))
    feats = {f['id']: f for f in run.load('features')['features']}
    z = run.state['z']
    ctx = run.load('context', {})
    dl = (ctx.get('derived') or {}).get('luminosity_distance_Mpc')
    factor = _aperture_factor(run)
    blends = {}
    for p in run.path.glob('blend_*.json'):
        b = json.loads(p.read_text()); blends[b['name']] = b
    out, problems = [], []
    for item in solution['lines']:
        fid = item['feature']
        if fid not in feats:
            problems.append(f'{fid}: no such feature'); continue
        if item.get('status') not in STATUSES:
            problems.append(f"{fid}: status must be one of {STATUSES}"); continue
        f = feats[fid]
        lines = ([item['main']] if item.get('main') else [])+list(item.get('blend', []))
        fit = blends.get(item.get('blend_fit')) if item.get('blend_fit') else None
        if item.get('blend_fit') and fit is None:
            problems.append(f"{fid}: blend fit {item['blend_fit']!r} not found")
        rows = []
        for i, line in enumerate(lines):
            rest = float(line['rest_GHz']); expected = rest/(1+z)
            v = C*(expected-f['freq_GHz'])/expected
            flux, err, share, verdict = (f['flux_Jykms'], f['flux_error_Jykms'], 1., None) if len(lines) == 1 else (None, None, None, None)
            tied_to = None
            if fit:
                q = next((q for q in fit['components'] if abs(q['row']['rest_GHz']-rest) < 5e-4
                          and fit['regions'][q['region']]['feature'] == fid), None)
                if q:
                    flux, err, share, verdict = q['flux'], q['flux_error'], q.get('fraction'), q['verdict']
                    if q.get('tie'):        # not measured: set from another line through the hypothesis' tie
                        key, ratio = q['tie']
                        tied_to = dict(line=fit['names'].get(key, key), feature=key.split('|')[0], ratio=ratio, key=key)
            line_snr = abs(flux)/err if flux is not None and err and not tied_to else None
            ap = flux*factor[f['spw']] if flux is not None else None
            lp, ls = luminosities(abs(ap) if ap is not None else None, rest, z, dl)
            if line.get('status', item['status']) not in STATUSES:
                problems.append(f"{fid} {line['species']}: status must be one of {STATUSES}")
            rows.append(dict(role='main' if i == 0 and item.get('main') else 'blended', species=line['species'], qns=line.get('qns', ''),
                             status=line.get('status', item['status']),
                             rest_GHz=rest, expected_GHz=expected, velocity_kms=float(v) if len(lines) == 1 else None,
                             flux_Jybkms=flux, flux_error_Jybkms=err, snr=line_snr, tied_to=tied_to, share=share, fit_verdict=verdict,
                             aperture_flux_Jykms=ap, Lprime_K_kms_pc2=lp, L_Lsun=ls))
        out.append(dict(feature=fid, spw=f['spw'], sign=f['sign'], status=item['status'], observed_GHz=f['freq_GHz'],
                        width_kms=f['width_kms'], sigma_kms=f.get('sigma_kms'), snr=f['snr'], reliability=f['reliability'],
                        blend_fit=item.get('blend_fit'), blend_fit_reason=item.get('blend_fit_reason', ''), lines=rows, evidence=item.get('evidence', ''),
                        alternatives=item.get('alternatives', [])))
    # name tie partners as the solution names them (the catalogue's quantum-number codes are cryptic)
    named = {(e['feature'], round(r['rest_GHz'], 4)): f"{r['species']} {r['qns']}".strip() for e in out for r in e['lines']}
    for e in out:
        for r in e['lines']:
            t = r.get('tied_to')
            if t and t.get('key'):
                _, _, rest = t['key'].split('|')
                t['line'] = named.get((t['feature'], round(float(rest), 4)), t['line'])
    result = dict(z=z, luminosity_distance_Mpc=dl, lines=out, problems=problems)
    run.save('solution', result)
    _overview(run, result, z)
    # flat table for spreadsheets
    import csv
    with (run.path/'products'/'line_solution.csv').open('w', newline='') as stream:
        w = csv.writer(stream)
        w.writerow(['feature', 'status', 'role', 'species', 'qns', 'rest_GHz', 'observed_GHz', 'velocity_kms', 'width_kms',
                    'flux_Jybeam_kms', 'flux_error', 'tied_to', 'share', 'aperture_flux_Jy_kms', 'Lprime_K_kms_pc2', 'L_Lsun',
                    'line_S/N', 'feature_S/N', 'reliability'])
        for e in out:
            for r in e['lines'] or [dict(status=e['status'], role='', species='unidentified', qns='', rest_GHz=None, velocity_kms=None, flux_Jybkms=None, flux_error_Jybkms=None,
                                        share=None, aperture_flux_Jykms=None, Lprime_K_kms_pc2=None, L_Lsun=None)]:
                w.writerow([e['feature'], r['status'], r['role'], r['species'], r['qns'], r['rest_GHz'], r.get('expected_GHz', e['observed_GHz']), r['velocity_kms'],
                            e['width_kms'], r['flux_Jybkms'], r['flux_error_Jybkms'],
                            (f"{r['tied_to']['line']} x {r['tied_to']['ratio']:g}" if r.get('tied_to') else ''), r['share'], r['aperture_flux_Jykms'],
                            r['Lprime_K_kms_pc2'], r['L_Lsun'], r.get('snr'), e['snr'], e['reliability']])
    run.log('solve', f'{len(out)} features in the solution', problems=problems)
    return result


def summary(result):
    lines = []
    for e in result['lines']:
        names = ' + '.join(r['species'] for r in e['lines']) or 'unidentified'
        lines.append(f"{e['feature']} [{e['status']}] {names}")
        for r in e['lines']:
            flux = 'n/a' if r['flux_Jybkms'] is None else f"{r['flux_Jybkms']:.3g}"
            error = '' if r['flux_error_Jybkms'] is None else f" ± {r['flux_error_Jybkms']:.2g}"
            if r.get('tied_to'):
                error += f" (not measured: tied to {r['tied_to']['line']})"
            lum = f", L' {r['Lprime_K_kms_pc2']:.3g} K km/s pc2" if r['Lprime_K_kms_pc2'] else ''
            tag = f" [{r['status']}]" if r['status'] != e['status'] else ''
            lines.append(f"    {r['role']:<8} {r['species']+tag:<28} {r['qns']:<16} flux {flux}{error} Jy/b km/s{lum}")
    lines += ['  ! '+p for p in result['problems']]
    return '\n'.join(lines)


def _overview(run, result, z):
    """plots/overview.png: every window's aperture spectrum with the adopted solution."""
    import numpy as np
    from . import figs
    windows = []
    for p in sorted((run.path/'cache').glob('spectrum_*.npz')):
        d = np.load(p)
        windows.append(dict(key=p.stem.replace('spectrum_', ''), frequency=d['frequency'], spectrum=d['spectrum'], noise=d['noise']))
    if windows:
        windows.sort(key=lambda w: float(np.nanmin(w['frequency'])))
        target = run.state.get('target') or ''
        figs.overview(run.plot('overview'), windows, result, z, f'{target}: adopted line solution'.strip(': '))
