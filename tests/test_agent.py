"""Line-identification agent toolkit: run folder, data discovery, and an offline end-to-end run on a
synthetic two-window cube (HCN 3-2 + HCO+ 3-2 + HCN-vib l=1f), from detection to the weblog."""
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from conlines.agent.__main__ import main
from conlines.agent.data import discover
from conlines.agent.rundir import Run
from conlines.agent.solution_step import luminosities
from conlines.core import C

Z = 0.0428
HCN, VIB_E, VIB_F, HCOP = 265.886434, 265.852709, 267.199283, 267.557626


def _sofia():
    exe = os.environ.get('SOFIA_EXECUTABLE')
    if exe and Path(exe).exists():
        return exe
    local = Path(__file__).resolve().parents[1]/'.local-tools'/'sofia'
    return str(local) if local.exists() else shutil.which('sofia')


def test_run_folder_round_trip_and_json_cleaning(tmp_path):
    run = Run(tmp_path/'r')
    run.save('x', dict(a=np.float32(1.5), b=np.array([1, 2]), c=float('nan'), d=Path('/tmp')))
    assert run.load('x') == dict(a=1.5, b=[1, 2], c=None, d='/tmp')
    assert run.load('missing', {'k': 1}) == {'k': 1}
    with pytest.raises(FileNotFoundError):
        run.load('missing')
    run.update(z=0.1, target='T')
    assert Run(tmp_path/'r').state['z'] == 0.1
    run.log('step', 'message', value=np.int64(3))
    entry = json.loads((run.path/'log.jsonl').read_text().splitlines()[-1])
    assert entry['step'] == 'step' and entry['value'] == 3
    for sub in ('plots', 'products', 'notes', 'cache'):
        assert (run.path/sub).is_dir()


def test_discover_groups_alma_products_by_spectral_window(tmp_path):
    base = 'member.uid___A001_X1_X2.NGC_4418_sci'
    for name in [f'{base}.spw25.cube.I.pbcor.fits', f'{base}.spw25.cube.I.pb.fits.gz', f'{base}.spw25.cube.I.mask.fits.gz',
                 f'{base}.spw27.cube.I.pbcor.fits', f'{base}.spw25_27.cont.I.tt0.pbcor.fits']:
        (tmp_path/name).write_bytes(b'')
    cubes, continuum, notes = discover(tmp_path)
    assert [c['key'] for c in cubes] == ['spw25', 'spw27']
    assert cubes[0]['pb'].endswith('pb.fits.gz') and cubes[0]['mask'].endswith('mask.fits.gz')
    assert cubes[1]['pb'] is None
    assert len(continuum) == 1


def test_luminosities_follow_solomon_and_vanden_bout():
    lp, l = luminosities(1., 115.271, 0.01, 43.)
    assert lp == pytest.approx(3.25e7*(115.271/1.01)**-2*43**2*1.01**-3)
    assert l == pytest.approx(1.04e-3*115.271*43**2)
    assert luminosities(1., 100., 0.01, None) == (None, None)


def _write_cube(path, f0_GHz, nchan, lines, rng, n=128, df=3.9e6):
    """Gaussian source (1.5 px = 0.15″) with double-peaked lines: lines = [(rest_GHz, peak Jy/beam)]."""
    freq = f0_GHz*1e9+df*np.arange(nchan)
    yy, xx = np.indices((n, n)); src = np.exp(-((xx-n//2)**2+(yy-n//2)**2)/(2*1.5**2))
    spec = np.zeros(nchan)
    for rest, peak in lines:
        e = rest*1e9/(1+Z); v = C*(e-freq)/e
        spec += peak*(np.exp(-(v-170)**2/(2*130**2))+.5*np.exp(-(v+200)**2/(2*130**2)))
    data = spec[:, None, None]*src[None]+rng.normal(0, 1.5e-3, (nchan, n, n))
    h = fits.Header()
    h.update(NAXIS=3, CTYPE1='RA---SIN', CRVAL1=260.8415, CDELT1=-0.1/3600, CRPIX1=n//2+1, CUNIT1='deg',
             CTYPE2='DEC--SIN', CRVAL2=-0.2836, CDELT2=0.1/3600, CRPIX2=n//2+1, CUNIT2='deg',
             CTYPE3='FREQ', CRVAL3=freq[0], CDELT3=df, CRPIX3=1, CUNIT3='Hz', SPECSYS='LSRK',
             BUNIT='Jy/beam', BMAJ=0.3/3600, BMIN=0.25/3600, BPA=50., OBJECT='SYNTH-CON', RESTFRQ=freq[nchan//2])
    fits.writeto(path, data.astype('float32'), h, overwrite=True)


@pytest.mark.skipif(_sofia() is None, reason='SoFiA 2 is not installed')
def test_agent_pipeline_end_to_end(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('SOFIA_EXECUTABLE', _sofia())
    rng = np.random.default_rng(3)
    data = tmp_path/'data'; data.mkdir()
    base = 'member.uid___A001_X1_X2.SYNTH-CON_sci'
    # spw27: HCN + HCN-vib l=1e; spw25: HCO+ + HCN-vib l=1f (20% of HCN), both 1.9 GHz wide, as ALMA spws
    _write_cube(data/f'{base}.spw27.cube.I.pbcor.fits', 254.05, 480, [(HCN, .03), (VIB_E, .006)], rng)
    _write_cube(data/f'{base}.spw25.cube.I.pbcor.fits', 255.65, 480, [(HCOP, .02), (VIB_F, .006)], rng)
    catalogue = tmp_path/'lines.json'
    catalogue.write_text(json.dumps([
        dict(species='HCN, v = 0', rest_GHz=HCN, qns='03-02', Eu_K=25.5, log10_A_s=-3.08, S_mu2=80.2, seen_in_space=True),
        dict(species='HCN, v2 = 1', rest_GHz=VIB_E, qns='0301-02-1', Eu_K=1050., log10_A_s=-3.15, S_mu2=72.4, seen_in_space=True),
        dict(species='HCN, v2 = 1', rest_GHz=VIB_F, qns='03-1-0201', Eu_K=1050., log10_A_s=-3.14, S_mu2=72.4, seen_in_space=True),
        dict(species='HCO+, v = 0', rest_GHz=HCOP, qns='03-02', Eu_K=25.7, log10_A_s=-2.84, S_mu2=45.6, seen_in_space=True),
        dict(species='SO2, v = 0', rest_GHz=267.537, qns='x', Eu_K=200., log10_A_s=-4., S_mu2=10., seen_in_space=True)]))
    r = str(tmp_path/'run')
    main(['data', '--run', r, '--local', str(data), '--target', 'SYNTH-CON'])
    main(['inspect', '--run', r, '--z', str(Z)])
    assert not any('continuum' in fl for c in Run(r).load('inventory')['cubes'] for fl in c['flags'])   # broad lines, no continuum
    main(['detect', '--run', r])
    main(['extract', '--run', r])
    run = Run(r)
    feats = run.load('features')['features']
    good = [f for f in feats if f['significant'] and f['reliable']]
    assert {f['spw'] for f in good} == {'spw25', 'spw27'}
    main(['catalog', '--run', r, '--z', str(Z), '--source', f'file:{catalogue}'])
    main(['identify', '--run', r])
    cands = run.load('candidates')['features']
    top = {cands[f['id']]['candidates'][0]['species'] for f in good}
    assert top == {'HCN, v = 0', 'HCO+, v = 0'}
    hcn = next(f['id'] for f in good if f['spw'] == 'spw27'); hcop = next(f['id'] for f in good if f['spw'] == 'spw25')
    hyps = [dict(name='no vib', shape='gauss2', features=[hcn, hcop],
                 lines=[dict(feature=hcn, species='HCN, v = 0', rest_GHz=HCN), dict(feature=hcop, species='HCO+, v = 0', rest_GHz=HCOP)]),
            dict(name='vib tied', shape='gauss2', features=[hcn, hcop],
                 lines=[dict(feature=hcn, species='HCN, v = 0', rest_GHz=HCN),
                        dict(feature=hcn, species='HCN, v2 = 1', rest_GHz=VIB_E, tie=2, ratio=1),
                        dict(feature=hcop, species='HCN, v2 = 1', rest_GHz=VIB_F), dict(feature=hcop, species='HCO+, v = 0', rest_GHz=HCOP)])]
    hfile = tmp_path/'h.json'; hfile.write_text(json.dumps(hyps))
    main(['blend', '--run', r, '--hypothesis', str(hfile)])
    fit = run.load('blend_vib_tied')
    vib = next(q for q in fit['components'] if abs(q['row']['rest_GHz']-VIB_F) < 1e-4)
    assert vib['verdict'] == 'Needed' and vib['flux'] > 5*vib['flux_error']
    assert run.load('blend_no_vib')['bic'] > fit['bic']
    main(['robust', '--run', r, '--hypotheses', str(hfile), '--n', '4', '--name', 'vib'])
    assert run.load('robust_vib')['win_fraction']['vib tied'] == 1.
    solution = dict(lines=[
        dict(feature=hcn, status='secure', main=dict(species='HCN, v = 0', rest_GHz=HCN, qns='J=3-2'),
             blend=[dict(species='HCN, v2 = 1', rest_GHz=VIB_E, qns='l=1e', status='tentative')], blend_fit='vib tied', evidence='test'),
        dict(feature=hcop, status='secure', main=dict(species='HCO+, v = 0', rest_GHz=HCOP, qns='J=3-2'),
             blend=[dict(species='HCN, v2 = 1', rest_GHz=VIB_F, qns='l=1f')], blend_fit='vib tied',
             alternatives=[dict(solution='no vib', why_not='worse BIC')])])
    sfile = tmp_path/'s.json'; sfile.write_text(json.dumps(solution))
    main(['solve', '--run', r, '--solution', str(sfile)])
    sol = run.load('solution')
    assert not sol['problems']
    assert [l['status'] for e in sol['lines'] for l in e['lines']] == ['secure', 'tentative', 'secure', 'secure']
    csv = (run.path/'products'/'line_solution.csv').read_text().splitlines()
    assert len(csv) == 5 and 'tentative' in csv[2]
    main(['maps', '--run', r])
    maps = run.load('maps')
    assert any(p.endswith('mom0.fits') for p in maps['products']) and any('blendfit' in p for p in maps['products'])
    (run.path/'notes'/'summary.md').write_text('**Synthetic** test run.')
    main(['report', '--run', r, '--title', 'Synthetic CON'])
    html = (run.path/'weblog'/'index.html').read_text()
    assert __import__('re').search(r'<(b|strong)>Synthetic</(b|strong)>', html) and 'HCN, v2 = 1' in html and 'tentative' in html
    for src in __import__('re').findall(r'src="(plots/[^"]+)"', html):
        assert (run.path/'weblog'/src).exists(), src
    main(['status', '--run', r])
    assert 'solution' in capsys.readouterr().out


def test_ads_standard_queries_avoid_unsupported_object_field():
    from conlines.agent.context import standard_queries
    queries = standard_queries('IRAS 17208-0014', ['IRAS F17207-0014'])
    assert all('object:' not in q for q in queries.values())
    assert all('abs:"IRAS 17208-0014"' in q for q in queries.values())


def test_ned_reads_current_column_names(monkeypatch):
    from astropy.table import Table
    from astroquery.ipac.ned import Ned
    from conlines.agent import context
    obj = Table(dict(zip(['Object Name', 'RA', 'Dec', 'Redshift (z)', 'Phys Type'],
                         [['X'], [260.84], [-0.28], [0.04281], ['G']])))
    phot = Table(dict(zip(['Bandpass', 'Frequency (Hz)', 'Flux Density', 'Reference'],
                          [['60 microns (IRAS)'], [5e12], [32.13], ['ref']])))
    monkeypatch.setattr(Ned, 'query_object', lambda name: obj)
    monkeypatch.setattr(Ned, 'get_table', lambda name, table: phot)
    out = context.ned('X')
    assert out['redshift'] == 0.04281 and out['type'] == 'G' and out['dec'] == -0.28
    assert context.iras_fluxes(out['photometry']) == {60: 32.13}
    assert 'photometry_error' not in out


def test_discover_keeps_tunings_that_reuse_spw_numbers(tmp_path):
    from conlines.agent.data import discover
    for mous in ('X2143', 'X2157'):
        for spw in (25, 27):
            (tmp_path/f'member.uid___A001_X2d20_{mous}.IC_860_sci.spw{spw}.cube.I.pbcor.fits').write_bytes(b'')
    cubes, _, _ = discover(tmp_path)
    assert sorted(c['key'] for c in cubes) == ['X2143-spw25', 'X2143-spw27', 'X2157-spw25', 'X2157-spw27']


def test_reliability_untestable_when_control_zones_are_tiny():
    from conlines.reliability import detection_reliability
    props = [dict(sign='Absorption', peak_snr=10., channels=8), dict(sign='Absorption', peak_snr=12., channels=9)]
    out = detection_reliability(props, [0.01, 0.30], 0.2, dict(inside=1000, near=50, far=0))
    assert out[0]['reliability'] is None          # 50-pixel zone cannot calibrate a 1000-pixel region
    out = detection_reliability(props, [0.01, 0.30], 0.2, dict(inside=1000, near=5000, far=0))
    assert out[0]['reliability'] is not None


def test_combined_weblog_merges_high_level_sections(tmp_path):
    from conlines.agent.rundir import Run
    from conlines.agent import report
    runs = []
    for name in ('a', 'b'):
        d = tmp_path/name
        for sub in ('plots', 'products', 'notes'):
            (d/sub).mkdir(parents=True)
        (d/'notes'/'summary.md').write_text(f'Own summary of {name}.')
        (d/'notes'/'detection.md').write_text(f'Detection notes of {name}.')
        runs.append(Run(d))
    nd = tmp_path/'merged'; nd.mkdir()
    (nd/'summary.md').write_text('Merged summary.'); (nd/'context.md').write_text('Merged context.')
    out = report.build_combined([(runs[0], 'Band A', 'ba'), (runs[1], 'Band B', 'bb')], tmp_path/'web', 'Combined', nd)
    page = out.read_text()
    assert page.count('Merged summary.') == 1 and 'Own summary of a.' not in page       # high-level text once, not per band
    assert 'id="summary"' in page and 'id="context"' in page and 'id="provenance"' in page
    assert 'id="ba-detection"' in page and 'id="bb-detection"' in page and 'id="ba-summary"' not in page
    assert 'Detection notes of b.' in page and (tmp_path/'web'/'bb'/'plots').is_dir()
