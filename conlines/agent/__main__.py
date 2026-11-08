"""Command line: python -m conlines.agent <step> --run RUN_FOLDER [options]   (or scripts/agent ...)

Steps, in the usual order:
  data       register cubes: --local DIR, or --alma PROJECT [--target NAME] [--list | --download --dest DIR]
  inspect    describe every cube (coverage, beam, noise, continuum, brightest emission)
  context    SIMBAD + NED + derived properties + ADS literature (--name NAME, --extra-query Q ...)
  ads        free-form ADS search (--query Q), stored with the literature
  detect     ESPADA/SoFiA on every cube
  extract    aperture spectra, significance and reliability of every feature
  zfrom      redshift implied by a feature and a rest frequency (--feature F --rest GHZ)
  catalog    candidate transitions at a redshift (--z Z)
  plan       species screen: which species to test and why (write notes/strategy.md after it)
  identify   candidates per feature with companion-line (LTE) tests
  blend      fit a blend hypothesis (--hypothesis FILE_OR_JSON)
  robust     compare hypotheses under noise realisations (--hypotheses FILE_OR_JSON --n 100)
  solve      record the adopted solution and derive line properties (--solution FILE_OR_JSON)
  maps       moment maps and per-line blend maps of the solution
  report     write the weblog (weblog/index.html)
  status     what has been done in this run
"""
import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _env():
    sofia = ROOT/'.local-tools'/'sofia'
    if not os.environ.get('SOFIA_EXECUTABLE') and sofia.exists():
        os.environ['SOFIA_EXECUTABLE'] = str(sofia)


def main(argv=None):
    _env()
    p = argparse.ArgumentParser(prog='agent', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='step', required=True)
    def cmd(name, **kw):
        s = sub.add_parser(name, **kw); s.add_argument('--run', required=True, help='run folder'); return s
    s = cmd('data'); s.add_argument('--local'); s.add_argument('--alma'); s.add_argument('--target'); s.add_argument('--list', action='store_true')
    s.add_argument('--download', action='store_true'); s.add_argument('--dest'); s.add_argument('--include-pb', action='store_true')
    s = cmd('inspect'); s.add_argument('--z', type=float)
    s = cmd('context'); s.add_argument('--name'); s.add_argument('--z', type=float, help='redshift to use if NED/SIMBAD give none')
    s.add_argument('--no-literature', action='store_true'); s.add_argument('--extra-query', action='append', default=[]); s.add_argument('--rows', type=int, default=15)
    s = cmd('ads'); s.add_argument('--query', required=True); s.add_argument('--rows', type=int, default=20); s.add_argument('--sort', default='citation_count desc')
    s = cmd('detect'); s.add_argument('--widths', type=int, nargs='+', default=[3, 5, 9]); s.add_argument('--clip', type=float, default=5.)
    s.add_argument('--mode', default='Both', choices=['Both', 'Emission', 'Absorption']); s.add_argument('--no-contsub', action='store_true')
    s.add_argument('--noise', default='local', choices=['local', 'spectral']); s.add_argument('--only', nargs='+')
    s = cmd('extract'); s.add_argument('--position', default='auto', choices=['auto', 'line', 'continuum', 'centre', 'radec'])
    s.add_argument('--ra', type=float); s.add_argument('--dec', type=float); s.add_argument('--radius', type=float, help='arcsec (default: one beam)')
    s.add_argument('--threshold', type=float, default=8.); s.add_argument('--region-beams', type=float, default=5.); s.add_argument('--min-reliability', type=float, default=.5)
    s = cmd('zfrom'); s.add_argument('--feature', required=True); s.add_argument('--rest', type=float, required=True)
    s = cmd('catalog'); s.add_argument('--z', type=float, required=True); s.add_argument('--filters', nargs='+', default=['extragalactic', 'ism-hot-core'])
    s.add_argument('--emax', type=float, default=1500.); s.add_argument('--min-log-a', type=float, default=-6.5); s.add_argument('--margin', type=float, default=600.)
    s.add_argument('--source', default='auto', help="auto (Splatalogue, else local CDMS/JPL), splatalogue, local, or file:PATH (JSON list of transitions)")
    s = cmd('plan'); s.add_argument('--expect', action='append', default=[], metavar='SPECIES',
                                    help="expected species, e.g. 'HCN' or 'HCN v2=1' (repeat, or separate with ';'); or @FILE with one per line")
    s.add_argument('--t-ref', type=float, default=150., help='excitation temperature (K) for line strengths in the screen')
    s.add_argument('--pmax', type=float, default=.3, help='chance-coincidence probability above which a species is a line forest')
    s.add_argument('--margin', type=float, default=100.)
    s = cmd('identify'); s.add_argument('--margin', type=float, default=100.); s.add_argument('--include', default='reliable', choices=['reliable', 'significant', 'all'])
    s.add_argument('--top', type=int, default=6)
    s.add_argument('--species', default='plan', choices=['plan', 'all'], help="plan: species from the `plan` screen (default); all: every catalogue species")
    s = cmd('blend'); s.add_argument('--hypothesis', required=True)
    s = cmd('robust'); s.add_argument('--hypotheses', required=True); s.add_argument('--n', type=int, default=100); s.add_argument('--name', default='robustness')
    s = cmd('solve'); s.add_argument('--solution', required=True)
    s = cmd('maps'); s.add_argument('--baseline', default='linear', choices=['linear', 'constant', 'none']); s.add_argument('--fixed-shape', action='store_true')
    s.add_argument('--zoom', type=float, default=4.)
    s = cmd('report'); s.add_argument('--title')
    s = cmd('combine'); s.add_argument('--part', action='append', required=True, metavar='RUN:LABEL[:ID]',
                                       help='a run folder to include, its heading and a short id (repeat)')
    s.add_argument('--title', required=True)
    s.add_argument('--notes', help='folder with merged notes: summary.md, data.md, context.md, provenance.md')
    cmd('status')
    a = p.parse_args(argv)

    from .rundir import Run
    if a.step == 'combine':                 # --run is the output folder here, not an analysis run
        from . import report
        parts = []
        for i, spec in enumerate(a.part):
            bits = spec.split(':')
            parts.append((Run(bits[0]), bits[1] if len(bits) > 1 else bits[0], bits[2] if len(bits) > 2 else f'part{i+1}'))
        print(report.build_combined(parts, a.run, a.title, a.notes))
        return
    run = Run(a.run)

    if a.step == 'data':
        from . import data
        if a.local:
            cubes, continuum, notes = data.discover(a.local)
            if not cubes:
                sys.exit(f'No cubes found in {a.local}.')
            targets = sorted({c['target'] for c in cubes if c['target']})
            run.update(cubes=cubes, continuum=continuum, data_source=dict(local=str(Path(a.local).resolve())), target=a.target or (targets[0] if targets else None))
            run.log('data', f'{len(cubes)} cubes registered', notes=notes)
            print('\n'.join(f"{c['key']}: {c['path']}" + (' (+pb)' if c['pb'] else '') for c in cubes))
            print('\n'.join('! '+n for n in notes))
        elif a.alma:
            obs = data.alma_observations(a.alma, a.target)
            print(f'{len(obs)} observation(s) of {a.alma}' + (f' matching {a.target}' if a.target else '') + ':')
            for o in obs:
                print(f"  {o['target']}  {o['member_ous']}  band {o['band']}  res {o['resolution_arcsec']:.2f}\"  released {o['release'][:10]}")
            products = []
            for o in obs:
                prods = [x for x in data.alma_products(o['member_ous']) if x['role'] in (('cube', 'pb') if a.include_pb else ('cube',))
                         and (not a.target or (x['target'] or '').replace(' ', '').lower() == a.target.replace(' ', '').replace('_', '').lower() or not x['target'])]
                products += prods
                for x in prods:
                    print(f"    {x['role']:<5} {x['spw']}  {x['size']/1e9 if x['size'] else float('nan'):.2f} GB  {x['name']}")
            run.update(alma=dict(project=a.alma, target=a.target, observations=obs, products=products))
            if a.download:
                if not a.dest:
                    sys.exit('--download needs --dest DIR')
                total = sum(x['size'] or 0 for x in products)/1e9
                print(f'Downloading {len(products)} files ({total:.1f} GB) to {a.dest} ...')
                data.alma_download([x['url'] for x in products], a.dest)
                cubes, continuum, notes = data.discover(a.dest)
                run.update(cubes=cubes, continuum=continuum, data_source=dict(alma=a.alma, local=str(Path(a.dest).resolve())), target=a.target)
                print('\n'.join(f"{c['key']}: {c['path']}" for c in cubes))
        else:
            sys.exit('Give --local DIR or --alma PROJECT.')
    elif a.step == 'inspect':
        from . import inspect_step
        if a.z is not None:
            run.update(z=a.z)
        print(inspect_step.summary(inspect_step.run(run, a.z if a.z is not None else run.state.get('z'))))
    elif a.step == 'context':
        from . import context, figs
        inv = run.load('inventory')['cubes']
        name = a.name or run.state.get('target') or next((r['object'] for r in inv if r.get('object')), None)
        out = dict(name=name)
        for label, fn in (('simbad', lambda: context.simbad(name) if name else context.simbad(ra=inv[0]['centre_ra'], dec=inv[0]['centre_dec'])),
                          ('ned', lambda: context.ned(name))):
            try:
                out[label] = fn()
            except Exception as e:
                out[label] = None; out[f'{label}_error'] = str(e)
        z = (out.get('ned') or {}).get('redshift') or (out.get('simbad') or {}).get('redshift')
        out['z_source'] = 'NED' if (out.get('ned') or {}).get('redshift') else 'SIMBAD' if z else None
        if not z and (a.z or run.state.get('z')):
            z = a.z or run.state.get('z'); out['z_source'] = 'given (--z)' if a.z else 'run state'
        out['derived'] = context.derived(z, (out.get('ned') or {}).get('photometry'), max(r['beam_major_arcsec'] for r in inv)) if z else {}
        if (out.get('ned') or {}).get('photometry'):
            figs.sed(run.plot('sed'), out['ned']['photometry'], f"{name}: NED photometry")
        run.save('context', out)
        names = [n for n in [name, (out.get('simbad') or {}).get('main_id')] if n]
        aliases = [i.strip() for i in (out.get('simbad') or {}).get('ids', []) if i.strip().startswith(('IRAS', 'NGC', 'IC ', 'UGC', 'Arp', 'Mrk', 'ESO', 'Zw', 'CGCG'))][:3]
        if not a.no_literature:
            lit = context.literature(list(dict.fromkeys(names+aliases)), a.extra_query, a.rows)
            run.save('literature', lit)
        else:
            lit = run.load('literature', {'queries': {}})
        run.update(target=name)
        run.log('context', f'{name}: z={z}', simbad_error=out.get('simbad_error'), ned_error=out.get('ned_error'))
        print(json.dumps({k: v for k, v in out.items() if k != 'ned'} | {'ned': {k: v for k, v in (out.get('ned') or {}).items() if k != 'photometry'}}, indent=1, default=str)[:3000])
        for label, q in lit['queries'].items():
            print(f"\n== {label} ({len(q.get('papers', []))}){' ERROR '+q['error'] if q.get('error') else ''}")
            for d in q.get('papers', [])[:a.rows]:
                print(f"  {d['year']} {d['first_author']}: {d['title'][:110]} [{d['bibcode']}] {d['arxiv_url'] or ''}")
    elif a.step == 'ads':
        from . import context
        docs = context.ads_search(a.query, a.rows, a.sort)
        lit = run.load('literature', {'queries': {}})
        lit.setdefault('queries', {})[f'agent: {a.query[:60]}'] = dict(query=a.query, papers=docs)
        run.save('literature', lit)
        for d in docs:
            print(f"{d['year']} {d['first_author']}: {d['title']} [{d['bibcode']}] cited {d['citations']} {d['arxiv_url'] or ''}\n    {d['abstract'][:600]}\n")
    elif a.step == 'detect':
        from . import detect_step
        print(detect_step.summary(detect_step.run(run, tuple(a.widths), a.clip, a.mode, not a.no_contsub, a.noise, a.only)))
    elif a.step == 'extract':
        from . import extract_step
        how = 'radec' if a.ra is not None else a.position
        print(extract_step.summary(extract_step.run(run, how, a.ra, a.dec, a.radius, a.threshold, a.region_beams, a.min_reliability)))
    elif a.step == 'zfrom':
        f = next(f for f in run.load('features')['features'] if f['id'] == a.feature)
        z = a.rest/f['freq_GHz']-1
        print(f"{a.feature} at {f['freq_GHz']:.6f} GHz as {a.rest} GHz: z = {z:.6f} (cz = {z*299792.458:.0f} km/s, optical); "
              f"±{f.get('centroid_error_kms') or 0:.0f} km/s centroid error")
    elif a.step == 'catalog':
        from . import catalog_step
        print(catalog_step.summary(catalog_step.run(run, a.z, a.filters, a.margin, a.emax, a.min_log_a, a.source)))
    elif a.step == 'plan':
        from . import plan_step
        expect = []
        for item in a.expect:
            if item.startswith('@'):
                expect += [x.split('#')[0].strip() for x in Path(item[1:]).expanduser().read_text().splitlines()]
            else:
                expect += [x.strip() for x in item.split(';')]
        print(plan_step.summary(plan_step.run(run, [x for x in expect if x], a.t_ref, a.pmax, a.margin)))
        if not (run.path/'notes'/'strategy.md').exists():
            print('\nNext: write notes/strategy.md (the identification strategy for reviewers) before running `identify`.')
    elif a.step == 'identify':
        from . import identify_step
        if a.species == 'plan' and not (run.path/'notes'/'strategy.md').exists():
            print('  ! notes/strategy.md is missing: explain the identification strategy before identifying (see the skill, step 8).')
        print(identify_step.summary(identify_step.run(run, a.margin, a.include, top=3, species=a.species), a.top))
    elif a.step == 'blend':
        from . import blend_step
        for h in blend_step.load_hypotheses(a.hypothesis):
            print(blend_step.summary_fit(blend_step.fit(run, h)))
    elif a.step == 'robust':
        from . import blend_step
        print(blend_step.summary_robust(blend_step.robustness(run, blend_step.load_hypotheses(a.hypotheses), a.n, name=a.name)))
    elif a.step == 'solve':
        from . import solution_step
        print(solution_step.summary(solution_step.run(run, a.solution)))
    elif a.step == 'maps':
        from . import maps_step
        print(maps_step.summary(maps_step.run(run, a.baseline, a.zoom, not a.fixed_shape)))
    elif a.step == 'report':
        from . import report
        print(report.build(run, a.title))
    elif a.step == 'status':
        state = run.state
        print(f"run: {run.path}\ntarget: {state.get('target')}  z: {state.get('z')}  cubes: {len(state.get('cubes', []))}")
        for name in ('inventory', 'context', 'literature', 'components', 'features', 'catalogue', 'candidates', 'solution', 'maps'):
            print(f"  {'✓' if run.has(name) else '·'} {name}")
        print('  blends: '+', '.join(p.stem for p in run.path.glob('blend_*.json')))
        print('  notes: '+', '.join(p.stem for p in (run.path/'notes').glob('*.md')))
        if (run.path/'weblog'/'index.html').exists():
            print(f"  weblog: {run.path/'weblog'/'index.html'}")


if __name__ == '__main__':
    main()
