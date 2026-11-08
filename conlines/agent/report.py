"""Weblog: a static HTML report of a run folder (open weblog/index.html in a browser).

Sections follow the analysis: summary -> data -> source context -> detection -> significance ->
identification strategy -> identification -> blends & robustness -> solution & properties -> maps -> provenance. Each
section shows the step's results and plots, and the agent's reasoning from notes/<section>.md.
"""
from pathlib import Path
from datetime import datetime, timezone
import html
import json
import re
import shutil

SECTIONS = [('summary', 'Summary'), ('data', 'Data'), ('context', 'Source context & literature'), ('detection', 'Line search'),
            ('significance', 'Significance'), ('strategy', 'Identification strategy'), ('identification', 'Identification'),
            ('blends', 'Blends & robustness'), ('solution', 'Line solution & properties'), ('maps', 'Maps'), ('provenance', 'Provenance')]

CSS = """
:root{--ink:#173047;--muted:#5b6b7a;--line:#dbe3ea;--bg:#f7f9fc;--card:#fff;--accent:#138c87;--good:#166534;--warn:#92400e;--bad:#991b1b}
*{box-sizing:border-box}body{margin:0;font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;color:var(--ink);background:var(--bg)}
nav{position:fixed;top:0;left:0;bottom:0;width:230px;background:#0f2a3d;color:#cfe3ea;padding:20px 14px;overflow:auto}
nav b{display:block;color:#fff;font-size:16px;margin-bottom:12px}nav a{display:block;color:#cfe3ea;text-decoration:none;padding:5px 8px;border-radius:6px;font-size:14px}
nav a:hover{background:#1d4258}main{margin-left:230px;padding:28px 40px;max-width:1250px}
section{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 24px;margin:0 0 22px}
h1{margin:0 0 4px;font-size:26px}h2{margin:0 0 12px;font-size:20px;border-bottom:2px solid var(--accent);padding-bottom:6px}h3{font-size:16px;margin:18px 0 8px}
.sub{color:var(--muted);margin-bottom:20px}table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0 14px}
th,td{border-bottom:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top}th{background:#eef3f7;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}img{max-width:100%;border:1px solid var(--line);border-radius:6px;margin:6px 0 12px;background:#fff}
.tag{display:inline-block;padding:1px 8px;border-radius:10px;font-size:12px;font-weight:600}
.secure,.good{background:#dcfce7;color:var(--good)}.tentative,.warn{background:#fef3c7;color:var(--warn)}.unidentified,.bad{background:#fee2e2;color:var(--bad)}
.notes{background:#f3f7fa;border-left:4px solid var(--accent);padding:8px 16px;border-radius:6px;margin:8px 0 14px}
.notes p:first-child{margin-top:4px}details{margin:8px 0}summary{cursor:pointer;font-weight:600}code{background:#eef3f7;padding:1px 4px;border-radius:4px}
.missing{color:var(--bad)}.small{font-size:12px;color:var(--muted)}
@media (max-width:900px){nav{position:static;width:auto}main{margin-left:0;padding:16px}}
"""


# ---------------------------------------------------------------- minimal Markdown
def _indent_lists(text):
    """Python-Markdown nests lists only at 4-space steps; notes are usually written with 2."""
    out, prev_blank, in_code = [], True, False
    for line in text.splitlines():
        if line.lstrip().startswith('```'):
            in_code = not in_code
        n = len(line)-len(line.lstrip(' '))
        if not in_code and n and not prev_blank:
            line = ' '*(2*n)+line.lstrip(' ')
        out.append(line); prev_blank = not line.strip()
    return '\n'.join(out)


def markdown(text):
    try:
        import markdown as md
        return md.markdown(_indent_lists(text), extensions=['tables', 'fenced_code', 'sane_lists'])
    except ImportError:
        pass
    out, para, in_list, in_code = [], [], None, False
    def inline(s):
        s = html.escape(s)
        s = re.sub(r'`([^`]+)`', r'<code>\1</code>', s)
        s = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', s)
        s = re.sub(r'\*([^*]+)\*', r'<i>\1</i>', s)
        s = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', s)
        return s
    def flush():
        nonlocal para
        if para:
            out.append('<p>'+inline(' '.join(para))+'</p>'); para = []
    table = []
    def flush_table():
        nonlocal table
        if table:
            rows = [[c.strip() for c in r.strip().strip('|').split('|')] for r in table if not re.match(r'^\s*\|?[\s:|-]+\|?\s*$', r)]
            head, body = rows[0], rows[1:]
            out.append('<table><thead><tr>'+''.join(f'<th>{inline(c)}</th>' for c in head)+'</tr></thead><tbody>'
                       + ''.join('<tr>'+''.join(f'<td>{inline(c)}</td>' for c in r)+'</tr>' for r in body)+'</tbody></table>')
            table = []
    for line in text.splitlines():
        if line.lstrip().startswith('|') and not in_code:
            flush(); table.append(line); continue
        flush_table()
        if in_list and line.startswith('  ') and line.strip() and not re.match(r'^\s*(?:[-*]|\d+\.)\s+', line):
            out[-1] = out[-1][:-5]+' '+inline(line.strip())+'</li>'; continue
        if line.startswith('```'):
            flush(); out.append('</pre>' if in_code else '<pre>'); in_code = not in_code; continue
        if in_code:
            out.append(html.escape(line)); continue
        m = re.match(r'^(#{1,4})\s+(.*)', line)
        item = re.match(r'^\s*(?:[-*]|\d+\.)\s+(.*)', line)
        if m:
            flush(); out.append(f'<h{len(m.group(1))+2}>{inline(m.group(2))}</h{len(m.group(1))+2}>')
        elif item:
            flush()
            if not in_list:
                out.append('<ul>'); in_list = True
            out.append('<li>'+inline(item.group(1))+'</li>')
        elif not line.strip():
            flush()
            if in_list:
                out.append('</ul>'); in_list = False
        else:
            if in_list:
                out.append('</ul>'); in_list = False
            para.append(line.strip())
    flush(); flush_table()
    if in_list:
        out.append('</ul>')
    return '\n'.join(out)


# ---------------------------------------------------------------- helpers
def e(v):
    return html.escape('' if v is None else str(v))


def num(v, fmt='{:.3g}'):
    if v is None:
        return ''
    try:
        return fmt.format(v)
    except (ValueError, TypeError):
        return e(v)


def table(rows, columns):
    """columns: list of (header, function(row) -> html, numeric?)"""
    head = ''.join(f'<th>{e(h)}</th>' for h, _, _ in columns)
    body = ''.join('<tr>'+''.join(f"<td{' class=num' if n else ''}>{f(r)}</td>" for _, f, n in columns)+'</tr>' for r in rows)
    return f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def img(plots, name, alt=''):
    p = plots/name
    return f'<img src="plots/{e(name)}" alt="{e(alt or name)}" loading="lazy">' if p.exists() else ''


_SKIP_NOTES = set()          # note keys left out (the combined weblog shows merged notes instead)


def notes(run, key):
    if key in _SKIP_NOTES:
        return ''
    p = run.path/'notes'/f'{key}.md'
    if not p.exists():
        return ''
    return f'<div class="notes">{markdown(p.read_text())}</div>'


def tag(status):
    return f'<span class="tag {e(status)}">{e(status)}</span>'


# ---------------------------------------------------------------- sections
def section_summary(run, out):
    state = run.state
    sol = run.load('solution', {'lines': []}) if run.has('solution') else {'lines': []}
    rows = []
    for ent in sol['lines']:
        for r in ent['lines'] or [dict(role='', species='unidentified', qns='', rest_GHz=None, flux_Jybkms=None, flux_error_Jybkms=None)]:
            rows.append(dict(ent, **{'line': r, 'status': r.get('status', ent['status'])}))
    out.append(notes(run, 'summary'))
    if rows:
        out.append('<h3>Adopted line solution</h3>'+table(rows, [
            ('Feature', lambda r: e(r['feature']), False), ('Status', lambda r: tag(r['status']), False),
            ('Line', lambda r: e(r['line']['species'])+(' <span class=small>(blended)</span>' if r['line'].get('role') == 'blended' else ''), False),
            ('Transition', lambda r: e(r['line'].get('qns')), False), ('Rest (GHz)', lambda r: num(r['line'].get('rest_GHz'), '{:.5f}'), True),
            ('Observed (GHz)', lambda r: num(r['observed_GHz'], '{:.5f}'), True),
            ('Flux (Jy/beam km/s)', lambda r: num(r['line'].get('flux_Jybkms'))+(' ± '+num(r['line'].get('flux_error_Jybkms'), '{:.2g}') if r['line'].get('flux_error_Jybkms') is not None else ''), True),
            ('S/N', lambda r: num(r['snr'], '{:.1f}'), True), ('Reliability', lambda r: num(r['reliability'], '{:.2f}'), True)]))
    target = state.get('target') or ''
    out.append(f"<p class=small>Target: {e(target)} · z = {num(state.get('z'), '{:.6f}')} · run folder {e(run.path)}</p>")


def section_data(run, out, plots):
    inv = run.load('inventory', {'cubes': []})
    out.append(notes(run, 'data'))
    out.append(table(inv['cubes'], [
        ('Window', lambda r: e(r['key']), False), ('Observed GHz', lambda r: f"{r['freq_min_GHz']:.3f}–{r['freq_max_GHz']:.3f}", False),
        ('Channels', lambda r: f"{r['nchan']} × {r['channel_kms']:.1f} km/s", False),
        ('Beam (″)', lambda r: f"{r['beam_major_arcsec']:.2f} × {r['beam_minor_arcsec']:.2f}", False),
        ('Noise / channel', lambda r: f"{r['noise_median']:.2e} {e(r['unit'])}", True),
        ('Brightest line S/N', lambda r: num(r['line_peak_snr'], '{:.0f}'), True), ('Continuum S/N', lambda r: num(r['continuum_snr'], '{:.0f}'), True),
        ('Flags', lambda r: '<br>'.join(e(f) for f in r['flags']), False)]))
    out.append(img(plots, 'coverage.png'))
    for r in inv['cubes']:
        out.append(f"<h3>{e(r['key'])}</h3>"+img(plots, f"inspect_{r['key']}.png"))
    # optional continuum-subtraction diagnostics made by the agent
    extra = sorted(p.name for p in plots.glob('contsub_*.png'))+['linefree_residual_maps.png', 'b7_binned_overview.png']
    shown = [img(plots, n) for n in extra if (plots/n).exists()]
    if shown:
        out.append('<h3>Continuum-subtraction diagnostics</h3>'+''.join(shown))


def section_context(run, out, plots):
    ctx = run.load('context', {})
    out.append(notes(run, 'context'))
    if ctx:
        s, n, d = ctx.get('simbad') or {}, ctx.get('ned') or {}, ctx.get('derived') or {}
        facts = [('SIMBAD name', s.get('main_id')), ('Object type', s.get('otype')), ('Morphology', s.get('morph')),
                 ('NED type', n.get('type')), ('Redshift (SIMBAD / NED)', f"{num(s.get('redshift'), '{:.5f}')} / {num(n.get('redshift'), '{:.5f}')}"),
                 ('Luminosity distance', num(d.get('luminosity_distance_Mpc'), '{:.1f} Mpc')), ('Scale', num(d.get('pc_per_arcsec'), '{:.0f} pc/″')),
                 ('Beam', num(d.get('beam_pc'), '{:.0f} pc')), ('L_IR (IRAS)', num(d.get('L_IR_Lsun'), '{:.2e} L☉')+(f" ({d.get('class_')})" if d.get('class_') else '')),
                 ('IRAS f25/f60', num(d.get('f25_f60'), '{:.2f}'))]
        out.append(table([dict(k=k, v=v) for k, v in facts if v not in (None, '', ' / ')], [('Property', lambda r: e(r['k']), False), ('Value', lambda r: e(r['v']), False)]))
        out.append(img(plots, 'sed.png'))
    lit = run.load('literature', {})
    for label, q in (lit.get('queries') or {}).items():
        papers = q.get('papers', [])
        out.append(f"<details><summary>{e(label)} ({len(papers)})</summary><p class=small>ADS query: <code>{e(q.get('query'))}</code>{' — '+e(q['error']) if q.get('error') else ''}</p>")
        out.append(table(papers, [('Year', lambda r: e(r['year']), False), ('Paper', lambda r: f"<a href=\"{e(r['url'])}\">{e(r['title'])}</a><br><span class=small>{e(r['first_author'])}{' et al.' if r['n_authors'] > 1 else ''}, {e(r['journal'])}</span>", False),
                                  ('Cited', lambda r: e(r['citations']), True)]))
        out.append('</details>')


def section_detection(run, out, plots):
    comps = run.load('components', {})
    out.append(notes(run, 'detection'))
    rows = [dict(key=k, **v) for k, v in comps.items()]
    out.append(table(rows, [('Window', lambda r: e(r['key']), False), ('Components', lambda r: str(len(r['components'])), True),
                            ('ESPADA status', lambda r: e(r['status']), False), ('Settings', lambda r: e(r['settings']), False)]))
    for k in comps:
        out.append(f'<h3>{e(k)}</h3>'+img(plots, f'detect_{k}.png'))


def section_significance(run, out, plots):
    feats = run.load('features', {})
    out.append(notes(run, 'significance'))
    if not feats:
        return
    p = feats['position']
    out.append(f"<p>Aperture radius {feats['radius_arcsec']:.3f}″ at {e(p['reason'])} (RA {p['ra']:.6f}, Dec {p['dec']:.6f}); "
               f"integrated S/N threshold {feats['threshold']:g}; reliability from the artefact rate outside {feats['region_radius_arcsec']:.2f}″; "
               f"minimum reliability {feats['min_reliability']:g}.</p>")
    def verdict(r):
        cls = 'good' if r['significant'] and r['reliable'] else 'warn' if r['significant'] or r['reliable'] else 'bad'
        return f"<span class='tag {cls}'>{e(r['verdict'])}</span>"
    out.append(table(sorted(feats['features'], key=lambda f: f['freq_GHz']), [
        ('Feature', lambda r: e(r['id']), False), ('Type', lambda r: e(r['sign']), False), ('GHz', lambda r: num(r['freq_GHz'], '{:.5f}'), True),
        ('Width (km/s)', lambda r: num(r['width_kms'], '{:.0f}'), True), ('S/N', lambda r: num(r['snr'], '{:.1f}'), True),
        ('Reliability', lambda r: num(r['reliability'], '{:.2f}'), True), ('Flux (Jy/b km/s)', lambda r: num(r['flux_Jykms'])+' ± '+num(r['flux_error_Jykms'], '{:.2g}'), True),
        ('Verdict', verdict, False), ('Flags', lambda r: e(', '.join(r['flags'])), False)]))
    t = feats.get('trials', {})
    if t:
        out.append('<p class=small>Independent elements searched (look-elsewhere): '+'; '.join(f"{e(k)}: {v['independent_elements']:.2e}" for k, v in t.items())+'</p>')
    for k in feats['spectra']:
        out.append(f'<h3>{e(k)}</h3>'+img(plots, f'spectrum_{k}.png'))
    out.append(img(plots, 'binned_spectra.png'))      # optional diagnostics, shown when the agent made them
    out.append(img(plots, 'bump_maps.png'))


def section_strategy(run, out, plots):
    plan = run.load('plan') if run.has('plan') else None
    out.append(notes(run, 'strategy'))
    if not plan:
        if not out[-1]:
            out.append('<p class=small>No species screen in this run (the `plan` step was not run): every catalogue species was tested.</p>')
        return
    c = plan['counts']
    out.append(f"<p>Species screen of the catalogue for the observed windows ({plan['bandwidth_GHz']:.2f} GHz): "
               f"<b>{c['test']}</b> species tested for every feature, <b>{c['reserve']}</b> in reserve (tried only for features that the tested "
               f"species cannot explain), <b>{c['excluded']}</b> excluded. P<sub>c</sub> is the chance that a feature at a random frequency "
               f"matches the species: the fraction of the band within ±{plan['half_width_kms']:.0f} km/s of one of its strong lines (within 10× of "
               f"its strongest, LTE {plan['T_ref']:g} K). Above {plan['pmax']:.2f} a single match says little: a <i>line forest</i> (many strong "
               f"lines, e.g. large organics) can only be claimed through a full multi-line model, and a <i>low-specificity</i> species (few lines, "
               f"but broad features) only through consistent companions or a strong prior.</p>")
    if plan.get('unmatched_expect'):
        out.append(f"<p class=missing>Expected species with no transition in the observed windows: {e(', '.join(plan['unmatched_expect']))}</p>")
    def label(r):
        t = ' <span class="tag warn">line forest</span>' if r['line_forest'] else ' <span class="tag warn">low specificity</span>' if r.get('low_specificity') else ''
        return e(r['species'])+(f"<br><span class=small>= {e(r['aliases'][0])}</span>" if r.get('aliases') else '')+t
    cols = [('Species', label, False),
            ('Why', lambda r: e(r['reason']), False), ('Lines in band', lambda r: str(r['n_lines']), True),
            ('P_c (chance match)', lambda r: num(r['p_coincidence'], '{:.2f}'), True),
            ('Strongest line', lambda r: f"{e(r['strongest']['qns'])} {r['strongest']['rest_GHz']:.4f} GHz (E<sub>u</sub> {num(r['strongest']['Eu_K'], '{:.0f}')} K)", False)]
    groups = {d: [r for r in plan['species'] if r['decision'] == d] for d in ('test', 'reserve', 'excluded')}
    out.append('<h3>Species tested</h3>'+table(groups['test'], cols))
    out.append(f"<details><summary>Reserve ({len(groups['reserve'])})</summary>{table(groups['reserve'], cols)}</details>")
    out.append(f"<details><summary>Excluded ({len(groups['excluded'])})</summary>{table(groups['excluded'], cols)}</details>")
    out.append('<p class=small>Full screen: <a href="products/species_screen.csv">species_screen.csv</a>.</p>')
    for p in sorted(plots.glob('strategy_*.png')):
        out.append(img(plots, p.name))


def section_identification(run, out, plots):
    cands = run.load('candidates', {})
    out.append(notes(run, 'identification'))
    for fid, d in (cands.get('features') or {}).items():
        f = d['feature']
        out.append(f"<h3>{e(fid)} · {e(f['sign'])} · {f['freq_GHz']:.5f} GHz · {f['width_kms']:.0f} km/s</h3>")
        out.append(notes(run, f'identification_{fid}'))
        if d.get('stage'):
            ex = d.get('excluded_coincident') or []
            out.append(f"<p class=small>Candidates from: {e(d['stage'])}."
                       + (f" Excluded by the species screen but with a line here: {len(ex)} ({e(', '.join(ex[:8]))}{', …' if len(ex) > 8 else ''})." if ex else '') + '</p>')
        out.append(img(plots, f'feature_{fid}.png'))
        out.append(table(d['candidates'][:10], [
            ('Score', lambda c: num(c['score'], '{:.2f}'), True), ('Species', lambda c: e(c['species'])+'<br><span class=small>'+e(', '.join(c['flags']+(['reserve'] if c.get('stage') == 'reserve' else [])+(['line forest'] if c.get('line_forest') else [])+(['low specificity'] if c.get('low_specificity') else [])))+'</span>', False),
            ('Transitions in feature', lambda c: '<br>'.join(f"{e(t['qns'])} {t['rest_GHz']:.5f} (E<sub>u</sub> {num(t['Eu_K'], '{:.0f}')} K)" for t in c['transitions'][:4]), False),
            ('Δv (km/s)', lambda c: num(c['offset_kms'], '{:+.0f}'), True), ('f_max', lambda c: num(c['fmax'], '{:.2f}'), True),
            ('T (K)', lambda c: num(c['best_T'], '{:.0f}'), True),
            ('Companions', lambda c: f"{c.get('companions_supported', 0)} seen / {c.get('companions_testable', 0)} testable" + (f", <span class=missing>{c['companions_missing']} missing</span>" if c['companions_missing'] else ''), False),
            ('Prior', lambda c: num(c['prior'], '{:.2f}'), True)]))
        for i in range(1, 4):
            out.append(img(plots, f'companions_{fid}_{i}.png'))


def section_blends(run, out, plots):
    out.append(notes(run, 'blends'))
    out.append(img(plots, 'hco_hcn_template.png'))    # optional diagnostic
    out.append(img(plots, 'f041_shape_tests.png'))
    for p in sorted(run.path.glob('blend_*.json')):
        b = json.loads(p.read_text())
        out.append(f"<h3>{e(b['name'])}</h3><p>Shared shape: {e(b['shape_text'])} · reduced χ² {b['reduced_chi2']:.2f} · BIC {b['bic']:.1f}"
                   + (f" · <span class=missing>shape at the edge of its range: {e(b['at_bounds'])}</span>" if b['at_bounds'] else '')+'</p>')
        out.append(img(plots, p.stem+'.png'))
        out.append(table(b['components'], [
            ('Line', lambda q: e(b['names'].get(q['key'], q['key'])), False), ('Feature', lambda q: e(b['regions'][q['region']].get('feature', q['region'])), False),
            ('Flux (Jy/b km/s)', lambda q: num(q['flux'])+' ± '+num(q['flux_error'], '{:.2g}'), True), ('Share', lambda q: num(q.get('fraction'), '{:.2f}'), True),
            ('ΔBIC', lambda q: num(q['delta_bic'], '{:.1f}'), True), ('Verdict', lambda q: e(q['verdict']), False),
            ('Trades flux with', lambda q: e(', '.join(b['names'].get(d, d) for d in q['degenerate_with'])), False),
            ('Tied', lambda q: e(f"{b['names'].get(q['tie'][0], q['tie'][0])} × {q['tie'][1]:g}") if q['tie'] else '', False)]))
    for p in sorted(run.path.glob('robust_*.json')):
        r = json.loads(p.read_text())
        out.append(f"<h3>Robustness: {e(r['name'])}</h3><p class=small>{e(r['note'])}; {r['n']} realisations.</p>")
        out.append(img(plots, p.stem+'.png'))
        out.append(table(r['lines'], [('Hypothesis', lambda q: e(q['hypothesis']), False), ('Line', lambda q: e(q['line']), False),
                                      ('Flux', lambda q: num(q['flux']), True), ('Median (noise)', lambda q: num(q['median']), True),
                                      ('16–84%', lambda q: num(q['p16'])+' – '+num(q['p84']), True), ('Same sign', lambda q: num(q['fraction_positive'], '{:.0%}'), True)]))


def section_solution(run, out, plots):
    sol = run.load('solution', {}) if run.has('solution') else {}
    out.append(notes(run, 'solution'))
    for ent in sol.get('lines', []):
        out.append(f"<h3>{e(ent['feature'])} {tag(ent['status'])}</h3>")
        if ent['lines']:
            out.append(table(ent['lines'], [
                ('Role', lambda r: e(r['role']), False), ('Status', lambda r: tag(r.get('status', ent['status'])), False), ('Line', lambda r: e(r['species']), False), ('Transition', lambda r: e(r['qns']), False),
                ('Rest (GHz)', lambda r: num(r['rest_GHz'], '{:.5f}'), True), ('Velocity (km/s)', lambda r: num(r['velocity_kms'], '{:+.0f}'), True),
                ('Flux (Jy/b km/s)', lambda r: num(r['flux_Jybkms'])+(' ± '+num(r['flux_error_Jybkms'], '{:.2g}') if r['flux_error_Jybkms'] is not None else ''), True),
                ('Share', lambda r: num(r['share'], '{:.2f}'), True), ('Aperture flux (Jy km/s)', lambda r: num(r['aperture_flux_Jykms']), True),
                ("L′ (K km/s pc²)", lambda r: num(r['Lprime_K_kms_pc2'], '{:.2e}'), True), ('L (L☉)', lambda r: num(r['L_Lsun'], '{:.2e}'), True)]))
        if ent.get('evidence'):
            out.append(f"<div class=notes>{markdown(ent['evidence'])}</div>")
        if ent.get('alternatives'):
            out.append(table(ent['alternatives'], [('Alternative', lambda a: e(a.get('solution')), False), ('Why not adopted', lambda a: e(a.get('why_not')), False)]))
    if sol.get('problems'):
        out.append('<p class=missing>'+'<br>'.join(e(p) for p in sol['problems'])+'</p>')
    if (run.path/'products'/'line_solution.csv').exists():
        out.append('<p><a href="products/line_solution.csv">line_solution.csv</a></p>')


def section_maps(run, out, plots):
    m = run.load('maps', {}) if run.has('maps') else {}
    out.append(notes(run, 'maps'))
    for p in m.get('plots', []):
        out.append(img(plots, p['plot']))
        for w in p.get('warnings', []):
            out.append(f'<p class=missing>{e(w)}</p>')
    if m.get('products'):
        out.append('<details><summary>FITS products</summary><ul>'+''.join(f'<li><a href="products/{e(p)}">{e(p)}</a></li>' for p in m['products'])+'</ul></details>')


def section_provenance(run, out, plots):
    out.append(notes(run, 'provenance'))
    log = run.path/'log.jsonl'
    if log.exists():
        rows = [json.loads(l) for l in log.read_text().splitlines() if l.strip()]
        out.append(table(rows, [('UTC', lambda r: e(r['utc']), False), ('Step', lambda r: e(r['step']), False), ('Result', lambda r: e(r['message']), False),
                                ('Settings', lambda r: '<span class=small>'+e(json.dumps({k: v for k, v in r.items() if k not in ('utc', 'step', 'message')}))[:400]+'</span>', False)]))
    out.append('<p class=small>Generated by conlines.agent (con-line-search). Detection: ESPADA/SoFiA 2; catalogues: Splatalogue (CDMS/JPL) or local CDMS/JPL; '
               'identification criteria after Snyder et al. (2005), Belloche et al. (2013) and ALCHEMI (Martín et al. 2021); fidelity in the spirit of '
               'ASPECS (González-López et al. 2019).</p>')


def _sections(run, plots, prefix=''):
    """HTML <section> blocks of one run; ids are prefixed (e.g. 'b5-') when several runs share a page."""
    funcs = dict(summary=lambda o: section_summary(run, o), data=lambda o: section_data(run, o, plots), context=lambda o: section_context(run, o, plots),
                 detection=lambda o: section_detection(run, o, plots), significance=lambda o: section_significance(run, o, plots),
                 strategy=lambda o: section_strategy(run, o, plots), identification=lambda o: section_identification(run, o, plots), blends=lambda o: section_blends(run, o, plots),
                 solution=lambda o: section_solution(run, o, plots), maps=lambda o: section_maps(run, o, plots),
                 provenance=lambda o: section_provenance(run, o, plots))
    body = []
    for key, label in SECTIONS:
        part = []
        try:
            funcs[key](part)
        except Exception as ex:      # one broken section must not lose the report
            part.append(f'<p class=missing>Could not build this section: {e(ex)}</p>')
        body.append(f'<section id="{prefix}{key}"><h2>{e(label)}</h2>{"".join(part)}</section>')
    extra = sorted(p for p in (run.path/'notes').glob('*.md') if p.stem not in {k for k, _ in SECTIONS} and not p.stem.startswith('identification_'))
    if extra:
        body.append(f'<section id="{prefix}other"><h2>Other notes</h2>'+''.join(f'<h3>{e(p.stem)}</h3><div class=notes>{markdown(p.read_text())}</div>' for p in extra)+'</section>')
    return body


def _copy_outputs(run, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for sub in ('plots', 'products'):
        target = dest/sub
        if target.is_symlink() or target.is_file():
            target.unlink()
        elif target.exists():
            shutil.rmtree(target)
        shutil.copytree(run.path/sub, target)


def build(run, title=None):
    out_dir = run.path/'weblog'
    _copy_outputs(run, out_dir)
    plots = out_dir/'plots'
    state = run.state
    title = title or f"Line identification: {state.get('target') or run.path.name}"
    body = _sections(run, plots)
    nav = ''.join(f'<a href="#{k}">{e(l)}</a>' for k, l in SECTIONS)
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{e(title)}</title><style>{CSS}</style></head><body><nav><b>{e(title)}</b>{nav}</nav><main>'
            f'<h1>{e(title)}</h1><div class=sub>Weblog generated {stamp}</div>{"".join(body)}</main></body></html>')
    (out_dir/'index.html').write_text(page)
    run.log('report', 'weblog written', path=str(out_dir/'index.html'))
    return out_dir/'index.html'


TOP = [('summary', 'Summary'), ('data', 'Data'), ('context', 'Source context & literature')]
BAND = [('detection', 'Line search'), ('significance', 'Significance'), ('strategy', 'Identification strategy'), ('identification', 'Identification'),
        ('blends', 'Blends & robustness'), ('solution', 'Line solution & properties'), ('maps', 'Maps')]


def build_combined(parts, out_dir, title, notes_dir=None):
    """One weblog for several runs of one source (e.g. one run per ALMA band).

    Structure: merged high-level sections first (Summary, Data, Source context & literature), then one
    sub-index per run with only its line-identification sections, then merged Provenance.
    parts: list of (run, label, short_id). notes_dir: folder with the merged notes summary.md, data.md,
    context.md and provenance.md; a merged note replaces the runs' own notes of that section, so the
    high-level text is not repeated per run. Plots and products of each run go to out_dir/<short_id>/.
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    notes_dir = Path(notes_dir) if notes_dir else None
    merged = lambda key: (notes_dir/f'{key}.md') if notes_dir and (notes_dir/f'{key}.md').exists() else None
    for run, _, sid in parts:
        _copy_outputs(run, out_dir/sid)
    def rel(html_, sid):
        return re.sub(r'(src|href)="(plots|products)/', lambda m: f'{m.group(1)}="{sid}/{m.group(2)}/', html_)
    def one(fn, run, sid, skip=()):
        global _SKIP_NOTES
        part, old = [], _SKIP_NOTES
        _SKIP_NOTES = set(skip)
        try:
            fn(run, part, out_dir/sid/'plots')
        except Exception as ex:
            part.append(f'<p class=missing>Could not build this section: {e(ex)}</p>')
        finally:
            _SKIP_NOTES = old
        return rel(''.join(part), sid)
    fns = dict(summary=lambda r, o, p: section_summary(r, o), data=section_data, context=section_context, detection=section_detection,
               significance=section_significance, strategy=section_strategy, identification=section_identification, blends=section_blends,
               solution=section_solution, maps=section_maps, provenance=section_provenance)
    nav, body = [], []
    def note_html(key):
        m = merged(key)
        return f'<div class=notes>{markdown(m.read_text())}</div>' if m else ''
    # 1. merged high-level sections
    for key, label in TOP:
        html_ = note_html(key)
        if key == 'context':                       # same source: show the first run's facts and literature once
            run, _, sid = parts[0]
            html_ += one(fns[key], run, sid, skip={'context'} if merged('context') else ())
        else:
            for run, lab, sid in parts:            # per-run results under a heading; the run's own note only if no merged one
                skip = {key} if (merged(key) and key == 'summary') else set()
                inner = one(fns[key], run, sid, skip)
                html_ += f'<h3>{e(lab)}</h3>{inner}'
        body.append(f'<section id="{key}"><h2>{e(label)}</h2>{html_}</section>')
        nav.append(f'<a href="#{key}">{e(label)}</a>')
    # 2. one sub-index per run: line identification only
    for run, lab, sid in parts:
        secs, links = [], []
        for key, label in BAND:
            inner = one(fns[key], run, sid)
            if not re.sub(r'<[^>]+>', '', inner).strip():
                continue                           # nothing to show (e.g. no blend fits in this run)
            secs.append(f'<section id="{sid}-{key}"><h2>{e(label)}</h2>{inner}</section>')
            links.append((key, label))
        contents = ' · '.join(f'<a href="#{sid}-{k}">{e(l)}</a>' for k, l in links)
        body.append(f'<div class=band id="{sid}"><h1 class=bandtitle>{e(lab)}</h1><p class=small>Contents: {contents}</p>{"".join(secs)}</div>')
        nav.append(f'<a class=bandlink href="#{sid}">{e(lab)}</a>'+''.join(f'<a class=subl href="#{sid}-{k}">{e(l)}</a>' for k, l in links))
    # 3. merged provenance, with each run's step log
    html_ = note_html('provenance')
    for run, lab, sid in parts:
        inner = one(fns['provenance'], run, sid, skip={'provenance'} if merged('provenance') else ())
        html_ += f'<details><summary>{e(lab)}: run notes and step log</summary>{inner}</details>'
    body.append(f'<section id="provenance"><h2>Provenance</h2>{html_}</section>')
    nav.append('<a class=bandlink href="#provenance">Provenance</a>')
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    css = CSS+'.bandlink{margin-top:14px;font-weight:700;color:#fff!important}.subl{padding-left:20px!important;font-size:13px!important}.bandtitle{margin:30px 0 12px;font-size:24px;border-bottom:3px solid var(--accent);padding-bottom:6px}'
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{e(title)}</title><style>{css}</style></head><body><nav><b>{e(title)}</b>{"".join(nav)}</nav><main>'
            f'<h1>{e(title)}</h1><div class=sub>Weblog generated {stamp}</div>{"".join(body)}</main></body></html>')
    (out_dir/'index.html').write_text(page)
    return out_dir/'index.html'
