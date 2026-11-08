"""Step: identification strategy, i.e. which species are worth testing, and why.

Before any feature is matched to the catalogue, every species with transitions in the observed
windows is screened, and the agent explains its choice in notes/strategy.md. Matching every
feature against every catalogue species gives dozens of "candidates" per feature, most of them
meaningless. The screen keeps the candidate lists short and the reasons explicit.

The screen uses three properties of each species:

* **Expectation.** Is it on the agent's expectation list for this source (from the literature and
  the source type: `--expect`)? Expected species are always tested.
* **Abundance tier.** Species that make bright lines in almost every galaxy (BRIGHT) are tested.
  Species that are widely detected but weaker (COMMON, including the common isotopologues), and
  species seen in other galaxies (Splatalogue 'extragalactic' flag), are kept in reserve: the
  agent tests them for every feature only by putting them on the expectation list, with a reason. Species only seen in space, or never, are not tested
  unless expected.
* **Chance coincidence.** P_c is the fraction of the observed bandwidth that lies within ±w of one
  of the species' strong transitions, i.e. lines predicted (LTE, T_ref) within 10x of its
  strongest line in the band. Here w = max(match margin, half the median feature width). A
  feature at a random frequency "matches" the species with probability P_c, so for P_c ≳ 0.3 a
  frequency match says almost nothing. With at least FOREST_LINES strong lines the species is a
  *line forest*; with fewer, P_c is high only because the features are broad (*low specificity*). This applies to large organics such as (CH3)2CO, whose
  hundreds of lines cover the band like a forest. Such species can only be claimed through a full
  multi-line (LTE) model that reproduces the whole spectrum, never from one feature. They are
  excluded unless expected, and even then the candidate is marked as a line forest.

Vibrationally excited states of BRIGHT/COMMON species, and their rarer isotopologues not on those
lists, go to reserve; the agent should test them once the parent is secure, or put them on the
expectation list when the source calls for them (e.g. HCN v2=1 and HC3N v7=1 in a CON).

Decisions:
* test: tried for every feature;
* reserve: tried only for features that no tested species explains;
* excluded: never tried automatically (the weblog lists them with the reason).
"""
import csv
import re
import numpy as np
from ..core import C
from .identify_step import BRIGHT, COMMON, base_species, kind_flags, unit_intensity, vib_aliases
from . import figs


FOREST_LINES = 10      # a 'line forest' has at least this many strong lines in the observed windows


def norm(name):
    return re.sub(r'[\s,_]+', '', str(name)).lower()


def matches(item, species):
    """Does an expectation-list item ('HCN', 'HCN v2=1', 'HC3N, v7=1') name this catalogue species?

    A bare molecule name means its ground state and hyperfine variants; excited states must be
    named explicitly (as in 'HCN v2=1')."""
    n, s = norm(item), norm(species)
    if n == s or n+'[hf]' == s or n+'v=0' == s or n+'v=0[hf]' == s or any(n == norm(a) for a in vib_aliases(species)):
        return True
    base = norm(base_species(species))
    return n == base and 'vibrationally excited' not in kind_flags(species)


ISOTOPES = [('13C', 'C'), ('15N', 'N'), ('18O', 'O'), ('17O', 'O'), ('34S', 'S'), ('33S', 'S'), ('36S', 'S'), ('29Si', 'Si'), ('30Si', 'Si')]


def parent_of(name):
    """Main isotopologue of a molecule name: 'H13CN' -> 'HCN', 'HC18O+' -> 'HCO+', 'O13CS' -> 'OCS'."""
    for iso, main in ISOTOPES:
        name = name.replace(iso, main)
    name = name.replace('HDO', 'H2O').replace('NH2D', 'NH3')
    return re.sub(r'D(?![a-z])', 'H', name)


def tier_of(species, row):
    b = base_species(species)
    if b in BRIGHT:
        return 'bright'
    if b in COMMON:
        return 'common'
    if row.get('seen_extragalactic'):
        return 'extragalactic'
    if row.get('seen_in_space'):
        return 'seen in space'
    return 'not known in space'


def _coverage(freqs, half_kms, windows):
    """Fraction of the summed window bandwidth within ±half_kms of any of the frequencies."""
    total = covered = 0.
    for lo, hi in windows:
        total += hi-lo
        iv = sorted((max(lo, f*(1-half_kms/C)), min(hi, f*(1+half_kms/C))) for f in freqs if lo-1 < f < hi+1)
        end = lo
        for a, b in iv:
            if b <= a or b <= end:        # outside this window, or already covered
                continue
            covered += b-max(a, end); end = b
    return covered/total if total else 0.


def screen_species(catalogue, windows_obs, z, half_kms, expect=(), T_ref=150., pmax=.3, strong_fraction=.1):
    """One entry per catalogue species with at least one transition inside the observed windows."""
    by_species = {}
    for key, w in catalogue['windows'].items():
        lo, hi = windows_obs[key]
        for r in w['rows']:
            f = r['rest_GHz']/(1+z)
            if lo <= f <= hi:
                by_species.setdefault(r['species'], []).append((key, f, r))
    # hyperfine duplicates: keep the plain entry when both exist (as identify does)
    for s in [s for s in by_species if s.endswith('[hf]') and s[:-4].strip() in by_species]:
        del by_species[s]
    expect = [x for x in expect if x.strip()]
    used_expect = set()
    info = {}
    for species, trans in by_species.items():
        u = np.array([unit_intensity(r, T_ref, 'Emission') or 0. for _, _, r in trans])
        strongest = int(np.argmax(u)) if len(u) else 0
        strong = [f for (k, f, r), x in zip(trans, u) if u.max() > 0 and x >= strong_fraction*u.max()]
        p_c = _coverage(strong, half_kms, list(windows_obs.values()))
        hits = [x for x in expect if matches(x, species)]
        used_expect.update(hits)
        info[species] = dict(trans=trans, strongest=strongest, strong=strong, p_c=p_c, hits=hits, flags=kind_flags(species),
                             tier=tier_of(species, trans[0][2]), base=base_species(species))

    def chance(i):
        """Why a frequency match with this species says little, or '' if it is specific enough."""
        if i['p_c'] <= pmax:
            return ''
        if len(i['strong']) >= FOREST_LINES:
            return f"line forest: {len(i['strong'])} strong lines cover {i['p_c']:.0%} of the band"
        return f"low specificity: at these line widths its lines cover {i['p_c']:.0%} of the band"

    decided = {}
    # pass 1: species in their main form (ground state, main isotopologue)
    for species, i in info.items():
        if i['flags']:
            continue
        c = chance(i); listed = 'bright' if i['base'] in BRIGHT else 'common' if i['base'] in COMMON else None
        if i['hits']:
            decided[species] = ('test', 'on the expectation list'+(f' ({c}: needs consistent companions, not a single match)' if c else ''))
        elif listed == 'bright':
            decided[species] = ('test', 'bright in galaxies'+(f' ({c})' if c else ''))
        elif listed == 'common':
            decided[species] = ('reserve', f'common in galaxies, but {c}') if c else ('reserve', 'common in galaxies; put it on the expectation list if this source calls for it')
        elif i['tier'] == 'extragalactic':
            decided[species] = ('excluded', f'{c}: a frequency match is not informative') if c else ('reserve', 'detected in other galaxies, usually weak')
        else:
            decided[species] = ('excluded', f"{c}; {i['tier']}") if c else ('excluded', f"{i['tier']}; not expected in this source")
    tested_bases = {info[sp]['base'] for sp, d in decided.items() if d[0] == 'test'}
    # pass 2: excited states and rare isotopologues, judged against their parent
    for species, i in info.items():
        if species in decided:
            continue
        c = chance(i); parent = parent_of(i['base']); listed = i['base'] in BRIGHT or i['base'] in COMMON
        deuterated = re.search(r'D(?![a-z])', i['base']) is not None
        if i['hits']:
            decided[species] = ('test', 'on the expectation list'+(f' ({c}: needs consistent companions, not a single match)' if c else ''))
        elif 'vibrationally excited' in i['flags'] and i['base'] in tested_bases:
            decided[species] = ('excluded', f'vibrationally excited, {c}') if c and len(i['strong']) >= FOREST_LINES else \
                               ('reserve', f"vibrationally excited {i['base']}: test once the ground state is secure")
        elif 'rare isotopologue' in i['flags'] and i['base'] in BRIGHT:
            decided[species] = ('test', 'bright in galaxies'+(f' ({c})' if c else ''))
        elif 'rare isotopologue' in i['flags'] and listed and not deuterated:
            decided[species] = ('reserve', 'common in galaxies'+(f', but {c}' if c else ''))
        elif deuterated:
            decided[species] = ('excluded', 'deuterated: rarely detectable in galaxies; add it to the expectation list if the source calls for it')
        elif 'rare isotopologue' in i['flags'] and (parent in tested_bases or parent in BRIGHT):
            decided[species] = ('excluded', f'rare isotopologue of {parent}, {c}') if c and len(i['strong']) >= FOREST_LINES else \
                               ('reserve', f'rare isotopologue of {parent}: test once {parent} is secure')
        elif 'rare isotopologue' in i['flags']:
            decided[species] = ('excluded', f'rare isotopologue of {parent}, which is not tested')
        else:
            decided[species] = ('excluded', f"excited state of {i['base']}, which is not tested")
    rows = []
    for species, i in info.items():
        decision, reason = decided[species]
        k, f, r = i['trans'][i['strongest']]
        rows.append(dict(species=species, base=i['base'], decision=decision, reason=reason, tier=i['tier'], flags=i['flags'],
                         expected=bool(i['hits']), n_lines=len(i['trans']), n_strong=len(i['strong']), p_coincidence=round(i['p_c'], 3),
                         line_forest=i['p_c'] > pmax and len(i['strong']) >= FOREST_LINES, low_specificity=i['p_c'] > pmax,
                         aliases=vib_aliases(species),
                         strongest=dict(window=k, obs_GHz=f, rest_GHz=r['rest_GHz'], qns=r.get('qns'), Eu_K=r.get('Eu_K'))))
    order = {'test': 0, 'reserve': 1, 'excluded': 2}
    rows.sort(key=lambda x: (order[x['decision']], not x['expected'], ['bright', 'common', 'extragalactic', 'seen in space', 'not known in space'].index(x['tier']), x['species']))
    unmatched = [x for x in expect if x not in used_expect]
    return rows, unmatched


def run(run, expect=(), T_ref=150., pmax=.3, margin_kms=100.):
    catalogue = run.load('catalogue'); z = catalogue['z']
    inv = {c['key']: c for c in run.load('inventory')['cubes']}
    windows = {k: (inv[k]['freq_min_GHz'], inv[k]['freq_max_GHz']) for k in catalogue['windows'] if k in inv}
    feats = run.load('features', {}).get('features', [])
    widths = [f['width_kms'] for f in feats if f.get('significant')] or [f['width_kms'] for f in feats]
    half = max(margin_kms, .5*float(np.median(widths)) if widths else 0.)
    rows, unmatched = screen_species(catalogue, windows, z, half, expect, T_ref, pmax)
    counts = {d: sum(r['decision'] == d for r in rows) for d in ('test', 'reserve', 'excluded')}
    out = dict(z=z, T_ref=T_ref, pmax=pmax, half_width_kms=half, expect=list(expect), unmatched_expect=unmatched,
               bandwidth_GHz=sum(hi-lo for lo, hi in windows.values()), counts=counts, species=rows)
    run.save('plan', out)
    with open(run.product('species_screen.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['decision', 'species', 'reason', 'tier', 'flags', 'expected', 'lines_in_band', 'strong_lines', 'P_coincidence',
                    'strongest_rest_GHz', 'strongest_qns', 'strongest_Eu_K'])
        for r in rows:
            s = r['strongest']
            w.writerow([r['decision'], r['species'], r['reason'], r['tier'], ';'.join(r['flags']), r['expected'], r['n_lines'], r['n_strong'],
                        r['p_coincidence'], s['rest_GHz'], s['qns'], s['Eu_K']])
    _plots(run, catalogue, windows, z, rows, T_ref)
    run.log('plan', f"{counts['test']} species to test, {counts['reserve']} in reserve, {counts['excluded']} excluded",
            expect=list(expect), unmatched=unmatched, half_width_kms=half)
    return out


def _plots(run, catalogue, windows, z, rows, T_ref):
    """Aperture spectrum of each window with the strong lines of the species to be tested."""
    tested = [r['species'] for r in rows if r['decision'] == 'test']
    colour = {s: figs.COLOURS[i % len(figs.COLOURS)] for i, s in enumerate(tested)}
    for key, (lo, hi) in windows.items():
        p = run.path/'cache'/f'spectrum_{key}.npz'
        if not p.exists():
            continue
        d = np.load(p)
        lines = []
        for s in tested:
            trans = [r for r in catalogue['windows'][key]['rows'] if r['species'] == s and lo <= r['rest_GHz']/(1+z) <= hi]
            if not trans:
                continue
            u = np.array([unit_intensity(r, T_ref, 'Emission') or 0. for r in trans])
            for r, x in zip(trans, u):
                if u.max() > 0 and x >= .3*u.max():
                    lines.append((r['rest_GHz']/(1+z), s, f"{base_species(s)} {r.get('qns') or ''}"[:26], colour[s]))
        feats = [f for f in run.load('features', {}).get('features', []) if f['spw'] == key]
        figs.strategy(run.plot(f'strategy_{key}'), d['frequency'], d['spectrum'], d['noise'], feats, lines,
                      f"{key}: strong lines (≥30% of each species' strongest, LTE {T_ref:g} K) of the species to be tested", z=z)


def summary(out, show=40):
    lines = [f"Species screen (z = {out['z']}, T_ref = {out['T_ref']:g} K, match half-width {out['half_width_kms']:.0f} km/s, "
             f"line forest if P_c > {out['pmax']}): {out['counts']['test']} test, {out['counts']['reserve']} reserve, {out['counts']['excluded']} excluded"]
    if out['unmatched_expect']:
        lines.append(f"  ! expectation-list items with no transition in the observed windows: {', '.join(out['unmatched_expect'])}")
    for d in ('test', 'reserve'):
        sel = [r for r in out['species'] if r['decision'] == d]
        lines.append(f"  {d.upper()} ({len(sel)}):")
        for r in sel[:show]:
            name = r['species']+(f" [{r['aliases'][0]}]" if r.get('aliases') else '')
            lines.append(f"    {name:<34} P_c {r['p_coincidence']:.2f}  {r['n_lines']:>4} lines  {r['reason']}")
    ex = [r for r in out['species'] if r['decision'] == 'excluded']
    forests = [r['species'] for r in ex if r['line_forest']]
    lines.append(f"  EXCLUDED ({len(ex)}): {len(forests)} line forests (e.g. {', '.join(forests[:6])}), {len(ex)-len(forests)} others (not expected in galaxies, deuterated, or forms of untested species)")
    return '\n'.join(lines)
