"""ADMIT-style line identification for ESPADA features.

A port of the identification rules in ADMIT's LineID_AT (ALMA Data Mining
Toolkit, University of Maryland, MIT licence; see conlines/data/ADMIT_LICENSE).
ADMIT itself requires CASA and its own pipeline objects, so the decision logic is
reproduced here on the features con-line-search has already detected:

1. Tier 1: ADMIT's list of bright, common transitions (admit_tier1_lines.list,
   copied unchanged) is matched first. For each Tier 1 transition in the
   coverage, the closest feature within +/-40 km/s (+/-300 km/s when
   |v_sys| > 150 km/s) is assigned to it. Its hyperfine components are reported
   as blends. Catalogue lines inside a claimed Tier 1 window are excluded from
   later searches.
2. Other features: catalogue transitions whose rest frequency lies inside the
   feature's own frequency extent are the candidates, after ADMIT's search
   filters (NRAO-recommended entries with a line strength; no exotic atoms
   unless allowed; recombination lines limited by `recomb_level`).
3. narrowpossibles: drop rare isotopologues of heavier molecules, take the
   lowest upper-state energy, prefer a much lighter species with a strong line,
   keep same-species entries at the same energy, and choose the highest line
   strength; the others are reported as blends.
4. No candidate: unidentified ("U").

Not ported: ADMIT's own segment/peak finding (ESPADA replaces it), its search
for pairs of equally spaced peaks from rotating sources, Gaussian refits, and
Splatalogue's 'atmospheric/potential/probable' species flags (not returned by
the live query).
"""
from pathlib import Path
import re
from .core import C

DATA = Path(__file__).resolve().parent/'data'/'admit_tier1_lines.list'
EXOTICS = ['Al', 'Cl', 'Mg', 'Mn', 'F', 'Li', 'Na', 'K', 'Ti']
RARE_ISOTOPES = ['D', '13C', '17O', '18O', '15N', '6Li', '7Li', '29Si', '30Si', '33S', '34S',
                 '46Ti', '47Ti', '49Ti', '50Ti', '54Fe', '57Fe', '58Fe']
ISOTOPE_MASS = {'13C': 13, '12C': 12, '16O': 16, '17O': 17, '18O': 18, '14N': 14, '15N': 15,
                '6Li': 6, '7Li': 7, '28Si': 28, '29Si': 29, '30Si': 30, '32S': 32, '33S': 33,
                '34S': 34, '36S': 36, '46Ti': 46, '47Ti': 47, '48Ti': 48, '49Ti': 49, '50Ti': 50,
                '54Fe': 54, '56Fe': 56, '57Fe': 57, '58Fe': 58, '24Mg': 24, '25Mg': 25, '26Mg': 26,
                '37Cl': 37, '35Cl': 35}
ELEMENT_MASS = {'H': 1, 'D': 2, 'He': 4, 'Li': 7, 'Be': 9, 'B': 11, 'C': 12, 'N': 14, 'O': 16,
                'F': 19, 'Ne': 20, 'Na': 23, 'Mg': 24, 'Al': 27, 'Si': 28, 'P': 31, 'S': 32,
                'Cl': 35, 'Ar': 40, 'K': 39, 'Ca': 40, 'Ti': 48, 'V': 51, 'Cr': 52, 'Mn': 55,
                'Fe': 56, 'Co': 59, 'Ni': 59, 'Cu': 64, 'Zn': 65}


def load_tier1(path=DATA):
    """Primary transitions and their hyperfine components (rest GHz)."""
    primary, hyperfine = [], {}
    section = 'primary'
    for line in Path(path).read_text().splitlines():
        if line.startswith('# Hyperfine'):
            section = 'hyperfine'
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split()
        freq, species = float(parts[0]), parts[1]
        if section == 'primary':
            primary.append(dict(species=species, rest_GHz=freq))
        else:
            group = re.search(r'#hfnum\s+(\d+)', line)
            hyperfine.setdefault((species, group.group(1) if group else freq), []).append(freq)
    # Attach each hyperfine group to the nearest primary transition of that species.
    for (species, _), freqs in hyperfine.items():
        centre = sum(freqs)/len(freqs)
        own = [p for p in primary if p['species'] == species]
        if own:
            nearest = min(own, key=lambda p: abs(p['rest_GHz']-centre))
            if abs(nearest['rest_GHz']-centre) < .5:
                nearest.setdefault('hyperfine_GHz', []).extend(sorted(freqs))
    return primary


def plain_formula(species):
    """Strip state labels and isomer prefixes: 'g-CH3CH2OH, v=0' -> 'CH3CH2OH'."""
    form = re.split(r',|\s', species.strip())[0]
    form = re.sub(r'v[t0-9]*\s*=.*$', '', form)
    pos = form.find('-')
    if pos != -1 and not any(-1 < form.find(a) < pos for a in 'CNOH'):
        form = form[pos+1:]
    return form


def isotope_count(species):
    form = plain_formula(species)
    count = 0
    for token in sorted(RARE_ISOTOPES, key=len, reverse=True):
        if token == 'D':
            n = len(re.findall(r'D(?![a-z])', form))
        else:
            n = form.count(token)
            form = form.replace(token, '')
        count += n
    return count


def rough_mass(species):
    """Integer molecular mass from a catalogue formula, as ADMIT's getmass."""
    form = plain_formula(species).replace('+', '').replace('-', '')
    mass = 0
    for token in sorted(ISOTOPE_MASS, key=len, reverse=True):
        while token in form:
            i = form.index(token)
            j = i+len(token)
            mult = re.match(r'\d+', form[j:])
            n = int(mult.group()) if mult else 1
            mass += n*ISOTOPE_MASS[token]
            form = form[:i]+form[j+(len(mult.group()) if mult else 0):]
    # Expand parenthesised groups such as (CH3)2CO.
    while '(' in form:
        form = re.sub(r'\(([^()]*)\)(\d*)', lambda m: m.group(1)*int(m.group(2) or 1), form)
    for symbol, count in re.findall(r'([A-Z][a-z]?)(\d*)', form):
        mass += ELEMENT_MASS.get(symbol, 0)*int(count or 1)
    return mass


def is_exotic(species):
    form = plain_formula(species)
    return any(e in form for e in EXOTICS)


def recombination_allowed(row, level):
    name = (row.get('chemical_name') or '')+' '+row['species']
    if 'recombination' not in name.lower():
        return True
    if level == 'off':
        return False
    if level == 'deep':
        return True
    species = row['species']
    return species.startswith(('H', 'He')) and any(s in species.lower() for s in ['α', 'β', 'alpha', 'beta'])


def search_filter(row, allow_exotics=False, recomb_level='shallow'):
    """ADMIT's generatepossibles screen applied to one catalogue row."""
    if str(row.get('recommended', '1')) not in ('1', 'True', 'true'):
        return 'not NRAO-recommended'
    if row.get('line_strength') in (None, ''):
        return 'no line strength'
    if not allow_exotics and is_exotic(row['species']):
        return 'exotic atom'
    if not recombination_allowed(row, recomb_level):
        return 'recombination line excluded'
    return None


def close_in_energy(e1, e2):
    return abs(e2-e1) < 1. and abs((e2+1e-4)/(e1+1e-4)-1) < .01


def narrow_possibles(possible):
    """ADMIT LineID_AT.narrowpossibles. Returns (winner, blends, notes)."""
    p = [dict(r, admit_mass=rough_mass(r['species']), admit_isotopes=isotope_count(r['species'])) for r in possible]
    notes = []
    for j in range(len(p)-1, -1, -1):
        r = p[j]
        if len(p) > 1 and ((r['admit_isotopes'] > 0 and r['admit_mass'] > 50) or
                           (r['admit_isotopes'] > 1 and r['admit_mass'] > 31) or
                           (r['admit_isotopes'] > 2 and r['admit_mass'] > 14)):
            notes.append(f"dropped isotopologue {r['species']}")
            del p[j]
    energy = lambda r: r['Eu_K'] if r.get('Eu_K') is not None else float('inf')
    strength = lambda r: r.get('line_strength') or 0.
    le = min(range(len(p)), key=lambda j: energy(p[j]))
    lm = min(range(len(p)), key=lambda j: p[j]['admit_mass'])
    rule = 'lowest upper-state energy'
    if p[le]['admit_mass'] > 1.6*p[lm]['admit_mass'] and strength(p[lm]) > 1. and '+' not in p[lm]['species']:
        le = lm
        rule = 'much lighter species with a strong line'
    chosen = p[le]
    keep = [r for j, r in enumerate(p) if j == le or (
        r['species'] == chosen['species'] and close_in_energy(energy(chosen), energy(r))
        and not (energy(r) == energy(chosen) and r.get('qns') == chosen.get('qns')))]
    winner = max(keep, key=strength)
    if winner is not chosen:
        rule += '; highest line strength within that species/energy'
    blends = [r for r in keep if r is not winner]
    notes.insert(0, rule)
    return winner, blends, notes


def identify(features, rows, z, frequencies_by_spw, systemic_kms=None, tier1_width_kms=None,
             allow_exotics=False, recomb_level='shallow', tier1=None):
    """Identify features with ADMIT's rules.

    features: list of dicts with id, spw, sign, centre_GHz (observed), lo_GHz, hi_GHz.
    rows: catalogue rows (rest frequencies) covering the loaded windows.
    frequencies_by_spw: spw -> (min, max) observed GHz, to select Tier 1 lines in coverage.
    Returns ident -> dict(status, rule, winner, candidates, blends, notes).
    """
    tier1 = tier1 if tier1 is not None else load_tier1()
    if systemic_kms is None:
        systemic_kms = C*z/(1+z)
    width = tier1_width_kms if tier1_width_kms else (300. if abs(systemic_kms) > 150 else 40.)
    rest = lambda f: f*(1+z)
    results = {}
    claimed_windows = []
    free = {f['id']: f for f in features}
    # 1. Tier 1, line by line.
    for t in sorted(tier1, key=lambda t: t['rest_GHz']):
        observed = t['rest_GHz']/(1+z)
        if not any(lo <= observed <= hi for lo, hi in frequencies_by_spw.values()):
            continue
        half = width/C*t['rest_GHz']
        near = [f for f in free.values() if abs(rest(f['centre_GHz'])-t['rest_GHz']) < half]
        if not near:
            continue
        best = min(near, key=lambda f: abs(rest(f['centre_GHz'])-t['rest_GHz']))
        claimed_windows.append((t['rest_GHz']-half, t['rest_GHz']+half))
        winner = enrich_tier1(t, rows)
        hf = [dict(winner, rest_GHz=h, qns=f"hyperfine component {i+1}", admit_note='Tier 1 hyperfine component')
              for i, h in enumerate(t.get('hyperfine_GHz', [])) if abs(h-t['rest_GHz']) > 1e-6]
        results[best['id']] = dict(status='assigned', rule=f'Tier 1 within ±{width:g} km/s', winner=winner,
                                   candidates=[winner]+hf, blends=hf, notes=[f"ADMIT Tier 1 list: {t['species']}"])
        del free[best['id']]
    # 2-4. Catalogue search for the remaining features.
    usable = []
    for r in rows:
        reason = search_filter(r, allow_exotics, recomb_level)
        if reason is None and not any(lo < r['rest_GHz'] < hi for lo, hi in claimed_windows):
            usable.append(r)
    for f in free.values():
        lo, hi = sorted([rest(f['lo_GHz']), rest(f['hi_GHz'])])
        possible = [r for r in usable if lo <= r['rest_GHz'] <= hi]
        if not possible:
            results[f['id']] = dict(status='unidentified', rule='no catalogue line inside the feature', winner=None,
                                    candidates=[], blends=[], notes=['ADMIT: U (no possibilities)'])
            continue
        winner, blends, notes = narrow_possibles(possible)
        key = lambda r: (r['species'], r['rest_GHz'], r.get('qns', ''))
        used = {key(winner), *(key(b) for b in blends)}
        ordered = [winner]+blends+[r for r in possible if key(r) not in used]
        results[f['id']] = dict(status='assigned', rule=notes[0], winner=winner, candidates=ordered,
                                blends=blends, notes=notes)
    return results


def enrich_tier1(t, rows):
    """Attach catalogue metadata (QNs, Eu, A) to a Tier 1 frequency when available."""
    key = re.sub(r'[^A-Za-z0-9+]', '', t['species']).lower()
    best = None
    for r in rows:
        if abs(r['rest_GHz']-t['rest_GHz']) < 1e-4 and key in re.sub(r'[^A-Za-z0-9+]', '', plain_formula(r['species'])).lower():
            best = r
            break
    row = dict(best) if best else dict(species=t['species'], rest_GHz=t['rest_GHz'], qns='', Eu_K=None,
                                        log10_A_s=None, catalogue='ADMIT Tier 1 list')
    row['admit_note'] = 'Tier 1'
    return row
