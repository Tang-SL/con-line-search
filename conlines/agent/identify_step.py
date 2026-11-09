"""Step: candidate identifications for every feature, tested against the rest of the spectrum.

For each feature and each candidate species with a transition predicted inside it:

* position: is the predicted frequency inside the feature, and how far from its centroid?
* companions (LTE, optically thin): the species' other transitions in all observed windows are
  predicted relative to the feature's flux, at several excitation temperatures. Their measured
  fluxes (same velocity, same width) give f_max(T): the largest fraction of the feature this
  species can account for without a predicted companion line being missing. This is the
  "every predicted line must be there" criterion of Snyder et al. (2005), Belloche et al.
  (2008/2013) and ALCHEMI (Martín et al. 2021), turned into a number.
* priors: tiers from extragalactic surveys (bright / common / seen in galaxies / seen in space),
  lowered for excited states and rare isotopologues.
* blending: a companion whose window falls on another detected feature cannot support the
  candidate (its emission belongs to that feature), but it can still be 'missing'.
* support: companions that are predicted bright enough to test and are actually seen. A
  temperature at which no companion is testable proves nothing; the reported temperature is
  the one where the species explains at least half the feature with the most companions seen.

Ranking score (a heuristic for ordering, not a probability):
    f_max(T_best) x prior x position x (0.6 + 0.1 x min(companions seen, 4))

Emission: integrated flux density ∝ ν³ S µ² exp(-E_u/T) (optically thin, Rayleigh-Jeans, same
beam). Absorption: optical depth ∝ S µ² exp(-E_l/T) (1 - exp(-hν/kT)), an approximation that
ignores how the continuum varies with frequency. Optically thick or self-absorbed lines
(e.g. HCN, HCO+, CO in CONs) are fainter than this predicts: their companions are only a
loose check, and the agent must say so.
"""
import re
import numpy as np
from ..core import C
from ..reliability import identification_probabilities
from .catalog_step import strength
from . import figs

T_REF = 150.     # K: reported when the companion test does not depend on temperature
T_GRID = (10., 20., 40., 80., 150., 300., 500.)
# Priors from extragalactic line surveys (e.g. Martín 2011 review; Aladro et al. 2015; ALCHEMI):
# BRIGHT species make the strongest features in galaxies; COMMON ones are widely detected but
# usually much weaker. Excited states and rare isotopologues rarely dominate a feature alone.
BRIGHT = ('CO', '13CO', 'HCN', 'HCO+', 'HNC', 'CS', 'CN', 'CCH', 'C2H', 'H2O', 'CH3OH', 'H2CO', 'C18O', 'SiO', 'N2H+')
COMMON = ('C17O', 'H13CN', 'HC15N', 'H13CO+', 'HC18O+', 'HN13C', '13CS', 'C34S', '29SiO', '30SiO', 'SO', 'SO2', 'HC3N', 'CH3CN',
          'H2S', 'HNCO', 'CH3CCH', 'c-C3H2', 'OCS', 'NS', 'HOC+', 'H2CS', 'CH2NH', 'NH2CHO', 'HDO', 'NO', 'CH3SH', 'HC5N', 'C3H',
          'l-C3H', 'H2 13CO', 'NH2CN', 'CH3CHO', 'HCOOH', 'CH3OCH3', 'HCOOCH3', 'CH2CO', 'SiS', 'PN', 'HCS+', 'CO+', 'HCO')


def prior_for(name, flags, seen_extragalactic, seen_in_space):
    p = 1. if name in BRIGHT else .6 if name in COMMON else .4 if seen_extragalactic else .25 if seen_in_space else .1
    if 'vibrationally excited' in flags or 'rare isotopologue' in flags:
        p *= .6
    return p


# CDMS writes some vibrational states as quantum-number tuples, e.g. HC3N (v4,v5,v6,v7):
# 'HC3N, (0,0,0,1)' is v7=1 and 'HC3N, (0,1,0,0)/(0,0,0,3)' is v5=1/v7=3.
TUPLE_MODES = {'HC3N': (4, 5, 6, 7), 'H13CCCN': (4, 5, 6, 7), 'HC13CCN': (4, 5, 6, 7), 'HCC13CN': (4, 5, 6, 7), 'HC3N-15': (4, 5, 6, 7)}
TUPLE = re.compile(r'\((\d+(?:,\d+)+)\)')


def vib_aliases(name):
    """Readable names of a tuple-notation state: 'HC3N, (0,0,0,1)' -> ['HC3N v7=1']."""
    base = re.split(r',\s*\(', name, maxsplit=1)[0].strip()
    modes = TUPLE_MODES.get(base)
    out = []
    for t in TUPLE.findall(name):
        vals = [int(x) for x in t.split(',')]
        if not modes or len(vals) != len(modes):
            continue
        parts = [f'v{m}={v}' for m, v in zip(modes, vals) if v]
        out.append(base+' '+(' '.join(parts) if parts else 'v=0'))
    return out


def base_species(name):
    """'HCN v2=1' -> 'HCN'; 'CH3OH vt=0-2' -> 'CH3OH'; 'HC3N, (0,0,0,1)' -> 'HC3N' (for priors)."""
    n = re.split(r'[,\s]+(?:v\d?|vt|v=|v\s*=|\(\d)', name, maxsplit=1)[0]
    return n.strip().rstrip(',').strip()


def kind_flags(name):
    flags = []
    if any(int(x) for t in TUPLE.findall(name) for x in t.split(',')) or \
            re.search(r'\bv\d*\s*=\s*[1-9]|vt\s*=\s*[1-9]', name):
        flags.append('vibrationally excited')
    if re.search(r'13C|15N|18O|17O|34S|33S|29Si|30Si|\bD\b|D(?![a-z])', name):
        flags.append('rare isotopologue')
    return flags


def unit_intensity(row, T, sign):
    s = row.get('S_mu2') or strength(row)
    if not s:
        return None
    nu = row['rest_GHz']
    if sign == 'Emission':
        return nu**3*s*np.exp(-(row.get('Eu_K') or 0.)/T)
    el = row.get('El_K')
    if el is None:
        el = (row.get('Eu_K') or 0.)-0.048*nu
    return s*np.exp(-el/T)*(1-np.exp(-0.048*nu/T))


def measure(freq, spec, noise, centre, width_kms, sign):
    """Signed integrated flux over |v - v_centre| <= width/2 (Jy/beam km/s) and its 1σ error."""
    dv = C*abs(freq[1]-freq[0])/centre
    v = C*(centre-freq)/centre
    sel = np.abs(v) <= width_kms/2
    if sel.sum() < 1:
        return None, None
    s = 1 if sign == 'Emission' else -1
    return float(s*np.sum(spec[sel])*dv), float(np.sqrt(np.sum(noise[sel]**2))*dv)


HF_TOL_KMS = 5.        # transitions of one species closer than this are one line (hyperfine structure)
HF_LABEL = re.compile(r'\bF\d?\s*=')


def merge_hyperfine(rows, tol_kms=HF_TOL_KMS):
    """One entry per resolvable line of a species.

    Catalogues list hyperfine components separately, often next to an entry for the whole
    line (e.g. HCN v2=1 J=3-2 l=1e and its F=4-3, 3-2, 2-1 components within 1 MHz). Counted
    separately they look like extra companion lines and double the predicted flux. Within each
    group closer than tol_kms: keep the entry without an F quantum number if there is one,
    otherwise sum the components (strength-weighted frequency, lowest E_u)."""
    rows = sorted(rows, key=lambda r: r['rest_GHz'])
    groups, out = [], []
    for r in rows:
        if groups and (r['rest_GHz']-groups[-1][-1]['rest_GHz'])/r['rest_GHz']*C <= tol_kms:
            groups[-1].append(r)
        else:
            groups.append([r])
    for g in groups:
        if len(g) == 1:
            out.append(g[0]); continue
        total = [r for r in g if not HF_LABEL.search(str(r.get('qns') or ''))]
        if total:
            out.append(max(total, key=lambda r: r.get('S_mu2') or strength(r) or 0)); continue
        w = np.array([(r.get('S_mu2') or strength(r) or 0.) for r in g])
        merged = dict(g[int(np.argmax(w))])
        if w.sum() > 0:
            merged['rest_GHz'] = float(np.sum(w*[r['rest_GHz'] for r in g])/w.sum())
            merged['S_mu2'] = float(w.sum())
        merged['Eu_K'] = min((r.get('Eu_K') or 0.) for r in g)
        merged['qns'] = HF_LABEL.split(str(g[0].get('qns') or ''))[0].rstrip(' ,')+f' ({len(g)} hf components)'
        out.append(merged)
    return out


def verdict(c):
    """Plain-language outcome of the companion test for one candidate."""
    n = len(c.get('companions') or [])
    if c['companions_missing'] and (c.get('fmax') or 0) >= .1:
        return 'partial', (f"it cannot make the whole feature: its other lines limit it to at most {c['fmax']:.0%} "
                           "(possible as a blend partner, to be tested with a blend fit)")
    if c['companions_missing']:
        return 'rejected', f"{c['companions_missing']} of its other lines should be clearly visible but are missing"
    if c['companions_supported']:
        return 'supported', f"{c['companions_supported']} of its other lines are seen at the predicted strength"
    if c['companions_testable']:
        return 'inconclusive', 'its testable lines fall on other features (blended), so they neither support nor reject it'
    if n:
        return 'untested', 'its other lines in these windows are predicted below the noise'
    return 'untested', 'it has no other line in the observed windows'


def candidates_for(feature, catalogue, spectra, z, margin_kms=100., sigma_kms=30., others=(), allowed=None):
    """Candidate species for one feature; `allowed` (a set of species names) restricts the list."""
    rows = catalogue['windows'][feature['spw']]['rows']
    pad = feature['freq_GHz']*margin_kms/C
    inside = [dict(r, expected_GHz=r['rest_GHz']/(1+z)) for r in rows
              if feature['freq_lo_GHz']-pad <= r['rest_GHz']/(1+z) <= feature['freq_hi_GHz']+pad]
    by_species = {}
    for r in inside:
        if allowed is None or r['species'] in allowed:
            by_species.setdefault(r['species'], []).append(r)
    by_species = {k: merge_hyperfine(v) for k, v in by_species.items()}
    out = []
    for species, trans in by_species.items():
        # all transitions of this species in every window, with their own spectra
        everywhere = []
        for key, w in catalogue['windows'].items():
            if key not in spectra:
                continue
            f = spectra[key]['frequency']; lo, hi = min(f), max(f)
            mine = [r for r in w['rows'] if r['species'] == species and lo <= r['rest_GHz']/(1+z) <= hi]
            everywhere += [(key, r) for r in merge_hyperfine(mine)]
        main = max(trans, key=lambda r: unit_intensity(r, 150., feature['sign']) or 0)
        offset = C*(main['expected_GHz']-feature['freq_GHz'])/main['expected_GHz']
        err = feature.get('centroid_error_kms')
        probabilities, _ = identification_probabilities([offset], sigma_kms, 0., max(margin_kms, 200.), err)
        in_feature = {round(r['rest_GHz'], 4) for r in trans}
        by_T = {}
        for T in T_GRID:
            units = [unit_intensity(r, T, feature['sign']) for r in trans]
            if any(u is None for u in units) or sum(units) <= 0:
                continue
            scale = feature['flux_Jykms']/sum(units)
            rows_T, fmax, detected, missing, testable, supported = [], 1., 0, 0, 0, 0
            for key, r in everywhere:
                if key == feature['spw'] and ((round(r['rest_GHz'], 4) in in_feature)
                                              or feature['freq_lo_GHz'] <= r['rest_GHz']/(1+z)*(1-offset/C) <= feature['freq_hi_GHz']):
                    continue          # part of this feature: not an independent check
                u = unit_intensity(r, T, feature['sign'])
                if u is None:
                    continue
                pred = scale*u
                centre = r['rest_GHz']/(1+z)*(1-offset/C)          # at the feature's velocity
                sp = spectra[key]
                meas, merr = measure(sp['frequency'], sp['spectrum'], sp['noise'], centre, feature['width_kms'], feature['sign'])
                if meas is None:
                    continue
                if pred > 3*merr:
                    fmax = min(fmax, max(0., (meas+3*merr)/pred))
                half = feature['width_kms']/2*centre/C
                overlap = next((o['id'] for o in others if o['spw'] == key and o['id'] != feature['id']
                                and o['freq_lo_GHz'] < centre+half and o['freq_hi_GHz'] > centre-half), None)
                status = ('detected' if meas > 3*merr else 'not detected') if pred <= 3*merr else \
                         ('consistent' if meas >= pred-3*merr-.3*pred else 'missing')
                if overlap and status in ('consistent', 'detected'):
                    status = f'blended with {overlap}'     # emission there belongs to another feature: no support
                detected += meas > 3*merr and pred > merr and not overlap
                missing += status == 'missing'
                testable += pred > 3*merr
                supported += status == 'consistent' and meas > 3*merr
                rows_T.append(dict(label=f"{r.get('qns', '')} {r['rest_GHz']:.4f}", qns=r.get('qns'), window=key, rest_GHz=r['rest_GHz'], Eu_K=r.get('Eu_K'),
                                   centre_GHz=centre, overlap=overlap,
                                   predicted=pred, measured=meas, error=merr, status=status))
            by_T[T] = dict(fmax=fmax, companions=rows_T, detected=int(detected), missing=int(missing),
                           testable=int(testable), supported=int(supported))
        if by_T:
            # A temperature at which no companion is testable proves nothing: prefer temperatures
            # where the species explains at least half the feature AND its companions are seen.
            # When the temperature makes no difference, report one near 150 K, not the edge of the grid.
            bestT = max(by_T, key=lambda T: (by_T[T]['fmax'] >= .5, by_T[T]['supported'], round(by_T[T]['fmax'], 2), -abs(np.log(T/T_REF))))
            best = by_T[bestT]
            good = [T for T in by_T if by_T[T]['fmax'] >= max(.5, best['fmax']-.05)] or [bestT]
            T_range = (min(good), max(good))
        else:
            bestT, best, T_range = None, dict(fmax=None, companions=[], detected=0, missing=0, testable=0, supported=0), None
        name = base_species(species)
        prior = prior_for(name, kind_flags(species), trans[0].get('seen_extragalactic'), trans[0].get('seen_in_space'))
        inside = any(feature['freq_lo_GHz'] <= r['expected_GHz'] <= feature['freq_hi_GHz'] for r in trans)
        position = 1. if inside else .3
        if feature['channels'] <= 15 and err is not None:
            position *= max(probabilities[0], .05)
        plaus = best['fmax'] if best['fmax'] is not None else .5
        support = .6+.1*min(best['supported'], 4)          # companions predicted, testable and seen
        out.append(dict(species=species, base=name, flags=kind_flags(species), prior=prior, position=position, inside=inside,
                        transitions=[dict(rest_GHz=r['rest_GHz'], qns=r.get('qns'), Eu_K=r.get('Eu_K'), log10_A_s=r.get('log10_A_s'),
                                          expected_GHz=r['expected_GHz'], catalogue=r.get('catalogue'),
                                          inside=bool(feature['freq_lo_GHz'] <= r['expected_GHz'] <= feature['freq_hi_GHz'])) for r in trans],
                        main_rest_GHz=main['rest_GHz'], offset_kms=float(offset), id_probability=float(probabilities[0]),
                        best_T=bestT, T_range=T_range, fmax=best['fmax'], fmax_by_T={str(int(T)): v['fmax'] for T, v in by_T.items()},
                        companions_detected=best['detected'], companions_missing=best['missing'],
                        companions_testable=best['testable'], companions_supported=best['supported'],
                        companions=best['companions'], score=float(plaus*prior*position*support),
                        seen_extragalactic=trans[0].get('seen_extragalactic'), seen_in_space=trans[0].get('seen_in_space')))
    for c in out:
        c['verdict'], c['verdict_reason'] = verdict(c)
    # '[hf]' entries repeat a species with hyperfine splitting: keep the plain entry when both exist.
    names = {c['species'] for c in out}
    out = [c for c in out if not (c['species'].endswith('[hf]') and c['species'][:-4].strip() in names)]
    return sorted(out, key=lambda c: -c['score'])


def coincident(feature, catalogue, z, species, margin_kms=100.):
    """Species from a list with a catalogue line inside the feature (± margin)."""
    pad = feature['freq_GHz']*margin_kms/C
    lo, hi = feature['freq_lo_GHz']-pad, feature['freq_hi_GHz']+pad
    return sorted({r['species'] for r in catalogue['windows'][feature['spw']]['rows'] if r['species'] in species and lo <= r['rest_GHz']/(1+z) <= hi})


def run(run, margin_kms=100., include='reliable', sigma_kms=30., top=3, species='plan'):
    """species: 'plan' = the species screen of the plan step (tested species for every feature; reserve
    species only for features that no tested species explains); 'all' = every catalogue species."""
    feats = run.load('features')['features']
    catalogue = run.load('catalogue'); z = catalogue['z']
    plan = run.load('plan') if species == 'plan' and run.has('plan') else None
    if species == 'plan' and plan is None:
        print('  ! no species screen (run `plan` first): testing every catalogue species; candidate lists will be long')
    groups = {d: {r['species'] for r in plan['species'] if r['decision'] == d} for d in ('test', 'reserve', 'excluded')} if plan else None
    forests = {r['species'] for r in plan['species'] if r['line_forest']} if plan else set()
    lowspec = {r['species'] for r in plan['species'] if r.get('low_specificity') and not r['line_forest']} if plan else set()
    spectra = {}
    for key in catalogue['windows']:
        p = run.path/'cache'/f'spectrum_{key}.npz'
        if p.exists():
            d = np.load(p)
            spectra[key] = dict(frequency=d['frequency'], spectrum=d['spectrum'], noise=d['noise'])
    chosen = [f for f in feats if include == 'all' or (include == 'reliable' and (f['significant'] or (f['reliable'] and (f['snr'] or 0) >= 3)))
              or (include == 'significant' and f['significant'] and f['reliable'])]
    out = {}
    for f in chosen:
        stage = 'all species'
        if groups is None:
            cands = candidates_for(f, catalogue, spectra, z, margin_kms, sigma_kms, feats)
        else:
            cands = candidates_for(f, catalogue, spectra, z, margin_kms, sigma_kms, feats, allowed=groups['test'])
            stage = 'tested species'
            # widen to the reserve only for a significant feature that no tested species can explain
            # (a line inside the feature that accounts for at least half of it)
            if f.get('significant') and not any(c['inside'] and (c['fmax'] is None or c['fmax'] >= .5) for c in cands):
                extra = candidates_for(f, catalogue, spectra, z, margin_kms, sigma_kms, feats, allowed=groups['reserve'])
                for c in extra:
                    c['stage'] = 'reserve'
                cands = sorted(cands+extra, key=lambda c: -c['score'])
                stage = 'tested + reserve species (no tested species explains it)'
        for c in cands:
            c.setdefault('stage', 'tested' if groups else 'all')
            c['line_forest'] = c['species'] in forests
            c['low_specificity'] = c['species'] in lowspec
        for c in cands[top:]:
            c['companions'] = []          # keep the file small: companion tables only for the top candidates
        excluded = coincident(f, catalogue, z, groups['excluded'], margin_kms) if groups else []
        out[f['id']] = dict(feature=f, candidates=cands, stage=stage, excluded_coincident=excluded)
        sp = spectra[f['spw']]
        pad = max(f['freq_hi_GHz']-f['freq_lo_GHz'], 0.05)
        lines = [(t['expected_GHz'], c['species'], f"{c['base']} {t['qns'] or ''}"[:26], figs.COLOURS[i % len(figs.COLOURS)])
                 for i, c in enumerate(cands[:6]) for t in c['transitions'][:4]]
        figs.strategy(run.plot(f'feature_{f["id"]}'), sp['frequency'], sp['spectrum'], sp['noise'], [dict(f, significant=f['significant'])], lines,
                      f"{f['id']}: lines of the leading candidates (top {min(6, len(cands))} of {len(cands)}), at the redshift", z=z,
                      xlim=(f['freq_lo_GHz']-pad, f['freq_hi_GHz']+pad))
        for i, c in enumerate(cands[:top]):
            figs.companions(run.plot(f'companions_{f["id"]}_{i+1}'), f, c, spectra)
    run.save('candidates', dict(z=z, margin_kms=margin_kms, include=include, species=species if groups else 'all', T_grid=T_GRID, features=out))
    run.log('identify', f'{len(out)} features tested', margin_kms=margin_kms, include=include, species=species if groups else 'all')
    return out


def summary(out, top=5):
    lines = []
    for ident, d in out.items():
        f = d['feature']
        lines.append(f"{ident} ({f['sign'][:3]}, {f['freq_GHz']:.5f} GHz, {f['width_kms']:.0f} km/s, S/N {f['snr'] if f['snr'] is None else round(f['snr'], 1)}, "
                     f"rel {f['reliability'] if f['reliability'] is None else round(f['reliability'], 2)}): {len(d['candidates'])} candidate species "
                     f"[{d.get('stage', 'all species')}]")
        for c in d['candidates'][:top]:
            fm = 'n/a' if c['fmax'] is None else f"{c['fmax']:.2f}"
            lines.append(f"    {c['score']:.2f}  {c['species']:<28} dv {c['offset_kms']:+6.0f} km/s  f_max {fm} at T={c['best_T']}  "
                         f"{c.get('verdict', ''):<12} companions: {c['companions_supported']} seen of {c['companions_testable']} testable, {c['companions_missing']} missing  prior {c['prior']:.2f} {' '.join(c['flags'])}"
                         + ('  [reserve]' if c.get('stage') == 'reserve' else '') + ('  [line forest]' if c.get('line_forest') else '') + ('  [low specificity]' if c.get('low_specificity') else ''))
        if d.get('excluded_coincident'):
            ex = d['excluded_coincident']
            lines.append(f"    (also has lines here, excluded by the species screen: {len(ex)}, e.g. {', '.join(ex[:5])})")
    return '\n'.join(lines)
