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


def base_species(name):
    """'HCN v2=1' -> 'HCN'; 'CH3OH vt=0-2' -> 'CH3OH' (for priors)."""
    n = re.split(r'[,\s]+(?:v\d?|vt|v=|v\s*=)', name, maxsplit=1)[0]
    return n.strip().rstrip(',').strip()


def kind_flags(name):
    flags = []
    if re.search(r'v\d?\s*=\s*[1-9]|vt\s*=\s*[1-9]|v\d=1|\bv\s*=\s*1', name):
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


def candidates_for(feature, catalogue, spectra, z, margin_kms=100., sigma_kms=30., others=()):
    rows = catalogue['windows'][feature['spw']]['rows']
    pad = feature['freq_GHz']*margin_kms/C
    inside = [dict(r, expected_GHz=r['rest_GHz']/(1+z)) for r in rows
              if feature['freq_lo_GHz']-pad <= r['rest_GHz']/(1+z) <= feature['freq_hi_GHz']+pad]
    by_species = {}
    for r in inside:
        by_species.setdefault(r['species'], []).append(r)
    out = []
    for species, trans in by_species.items():
        # all transitions of this species in every window, with their own spectra
        everywhere = []
        for key, w in catalogue['windows'].items():
            if key not in spectra:
                continue
            f = spectra[key]['frequency']; lo, hi = min(f), max(f)
            for r in w['rows']:
                e = r['rest_GHz']/(1+z)
                if r['species'] == species and lo <= e <= hi:
                    everywhere.append((key, r))
        main = max(trans, key=lambda r: unit_intensity(r, 150., feature['sign']) or 0)
        offset = C*(main['expected_GHz']-feature['freq_GHz'])/main['expected_GHz']
        err = feature.get('centroid_error_kms')
        probabilities, _ = identification_probabilities([offset], sigma_kms, 0., max(margin_kms, 200.), err)
        in_feature = {(r['rest_GHz'], r.get('qns')) for r in trans}
        by_T = {}
        for T in T_GRID:
            units = [unit_intensity(r, T, feature['sign']) for r in trans]
            if any(u is None for u in units) or sum(units) <= 0:
                continue
            scale = feature['flux_Jykms']/sum(units)
            rows_T, fmax, detected, missing, testable, supported = [], 1., 0, 0, 0, 0
            for key, r in everywhere:
                if key == feature['spw'] and (((r['rest_GHz'], r.get('qns')) in in_feature)
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
                rows_T.append(dict(label=f"{r.get('qns', '')} {r['rest_GHz']:.4f}", window=key, rest_GHz=r['rest_GHz'], Eu_K=r.get('Eu_K'), overlap=overlap,
                                   predicted=pred, measured=meas, error=merr, status=status))
            by_T[T] = dict(fmax=fmax, companions=rows_T, detected=int(detected), missing=int(missing),
                           testable=int(testable), supported=int(supported))
        if by_T:
            # A temperature at which no companion is testable proves nothing: prefer temperatures
            # where the species explains at least half the feature AND its companions are seen.
            bestT = max(by_T, key=lambda T: (by_T[T]['fmax'] >= .5, by_T[T]['supported'], round(by_T[T]['fmax'], 2), -T))
            best = by_T[bestT]
        else:
            bestT, best = None, dict(fmax=None, companions=[], detected=0, missing=0, testable=0, supported=0)
        name = base_species(species)
        prior = prior_for(name, kind_flags(species), trans[0].get('seen_extragalactic'), trans[0].get('seen_in_space'))
        position = 1. if any(feature['freq_lo_GHz'] <= r['expected_GHz'] <= feature['freq_hi_GHz'] for r in trans) else .3
        if feature['channels'] <= 15 and err is not None:
            position *= max(probabilities[0], .05)
        plaus = best['fmax'] if best['fmax'] is not None else .5
        support = .6+.1*min(best['supported'], 4)          # companions predicted, testable and seen
        out.append(dict(species=species, base=name, flags=kind_flags(species), prior=prior,
                        transitions=[dict(rest_GHz=r['rest_GHz'], qns=r.get('qns'), Eu_K=r.get('Eu_K'), log10_A_s=r.get('log10_A_s'),
                                          expected_GHz=r['expected_GHz'], catalogue=r.get('catalogue'),
                                          inside=bool(feature['freq_lo_GHz'] <= r['expected_GHz'] <= feature['freq_hi_GHz'])) for r in trans],
                        main_rest_GHz=main['rest_GHz'], offset_kms=float(offset), id_probability=float(probabilities[0]),
                        best_T=bestT, fmax=best['fmax'], fmax_by_T={str(int(T)): v['fmax'] for T, v in by_T.items()},
                        companions_detected=best['detected'], companions_missing=best['missing'],
                        companions_testable=best['testable'], companions_supported=best['supported'],
                        companions=best['companions'], score=float(plaus*prior*position*support),
                        seen_extragalactic=trans[0].get('seen_extragalactic'), seen_in_space=trans[0].get('seen_in_space')))
    # '[hf]' entries repeat a species with hyperfine splitting: keep the plain entry when both exist.
    names = {c['species'] for c in out}
    out = [c for c in out if not (c['species'].endswith('[hf]') and c['species'][:-4].strip() in names)]
    return sorted(out, key=lambda c: -c['score'])


def run(run, margin_kms=100., include='reliable', sigma_kms=30., top=3):
    feats = run.load('features')['features']
    catalogue = run.load('catalogue'); z = catalogue['z']
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
        cands = candidates_for(f, catalogue, spectra, z, margin_kms, sigma_kms, feats)
        out[f['id']] = dict(feature=f, candidates=cands)
        sp = spectra[f['spw']]
        pad = max(f['freq_hi_GHz']-f['freq_lo_GHz'], 0.05)
        marks = [(t['expected_GHz'], f"{c['species']} {t['qns'] or ''}"[:28], figs.COLOURS[i % len(figs.COLOURS)])
                 for i, c in enumerate(cands[:6]) for t in c['transitions']]
        figs.spectrum(run.plot(f'feature_{f["id"]}'), sp['frequency'], sp['spectrum'], sp['noise'], [dict(f, significant=f['significant'])],
                      title=f"{f['id']}: candidates (top {min(6, len(cands))})", z=z, marks=marks,
                      xlim=(f['freq_lo_GHz']-pad, f['freq_hi_GHz']+pad))
        for i, c in enumerate(cands[:top]):
            if c['companions']:
                figs.companions(run.plot(f'companions_{f["id"]}_{i+1}'), c['companions'],
                                f"{f['id']} as {c['species']}: companions at T = {c['best_T']:g} K (f_max {c['fmax']:.2f})")
    run.save('candidates', dict(z=z, margin_kms=margin_kms, include=include, T_grid=T_GRID, features=out))
    run.log('identify', f'{len(out)} features tested', margin_kms=margin_kms, include=include)
    return out


def summary(out, top=5):
    lines = []
    for ident, d in out.items():
        f = d['feature']
        lines.append(f"{ident} ({f['sign'][:3]}, {f['freq_GHz']:.5f} GHz, {f['width_kms']:.0f} km/s, S/N {f['snr'] if f['snr'] is None else round(f['snr'], 1)}, "
                     f"rel {f['reliability'] if f['reliability'] is None else round(f['reliability'], 2)}): {len(d['candidates'])} candidate species")
        for c in d['candidates'][:top]:
            fm = 'n/a' if c['fmax'] is None else f"{c['fmax']:.2f}"
            lines.append(f"    {c['score']:.2f}  {c['species']:<28} dv {c['offset_kms']:+6.0f} km/s  f_max {fm} at T={c['best_T']}  "
                         f"companions: {c['companions_supported']} seen of {c['companions_testable']} testable, {c['companions_missing']} missing  prior {c['prior']:.2f} {' '.join(c['flags'])}")
    return '\n'.join(lines)
