import json
from pathlib import Path
import pytest
from conlines import admit_lineid as A
from conlines.core import parse_splatalogue, C

ROOT = Path(__file__).resolve().parents[1]
Z = 3880/C


def row(species, rest, eu, ls, qns='1-0'):
    return dict(species=species, rest_GHz=rest, Eu_K=eu, line_strength=ls, qns=qns, recommended='1')


@pytest.mark.parametrize('species,mass,iso', [('CH3OH vt = 0-2', 32, 0), ('30SiO v = 0', 46, 1), ('(CH3)2CO', 58, 0),
                                              ('13CH3OH', 33, 1), ('HC3N', 51, 0), ('g-CH3CH2OH', 46, 0),
                                              ('DCN', 28, 1), ('H2C34S', 48, 1), ('HDO', 19, 1)])
def test_mass_and_isotope_count(species, mass, iso):
    assert A.rough_mass(species) == mass and A.isotope_count(species) == iso


def test_narrowpossibles_rules():
    # Heavy isotopologue dropped; then lowest upper energy wins.
    w, blends, notes = A.narrow_possibles([row('CH3OCHO', 100., 50., 5.), row('13CH3OCHO', 100.001, 10., 5.), row('HNCO', 100.002, 30., 5.)])
    assert w['species'] == 'HNCO' and any('13CH3OCHO' in n for n in notes)
    # A much lighter species with a strong line overrides the lowest energy...
    w, _, notes = A.narrow_possibles([row('CH3CH2CN', 100., 10., 5.), row('CO', 100.001, 30., 2.)])
    assert w['species'] == 'CO' and 'lighter' in notes[0]
    # ...but not when it is an ion or its line is weak.
    assert A.narrow_possibles([row('CH3CH2CN', 100., 10., 5.), row('HCO+', 100.001, 30., 2.)])[0]['species'] == 'CH3CH2CN'
    assert A.narrow_possibles([row('CH3CH2CN', 100., 10., 5.), row('CO', 100.001, 30., .5)])[0]['species'] == 'CH3CH2CN'
    # Same species at (almost) the same energy: highest line strength wins, the other is a blend.
    w, blends, _ = A.narrow_possibles([row('HCOOH', 100., 20., 1., 'a'), row('HCOOH', 100.001, 20.05, 9., 'b')])
    assert w['qns'] == 'b' and [b['qns'] for b in blends] == ['a']


def test_tier1_list_and_hyperfine():
    tier1 = A.load_tier1()
    cch = [t for t in tier1 if t['species'] == 'CCH' and abs(t['rest_GHz']-174.6632) < 1e-4][0]
    assert len(cch['hyperfine_GHz']) >= 3
    assert any(t['species'] == 'SiO' and abs(t['rest_GHz']-173.68831) < 1e-5 for t in tier1)


def test_tier1_claims_before_catalogue_and_uses_300_kms_for_galaxies():
    sio = 173.68831
    f = sio/(1+Z)*(1+150/C)          # 150 km/s from SiO 4-3
    feats = [dict(id='A', spw='s', sign='Absorption', centre_GHz=f, lo_GHz=f-.01, hi_GHz=f+.01)]
    rows = [row('CH3OH', f*(1+Z), 20., 5.)]
    r = A.identify(feats, rows, Z, {'s': (f-1, f+1)})['A']
    assert r['winner']['species'] == 'SiO' and 'Tier 1' in r['rule']
    # A Galactic source (|v| < 150 km/s) uses +/-40 km/s, so the catalogue line wins instead.
    z0 = 5/C
    g = sio/(1+z0)*(1+150/C)
    feats = [dict(id='A', spw='s', sign='Emission', centre_GHz=g, lo_GHz=g-.002, hi_GHz=g+.002)]
    r = A.identify(feats, [row('CH3OH', g*(1+z0), 20., 5.)], z0, {'s': (g-1, g+1)})['A']
    assert r['winner']['species'] == 'CH3OH'


def test_ic860_spw0_absorbers_on_saved_catalogue():
    rows = {}
    for f in ['extragalactic', 'ism-hot-core']:
        parse_splatalogue(json.loads((ROOT/f'targets/ic860/alma/catalogue_spw0/catalogue/{f}_raw.json').read_text()), f, rows)
    rows = list(rows.values())
    ch = .011112
    feats = [dict(id='F02', spw='0', sign='Absorption', centre_GHz=166.5999305, lo_GHz=166.5999305-1.5*ch, hi_GHz=166.5999305+1.5*ch),
             dict(id='F05', spw='0', sign='Absorption', centre_GHz=167.3269134, lo_GHz=167.3269134-3*ch, hi_GHz=167.3269134+3*ch)]
    out = A.identify(feats, rows, Z, {'0': (166.0238, 167.8686)})
    names = {k: [c['species'] for c in v['candidates']] for k, v in out.items()}
    assert 'H2S' in names['F02'] and any('30SiO' in n for n in names['F05'])
    # ADMIT's rule: lowest upper-state energy among the filtered candidates.
    for k, v in out.items():
        w = v['winner']
        assert v['rule'].startswith('lowest upper-state energy')
        assert w['Eu_K'] == min(c['Eu_K'] for c in v['candidates'] if A.isotope_count(c['species']) == 0 or A.rough_mass(c['species']) <= 50)
    assert out['F05']['winner']['species'] == 'CH3CH2CN v = 0'


def test_no_candidate_is_unidentified():
    feats = [dict(id='X', spw='s', sign='Emission', centre_GHz=100., lo_GHz=99.999, hi_GHz=100.001)]
    assert A.identify(feats, [], 0., {'s': (99., 101.)}, tier1=[])['X']['status'] == 'unidentified'
