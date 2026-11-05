"""What kind of source is this? SIMBAD, NED (redshift, type, photometry) and NASA ADS literature.

Network services are only reached from here; every result is cached in the run folder so the
report can show exactly what was used. All functions return plain dicts.
"""
from pathlib import Path
import math
import os
import re
import requests
from ..core import C

ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------- SIMBAD / NED
def simbad(name=None, ra=None, dec=None, radius_arcsec=10.):
    from astroquery.simbad import Simbad
    s = Simbad()
    s.add_votable_fields('otype', 'rvz_redshift', 'ids', 'morphtype')
    if name:
        table = s.query_object(name)
    else:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        table = s.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius_arcsec*u.arcsec)
    if table is None or len(table) == 0:
        return None
    r = table[0]
    def get(*names):
        for n in names:
            if n in table.colnames:
                v = r[n]
                return None if getattr(v, 'mask', False) else (v.item() if hasattr(v, 'item') else v)
        return None
    return dict(main_id=str(get('main_id', 'MAIN_ID')), otype=str(get('otype', 'OTYPE')),
                ra=float(get('ra', 'RA_d')), dec=float(get('dec', 'DEC_d')),
                redshift=get('rvz_redshift', 'RVZ_REDSHIFT'), morph=str(get('morph_type', 'morphtype', 'MORPH_TYPE') or ''),
                ids=str(get('ids', 'IDS') or '').split('|'))


def ned(name):
    from astroquery.ipac.ned import Ned
    table = Ned.query_object(name)
    if table is None or len(table) == 0:
        return None
    r = table[0]
    out = dict(name=str(r['Object Name']), type=str(r['Type']), redshift=float(r['Redshift']) if r['Redshift'] is not None and str(r['Redshift']) not in ('--', 'nan') else None,
               ra=float(r['RA']), dec=float(r['DEC']))
    try:
        phot = Ned.get_table(name, table='photometry')
        rows = []
        for p in phot:
            freq = float(p['Frequency']) if str(p['Frequency']) not in ('--', 'nan') else None
            flux = float(p['Flux Density']) if str(p['Flux Density']) not in ('--', 'nan') else None
            if freq and flux:
                rows.append(dict(band=str(p['Observed Passband']), frequency_Hz=freq, wavelength_um=2.99792458e14/freq,
                                 flux_Jy=flux, reference=str(p['Refcode']) if 'Refcode' in phot.colnames else None))
        out['photometry'] = rows
    except Exception as e:     # photometry is optional
        out['photometry'] = []
        out['photometry_error'] = str(e)
    return out


def iras_fluxes(photometry):
    """IRAS 12/25/60/100 µm flux densities (Jy) from NED photometry, if present."""
    out = {}
    for p in photometry or []:
        band = p['band'].lower()
        if 'iras' in band:
            for w in (12, 25, 60, 100):
                if re.search(rf'\b{w}\b', band) and w not in out:
                    out[w] = p['flux_Jy']
    return out


def derived(z, photometry, beam_arcsec=None):
    """Distance, scale and IR luminosity (Sanders & Mirabel 1996) where possible."""
    from astropy.cosmology import Planck18
    import astropy.units as u
    out = {}
    if z is None or z <= 0:
        return out
    dl = Planck18.luminosity_distance(z).to_value(u.Mpc)
    scale = Planck18.kpc_proper_per_arcmin(z).to_value(u.pc/u.arcsec)
    out.update(redshift=z, cz_kms=C*z, luminosity_distance_Mpc=dl, pc_per_arcsec=scale)
    if beam_arcsec:
        out['beam_pc'] = beam_arcsec*scale
    f = iras_fluxes(photometry)
    if all(w in f for w in (12, 25, 60, 100)):
        fir = 1.8e-14*(13.48*f[12]+5.16*f[25]+2.58*f[60]+f[100])     # W m^-2
        lir = 4*math.pi*(dl*3.0857e22)**2*fir/3.828e26
        out.update(iras_Jy=f, L_IR_Lsun=lir, log_L_IR=math.log10(lir), f25_f60=f[25]/f[60] if f[60] else None,
                   class_=('ULIRG' if lir >= 1e12 else 'LIRG' if lir >= 1e11 else 'sub-LIRG'))
    elif 60 in f and 25 in f:
        out.update(iras_Jy=f, f25_f60=f[25]/f[60])
    return out


# ---------------------------------------------------------------- NASA ADS
def ads_token():
    for env in ('ADS_API_TOKEN', 'ADS_DEV_KEY'):
        if os.environ.get(env):
            return os.environ[env].strip()
    for p in (Path.home()/'.ads'/'dev_key', ROOT/'.local-tools'/'ads_token'):
        if p.exists():
            return p.read_text().strip()
    return None


def ads_search(query, rows=25, sort='citation_count desc', token=None):
    token = token or ads_token()
    if not token:
        raise ValueError('No ADS token: set ADS_API_TOKEN, or put the token in ~/.ads/dev_key or .local-tools/ads_token.')
    r = requests.get('https://api.adsabs.harvard.edu/v1/search/query', timeout=60,
                     headers=dict(Authorization=f'Bearer {token}'),
                     params=dict(q=query, rows=rows, sort=sort, fl='bibcode,title,author,year,abstract,keyword,identifier,citation_count,doi,pub'))
    r.raise_for_status()
    docs = r.json()['response']['docs']
    out = []
    for d in docs:
        arxiv = next((i.split(':', 1)[1] for i in d.get('identifier', []) if i.lower().startswith('arxiv:')), None)
        out.append(dict(bibcode=d['bibcode'], title=(d.get('title') or [''])[0], year=d.get('year'),
                        first_author=(d.get('author') or [''])[0], n_authors=len(d.get('author') or []),
                        journal=d.get('pub'), citations=d.get('citation_count'), keywords=d.get('keyword', []),
                        abstract=d.get('abstract', ''), arxiv=arxiv,
                        url=f"https://ui.adsabs.harvard.edu/abs/{d['bibcode']}/abstract",
                        arxiv_url=f'https://arxiv.org/abs/{arxiv}' if arxiv else None))
    return out


def standard_queries(name, aliases=()):
    """Literature searches the agent runs for every source (it can add its own)."""
    names = [name, *[a for a in aliases if a and a != name]][:4]
    obj = ' OR '.join(f'object:"{n}"' for n in names)
    title = ' OR '.join(f'abs:"{n}"' for n in names)
    return {
        'molecular lines in this source': f'({obj}) AND abs:(molecular OR molecule OR "line survey" OR HCN OR HCO+ OR CO OR absorption)',
        'ALMA / interferometer studies': f'({title}) AND abs:(ALMA OR NOEMA OR SMA OR PdBI OR interferometer)',
        'nucleus: obscuration, AGN, starburst': f'({obj}) AND abs:(nucleus OR AGN OR obscured OR "compact obscured" OR starburst OR outflow)',
        'most recent papers': f'({obj})',
    }


def literature(names, extra=(), rows=15):
    queries = standard_queries(names[0], names[1:])
    for i, q in enumerate(extra):
        queries[f'extra {i+1}'] = q
    results, seen = {}, set()
    for label, q in queries.items():
        try:
            docs = ads_search(q, rows=rows, sort='date desc' if label == 'most recent papers' else 'citation_count desc')
        except Exception as e:
            results[label] = dict(query=q, error=str(e), papers=[])
            continue
        results[label] = dict(query=q, papers=docs)
        seen.update(d['bibcode'] for d in docs)
    return dict(queries=results, n_unique=len(seen))
