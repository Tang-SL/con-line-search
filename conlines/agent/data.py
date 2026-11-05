"""Find the cubes to analyse: a local folder, or an ALMA archive project.

ALMA pipeline products follow a naming scheme (member.uid___..._sci.spwNN.cube.I.pbcor.fits,
.pb.fits(.gz), .mask.fits(.gz), and *.cont.I.pbcor.fits / *.mfs.* continuum images).
For each spectral window we keep the primary-beam-corrected cube for the line search,
plus its primary beam and mask when present.
"""
from pathlib import Path
import re

CUBE_PATTERNS = ('*.fits', '*.fits.gz', '*.fit')


def _role(name):
    n = name.lower()
    if '.cube.' in n or 'cube' in n:
        if '.pbcor.' in n:
            return 'cube'
        if re.search(r'\.pb\.fits', n):
            return 'pb'
        if '.mask.' in n:
            return 'mask'
        if '.residual.' in n or '.psf.' in n or '.model.' in n or '.sumwt.' in n:
            return 'other'
        if n.endswith(('.image.fits', '.fits')) and '.pbcor' not in n:
            return 'cube_flat'      # not primary-beam corrected
    if '.cont.' in n or '.mfs.' in n or 'cont' in n:
        if '.pbcor.' in n:
            return 'continuum'
        return 'other'
    return 'unknown'


def _spw(name):
    m = re.search(r'spw(\d+)', name)
    return f'spw{m.group(1)}' if m else None


def _target(name):
    m = re.search(r'\.([^.]+)_sci\.', name)
    return m.group(1).replace('_', ' ') if m else None


def discover(folder):
    """Group FITS files in a folder by spectral window. Returns (cubes, continuum, notes)."""
    folder = Path(folder).expanduser().resolve()
    if not folder.is_dir():
        raise ValueError(f'{folder} is not a folder.')
    files = sorted({p for pattern in CUBE_PATTERNS for p in folder.rglob(pattern)})
    groups, continuum, unknown = {}, [], []
    for p in files:
        role = _role(p.name)
        if role == 'continuum':
            continuum.append(str(p)); continue
        if role in ('other',):
            continue
        spw = _spw(p.name)
        if role == 'unknown' or spw is None:
            unknown.append(p); continue
        g = groups.setdefault((_target(p.name) or '', spw), {})
        g.setdefault(role, str(p))
    cubes, notes = [], []
    for (target, spw), g in sorted(groups.items(), key=lambda kv: (kv[0][0], int(kv[0][1][3:]))):
        path = g.get('cube') or g.get('cube_flat')
        if not path:
            continue
        if 'cube' not in g:
            notes.append(f'{spw}: only a non-primary-beam-corrected cube; fluxes away from the centre are underestimated.')
        if str(path).endswith('.gz'):
            notes.append(f'{spw}: the cube is gzip-compressed; reading it is slow and needs memory for the full cube.')
        key = spw if len({t for t, _ in groups}) <= 1 else f'{target.replace(" ", "")}-{spw}'
        cubes.append(dict(key=key, spw=spw, target=target or None, path=path, pb=g.get('pb'), mask=g.get('mask')))
    # Plain FITS cubes that do not follow ALMA naming: one per file.
    if not cubes:
        for p in unknown:
            if p.suffix.lower() in ('.fits', '.fit'):
                cubes.append(dict(key=p.stem, spw=None, target=None, path=str(p), pb=None, mask=None))
        if cubes:
            notes.append('Files do not follow ALMA pipeline naming; each FITS file is treated as one cube. Check that they are 3D spectral cubes.')
    return cubes, continuum, notes


# ---------------------------------------------------------------- ALMA archive
def alma_observations(project_code, target=None):
    """Observations of an ALMA project (needs network; uses astroquery)."""
    from astroquery.alma import Alma
    alma = Alma()
    table = alma.query(payload=dict(proposal_id=project_code), public=None)
    rows = []
    for r in table:
        name = str(r['target_name'])
        if target and target.replace(' ', '').lower() not in name.replace(' ', '').replace('_', '').lower():
            continue
        rows.append(dict(target=name, member_ous=str(r['member_ous_uid']), band=str(r['band_list']),
                         ra=float(r['s_ra']), dec=float(r['s_dec']), resolution_arcsec=float(r['s_resolution']),
                         frequency_support=str(r['frequency_support']), release=str(r['obs_release_date']),
                         public=bool(r['data_rights'] == 'Public') if 'data_rights' in r.colnames else None,
                         velocity_resolution_ms=float(r['velocity_resolution']) if 'velocity_resolution' in r.colnames else None))
    # one row per member OUS (the archive returns one per spectral window)
    unique = {}
    for r in rows:
        unique.setdefault(r['member_ous'], r)
    return list(unique.values())


def alma_products(member_ous):
    """Pipeline image products of one member OUS: name, URL, size (bytes), role, spw."""
    from astroquery.alma import Alma
    info = Alma().get_data_info(member_ous, expand_tarfiles=True)
    out = []
    for r in info:
        url = str(r['access_url'])
        name = url.rsplit('/', 1)[-1]
        if not name.lower().endswith(('.fits', '.fits.gz')):
            continue
        out.append(dict(name=name, url=url, size=int(r['content_length']) if r['content_length'] else None,
                        role=_role(name), spw=_spw(name), target=_target(name)))
    return out


def alma_download(urls, destination):
    from astroquery.alma import Alma
    destination = Path(destination).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    return Alma().download_files(list(urls), savedir=str(destination))
