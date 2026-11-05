"""Step: candidate transitions for every spectral window, from Splatalogue (or a local CDMS/JPL copy).

Each transition keeps what the companion-line (LTE) test needs: rest frequency, upper-state
energy, Einstein A and line strength S·µ², plus whether the species has been seen in other
galaxies (Splatalogue 'Extragalactic' list) or in space (CDMS/JPL lists).
"""
from pathlib import Path
import csv
import glob
import json
import os
import pickle
import numpy as np
from ..core import query_catalogue, C

ROOT = Path(__file__).resolve().parents[2]
VACUUM = 1.16395e-20          # A = VACUUM * nu(MHz)^3 * S mu^2(D^2) / g_u


def splatalogue(low, high, filters):
    rows = query_catalogue(low, high, tuple(filters))
    for r in rows:
        r['source'] = 'Splatalogue'
        r['seen_extragalactic'] = 'extragalactic' in r.get('filters', [])
    return rows


# ---------------------------------------------------------------- local CDMS/JPL (Astro AMASE database)
class _Any:
    def __init__(self, *a, **k):
        pass

    def __setstate__(self, s):
        self.__dict__.update(s if isinstance(s, dict) else {'state': s})


class _Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            return super().find_class(module, name)
        except Exception:
            return _Any


def amase_database():
    for p in (os.environ.get('AMASE_DATA'),):
        if p and Path(p).exists():
            return Path(p)
    saved = ROOT/'.local-tools'/'amase_paths.json'
    if saved.exists():
        d = Path(json.loads(saved.read_text()).get('data_dir', ''))
        if d.exists():
            return d
    return None


def local_cdms(low, high, folder=None):
    """Transitions between low and high GHz (rest) from the molsim pickles of the Astro AMASE database."""
    folder = Path(folder) if folder else amase_database()
    if folder is None:
        raise ValueError('No local CDMS/JPL copy found (install Astro AMASE, or set AMASE_DATA).')
    rows = []
    for cat, sub in (('CDMS', 'cdms'), ('JPL', 'jpl')):
        names = {}
        listing = folder/f'all_{sub}_final_official.csv'
        if listing.exists():
            for r in csv.DictReader(listing.open()):
                try:
                    names[int(r['tag'])] = (r['mol'], r.get('detected_in_ism', '') == 'y')
                except (KeyError, ValueError):
                    pass
        for path in glob.glob(str(folder/f'{sub}_pkl'/'*.pkl')):
            tag = int(Path(path).stem)
            try:
                obj = _Unpickler(open(path, 'rb')).load()
                c = obj.catalog.__dict__
            except Exception:
                continue
            freq = np.asarray(c['frequency'], float)/1e3
            sel = np.flatnonzero((freq >= low) & (freq <= high))
            name, in_space = names.get(tag, (getattr(obj, 'name', str(tag)), False))
            for i in sel:
                a = float(c['aij'][i]) if c.get('aij') is not None else None
                rows.append(dict(species=name, species_id=tag, catalogue=cat, rest_GHz=float(freq[i]),
                                 qns=f"{c['qnup_str'][i]}-{c['qnlow_str'][i]}", Eu_K=float(c['eup'][i]), El_K=float(c['elow'][i]),
                                 log10_A_s=float(np.log10(a)) if a and a > 0 else None,
                                 line_strength=float(c['sijmu'][i]) if c.get('sijmu') is not None else None,
                                 g_up=float(c['gup'][i]) if c.get('gup') is not None else None,
                                 log10_intensity_300K=float(c['logint'][i]), source=f'{cat} (local copy)',
                                 seen_in_space=in_space, seen_extragalactic=None))
    return rows


# ---------------------------------------------------------------- screening
def screen(rows, emax=1500., min_log_a=-6.5, recommended_only=True):
    """Keep transitions that can matter. Excited states are kept (CONs need them)."""
    out, reasons = [], {}
    for r in rows:
        why = None
        if r.get('Eu_K') is not None and r['Eu_K'] > emax:
            why = f'E_up > {emax:g} K'
        elif r.get('log10_A_s') is not None and r['log10_A_s'] < min_log_a:
            why = f'log A < {min_log_a:g}'
        elif recommended_only and r.get('recommended') in (False, 0, '0', 'false', 'False'):
            why = 'not the recommended catalogue entry'
        if why:
            reasons[why] = reasons.get(why, 0)+1
        else:
            out.append(r)
    return out, reasons


def strength(r):
    """S·µ² (D²); derived from Einstein A and g_u when only those are given."""
    if r.get('line_strength'):
        return float(r['line_strength'])
    if r.get('log10_A_s') is not None and r.get('g_up'):
        return 10**r['log10_A_s']*r['g_up']/(VACUUM*(r['rest_GHz']*1e3)**3)
    return None


def run(run, z, filters=('extragalactic', 'ism-hot-core'), margin_kms=600., emax=1500., min_log_a=-6.5, source='auto'):
    inv = run.load('inventory')['cubes']
    out, notes = {}, []
    for r in inv:
        low = r['freq_min_GHz']*(1+z)*(1-margin_kms/C); high = r['freq_max_GHz']*(1+z)*(1+margin_kms/C)
        rows, used = None, None
        if source.startswith('file:'):
            rows = [dict(x) for x in json.loads(Path(source[5:]).expanduser().read_text()) if low <= x['rest_GHz'] <= high]
            used = f'line list {Path(source[5:]).name}'
        if rows is None and source in ('auto', 'splatalogue'):
            try:
                rows = splatalogue(low, high, filters); used = 'Splatalogue'
            except Exception as e:
                if source == 'splatalogue':
                    raise
                notes.append(f"{r['key']}: Splatalogue unavailable ({e}); using the local CDMS/JPL copy.")
        if rows is None:
            rows = local_cdms(low, high); used = 'local CDMS/JPL'
        for x in rows:
            x.setdefault('seen_extragalactic', None)
        kept, reasons = screen(rows, emax, min_log_a, recommended_only=(used == 'Splatalogue'))
        for k in kept:
            k['S_mu2'] = strength(k)
        out[r['key']] = dict(source=used, rest_range_GHz=[low, high], n_raw=len(rows), n_kept=len(kept), removed=reasons, rows=kept)
    result = dict(z=z, filters=list(filters), emax_K=emax, min_log_a=min_log_a, notes=notes, windows=out)
    run.save('catalogue', result)
    run.update(z=z)
    run.log('catalog', f'{sum(w["n_kept"] for w in out.values())} transitions kept', z=z, filters=list(filters), emax=emax, source=source)
    return result


def summary(result):
    lines = [f"z = {result['z']:.6f}; filters {result['filters']}; E_up <= {result['emax_K']:g} K, log A >= {result['min_log_a']:g}"]
    for k, w in result['windows'].items():
        lines.append(f"  {k}: {w['source']}, {w['n_kept']} of {w['n_raw']} transitions kept ({w['rest_range_GHz'][0]:.3f}-{w['rest_range_GHz'][1]:.3f} GHz rest)")
    lines += ['  ! '+n for n in result['notes']]
    return '\n'.join(lines)
