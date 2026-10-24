"""Run Astro AMASE (Fried et al., arXiv:2608.18221) on aperture spectra.

Astro AMASE pins older numpy/scipy/astropy than this app, so it runs in its own
Python environment as a separate process (scripts/amase_runner.py). This module
prepares its input, runs it, and maps its per-line assignments onto ESPADA
features.

Absorption: AMASE models emission. For absorption features the spectrum is
sign-flipped and run separately. For optically thin lines from one excitation
temperature, absorption depths and emission intensities within one band scale
with the same optical depths, so relative line strengths, and hence AMASE's
consistency checks, carry over. Column densities and temperatures from the
flipped run are not physical.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import tempfile
import numpy as np
from .core import C

RUNNER = Path(__file__).resolve().parents[1]/'scripts'/'amase_runner.py'


def radio_velocity(z):
    """Radio-convention velocity (km/s) giving the same frequency shift as 1+z.

    AMASE/molsim shift frequencies as f_obs = f_rest (1 - v/c); with
    f_obs = f_rest/(1+z) that is v = c z/(1+z) exactly.
    """
    return C*z/(1+z)


def jy_per_beam_to_kelvin(values, unit, frequency_GHz, bmaj_arcsec, bmin_arcsec):
    unit = (unit or '').strip()
    values = np.asarray(values, float)
    if unit.lower() == 'k':
        return values
    if unit.lower() in ('jy/beam', 'mjy/beam'):
        jy = values/1000. if unit.lower() == 'mjy/beam' else values
        return 1.222e6*jy/(np.asarray(frequency_GHz, float)**2*bmaj_arcsec*bmin_arcsec)
    raise ValueError(f'Astro AMASE needs K or Jy/beam spectra; this cube is in {unit}.')


def build_spectrum(parts, sign=1.):
    """Concatenate aperture spectra (one per SPW) into AMASE's input.

    parts: dicts with frequency_GHz, spectrum, noise, unit, bmaj_arcsec, bmin_arcsec.
    Returns (frequency_MHz increasing, intensity_K, rms dict {(fmin, fmax) MHz: K}).
    """
    freqs, ints, rms = [], [], {}
    for p in parts:
        f = np.asarray(p['frequency_GHz'], float)
        t = sign*jy_per_beam_to_kelvin(p['spectrum'], p['unit'], f, p['bmaj_arcsec'], p['bmin_arcsec'])
        n = jy_per_beam_to_kelvin(p['noise'], p['unit'], f, p['bmaj_arcsec'], p['bmin_arcsec'])
        good = np.isfinite(f) & np.isfinite(t)
        freqs.append(f[good]*1000.); ints.append(t[good])
        rms[(float(f.min()*1000.), float(f.max()*1000.))] = float(np.nanmedian(np.abs(n)))
    f = np.concatenate(freqs); t = np.concatenate(ints)
    order = np.argsort(f)
    return f[order], t[order], dict(sorted(rms.items()))


def run(python, data_dir, frequency_MHz, intensity_K, params, workdir=None, timeout=7200):
    """Write the spectrum, call the runner in AMASE's environment, return its summary."""
    if not python or not Path(python).exists():
        raise ValueError('Set the Astro AMASE Python executable (the python of the environment where astro_amase is installed).')
    if not data_dir or not Path(data_dir).is_dir():
        raise ValueError('Set the Astro AMASE database folder (the files from its Dropbox link).')
    workdir = Path(workdir or tempfile.mkdtemp(prefix='conlines-amase-'))
    workdir.mkdir(parents=True, exist_ok=True)
    spectrum = workdir/'spectrum.txt'
    np.savetxt(spectrum, np.c_[frequency_MHz, intensity_K], fmt='%.6f %.8e')
    output = workdir/'amase_summary.json'
    params = dict(params)
    if isinstance(params.get('rms_noise'), dict):
        # JSON keys must be strings; the runner restores the tuples.
        params['rms_noise'] = {f'{a}|{b}': v for (a, b), v in params['rms_noise'].items()}
    config = dict(spectrum_path=str(spectrum), directory_path=str(data_dir), params=params, output_json=str(output))
    (workdir/'config.json').write_text(json.dumps(config, indent=1))
    proc = subprocess.run([str(python), str(RUNNER), str(workdir/'config.json')], cwd=str(workdir),
                          capture_output=True, text=True, timeout=timeout)
    (workdir/'amase_console.log').write_text(proc.stdout+'\n'+proc.stderr)
    if not output.exists():
        raise ValueError('Astro AMASE produced no summary: '+(proc.stderr or proc.stdout)[-1500:])
    summary = json.loads(output.read_text())
    if not summary.get('ok'):
        raise ValueError('Astro AMASE failed: '+summary.get('error', 'unknown error'))
    summary['workdir'] = str(workdir)
    return summary


def candidate_row(c, vlsr):
    rest = c['catalog_frequency']/1000./(1-vlsr/C)
    label = (f"AMASE global {c.get('global_score') or 0:.1f} · local {c.get('combined_score') or 0:.2f}")
    return dict(species=c['formula'], rest_GHz=float(rest), qns=str(c.get('quantum_number') or ''),
                catalogue=c.get('linelist') or 'AMASE', Eu_K=None, log10_A_s=None, smiles=c.get('smiles'),
                amase_global=c.get('global_score'), amase_local=c.get('combined_score'),
                amase_structural=c.get('structural_score'), amase_frequency=c.get('frequency_match_score'),
                amase_penalties=c.get('penalties', []), method_label=label)


def map_to_features(summary, features, sign, channel_GHz):
    """Attach AMASE's line results to features of one sign.

    features: dicts with id, sign, lo_GHz, hi_GHz (observed). An AMASE line belongs
    to a feature if its frequency lies within the feature's interval (+/- half a
    channel); the strongest such line is used.
    """
    vlsr = float(summary['vlsr'])
    kept = {m['molecule'] for m in summary.get('kept_molecules', [])}
    out = {}
    for f in features:
        if f['sign'] != sign:
            continue
        lo, hi = sorted([f['lo_GHz'], f['hi_GHz']])
        half = abs(channel_GHz.get(f['spw'], 0.))/2
        lines = [l for l in summary['lines'] if lo-half <= l['frequency_MHz']/1000. <= hi+half]
        if not lines:
            out[f['id']] = dict(status='no AMASE line', candidates=[], suggestion=None,
                                note='AMASE found no ≥ threshold peak inside this feature.')
            continue
        line = max(lines, key=lambda l: l['intensity'])
        candidates = [candidate_row(c, vlsr) for c in line['candidates']]
        best = candidate_row(line['best'], vlsr) if line.get('best') else None
        status = line['status']
        suggestion = None
        note = {'assigned': 'AMASE assigned this line.', 'multiple_carriers': 'AMASE: several possible carriers.',
                'unidentified': 'AMASE: unidentified.'}.get(status, status)
        if status == 'assigned' and best is not None:
            if kept and best['species'] not in kept:
                note = 'AMASE assigned this line, but the molecule was removed by its best-fit model.'
            else:
                suggestion = best
        if best is not None and not any(c['species'] == best['species'] and abs(c['rest_GHz']-best['rest_GHz']) < 1e-6 for c in candidates):
            candidates.insert(0, best)
        out[f['id']] = dict(status=status, candidates=candidates, suggestion=suggestion, note=note,
                            line_MHz=line['frequency_MHz'], line_intensity_K=line['intensity'])
    return out


def signature(*items):
    return hashlib.sha256(json.dumps(items, sort_keys=True, default=str).encode()).hexdigest()


SETTINGS = Path(__file__).resolve().parents[1]/'.local-tools'/'amase_paths.json'


def _saved():
    try:
        return json.loads(SETTINGS.read_text())
    except Exception:
        return {}


def default_python():
    return os.environ.get('AMASE_PYTHON') or _saved().get('python', '')


def default_data_dir():
    return os.environ.get('AMASE_DATA') or _saved().get('data_dir', '')


def remember_paths(python, data_dir, path=None):
    """Store the paths locally (.local-tools is git-ignored) for the next session."""
    path = Path(path or SETTINGS)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(dict(python=str(python), data_dir=str(data_dir)), indent=1))
    except OSError:
        pass
