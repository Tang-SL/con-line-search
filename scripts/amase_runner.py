"""Run Astro AMASE on one spectrum and write a JSON summary for con-line-search.

This script runs inside the Python environment where Astro AMASE is installed
(https://github.com/zfried/astro_amase); it must not import con-line-search.

    python amase_runner.py config.json

config.json: {"spectrum_path", "directory_path", "params": {...}, "output_json"}
'params' are passed unchanged to astro_amase.assign_observations().
"""
import json
import sys
import traceback


def candidate_record(c):
    keys = ['formula', 'smiles', 'isotope', 'quantum_number', 'catalog_tag', 'linelist',
            'catalog_frequency', 'observed_frequency', 'observed_intensity', 'simulated_intensity',
            'max_simulated_intensity', 'structural_score', 'frequency_match_score', 'global_score',
            'softmax_score', 'combined_score']
    out = {}
    for k in keys:
        v = getattr(c, k, None)
        try:
            out[k] = v.item() if hasattr(v, 'item') else v
        except Exception:
            out[k] = str(v)
    out['penalties'] = list(getattr(c, 'penalties', []) or [])
    return out


def summarise(results):
    assigner = results['assigner']
    lines = []
    for line in assigner.lines:
        status = getattr(line, 'assignment_status', None)
        ranked = sorted(line.candidates, key=lambda c: getattr(c, 'combined_score', 0.) or 0., reverse=True)
        lines.append(dict(frequency_MHz=float(line.frequency), intensity=float(line.intensity),
                          status=getattr(status, 'value', str(status)),
                          assigned_molecule=getattr(line, 'assigned_molecule', None),
                          best=candidate_record(line.best_candidate) if getattr(line, 'best_candidate', None) else None,
                          top=[candidate_record(c) for c in getattr(line, 'top_candidates', []) or []],
                          candidates=[candidate_record(c) for c in ranked[:10]]))
    kept = []
    cd = results.get('column_density_results')
    if cd is not None and hasattr(cd, 'to_dict'):
        kept = [dict(molecule=str(r.get('molecule')), column_density=float(r.get('column_density')))
                for r in cd.to_dict('records')]
    return dict(ok=True, vlsr=results.get('vlsr'), temperature=results.get('temperature'),
                linewidth=results.get('linewidth'), rms=None, lines=lines, kept_molecules=kept,
                output_files=results.get('output_files', {}),
                execution_time=results.get('execution_time'))


def main(path):
    config = json.load(open(path))
    try:
        import astro_amase
        params = dict(config.get('params', {}))
        if isinstance(params.get('rms_noise'), dict):
            params['rms_noise'] = {tuple(float(x) for x in k.split('|')): v for k, v in params['rms_noise'].items()}
        results = astro_amase.assign_observations(spectrum_path=config['spectrum_path'],
                                                  directory_path=config['directory_path'], **params)
        summary = summarise(results)
    except Exception as error:
        summary = dict(ok=False, error=f'{type(error).__name__}: {error}', traceback=traceback.format_exc())
    with open(config['output_json'], 'w') as stream:
        json.dump(summary, stream, indent=1, default=str)
    return 0 if summary['ok'] else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
