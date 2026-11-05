"""Step: find line features in every full cube with ESPADA/SoFiA (no aperture, no catalogue)."""
import numpy as np
from ..core import read_cube
from ..espada_backend import run_espada
from . import figs


def run(run, widths=(3, 5, 9), clip=5., mode='Both', contsub=True, noise='local', only=None):
    out = run.load('components', {})
    for key, record in run.cubes().items():
        if only and key not in only:
            continue
        cube = read_cube(record['path'])
        result = run_espada(cube, widths, mode, clip, contsub, noise)
        comps = sorted(result['components'], key=lambda c: (c['frequency_GHz'], c['sign']))
        for n, c in enumerate(comps):
            c['id'] = f'{key}-F{n+1:03}'
        out[key] = dict(components=comps, status=result['status'], settings=result['settings'],
                        searched_box=result.get('searched_box'))
        image = np.nanmax(np.abs(cube.data-np.nanmedian(cube.data, axis=0)), axis=0)
        figs.components(run.plot(f'detect_{key}'), image, comps, key)
        (run.path/'cache'/f'espada_{key}.log').write_text(result.get('console_log', ''))
        run.log('detect', f'{key}: {len(comps)} components', status=result['status'], settings=result['settings'])
        del cube
    run.save('components', out)
    return out


def summary(out):
    lines = []
    for key, r in out.items():
        signs = {}
        for c in r['components']:
            signs[c['sign']] = signs.get(c['sign'], 0)+1
        lines.append(f"{key}: {len(r['components'])} components ({', '.join(f'{v} {k.lower()}' for k, v in signs.items()) or 'none'}); {r['status']}")
    return '\n'.join(lines)
