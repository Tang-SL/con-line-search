"""Minimal stand-in for Astro AMASE used by the tests (same result structure)."""
from enum import Enum
from types import SimpleNamespace
import numpy as np


class Status(Enum):
    ASSIGNED = 'assigned'
    UNIDENTIFIED = 'unidentified'


def _cand(formula, freq, g, s):
    return SimpleNamespace(formula=formula, smiles='C', isotope=0, quantum_number='J=2-1', catalog_tag=1, linelist='CDMS',
                           catalog_frequency=freq, observed_frequency=freq, observed_intensity=1., simulated_intensity=1.,
                           max_simulated_intensity=1., structural_score=99., frequency_match_score=.99,
                           global_score=g, softmax_score=s, combined_score=s, penalties=['example penalty'])


def assign_observations(spectrum_path, directory_path, **params):
    data = np.loadtxt(spectrum_path)
    f, t = data[:, 0], data[:, 1]
    peak = f[np.argmax(t)]
    best = _cand('TestMol', peak, 98.0, .95)
    lines = [SimpleNamespace(frequency=float(peak), intensity=float(t.max()), assignment_status=Status.ASSIGNED,
                             assigned_molecule='TestMol', best_candidate=best, top_candidates=[best],
                             candidates=[best, _cand('Other', peak+1., 50., .05)])]
    import pandas as pd
    return dict(assigner=SimpleNamespace(lines=lines), vlsr=params.get('vlsr', 0.), temperature=params['temperature'],
                linewidth=10., column_density_results=pd.DataFrame(dict(molecule=['TestMol'], column_density=[1e14])),
                output_files={}, execution_time=1.)
