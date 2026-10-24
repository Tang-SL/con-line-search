import sys
from pathlib import Path
import numpy as np
import pytest
from conlines import amase_bridge as B
from conlines.core import C

STUB = Path(__file__).resolve().parent/'stub_amase'


def test_units_velocity_and_spectrum():
    assert B.jy_per_beam_to_kelvin([1.], 'Jy/beam', [100.], 1., 1.)[0] == pytest.approx(122.2)
    assert B.jy_per_beam_to_kelvin([1000.], 'mJy/beam', [100.], 1., 1.)[0] == pytest.approx(122.2)
    with pytest.raises(ValueError):
        B.jy_per_beam_to_kelvin([1.], 'Jy/pixel', [100.], 1., 1.)
    z = .01294
    v = B.radio_velocity(z)
    assert 100*(1-v/C) == pytest.approx(100/(1+z))
    parts = [dict(frequency_GHz=np.array([101., 100.]), spectrum=np.array([0., -1.]), noise=np.array([.1, .1]),
                  unit='K', bmaj_arcsec=1., bmin_arcsec=1.)]
    f, t, rms = B.build_spectrum(parts, sign=-1.)
    assert list(f) == [100000., 101000.] and list(t) == [1., 0.] and rms == {(100000., 101000.): .1}


def test_runner_end_to_end_with_stub(tmp_path, monkeypatch):
    monkeypatch.setenv('PYTHONPATH', str(STUB))
    data = tmp_path/'db'; data.mkdir()
    f = np.linspace(100000., 100100., 101); t = np.exp(-((f-100050.)/3)**2)
    summary = B.run(sys.executable, data, f, t, dict(temperature=150., vlsr=0., rms_noise={(100000., 100100.): .01}), tmp_path/'run')
    assert summary['ok'] and summary['lines'][0]['status'] == 'assigned'
    assert summary['lines'][0]['best']['formula'] == 'TestMol' and summary['kept_molecules'][0]['molecule'] == 'TestMol'
    feats = [dict(id='F1', spw='s', sign='Emission', lo_GHz=100.049, hi_GHz=100.051),
             dict(id='F2', spw='s', sign='Emission', lo_GHz=100.08, hi_GHz=100.09),
             dict(id='F3', spw='s', sign='Absorption', lo_GHz=100.049, hi_GHz=100.051)]
    mapped = B.map_to_features(summary, feats, 'Emission', {'s': .001})
    assert mapped['F1']['suggestion']['species'] == 'TestMol' and mapped['F1']['candidates'][0]['amase_global'] == 98.
    assert mapped['F2']['status'] == 'no AMASE line' and 'F3' not in mapped
    # A molecule removed by AMASE's best-fit model is not suggested.
    summary['kept_molecules'] = [dict(molecule='Something else', column_density=1.)]
    assert B.map_to_features(summary, feats, 'Emission', {'s': .001})['F1']['suggestion'] is None


def test_missing_environment_is_reported(tmp_path):
    with pytest.raises(ValueError, match='Python executable'):
        B.run(tmp_path/'nope', tmp_path, [1.], [1.], {})
