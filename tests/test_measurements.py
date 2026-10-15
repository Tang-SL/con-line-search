from types import SimpleNamespace
import numpy as np
import pytest
from conlines.core import C, assess
from conlines.measurements import measure_feature, matching_frequency
from conlines.solutions import feature_fingerprint
from conlines.espada_backend import match_features


@pytest.mark.parametrize('sign', ['Emission','Absorption'])
@pytest.mark.parametrize('descending', [True,False])
def test_aperture_centroid_width_errors_and_stable_identity(sign,descending):
    rng=np.random.default_rng(41)
    f=np.linspace(99.95,100.05,101);centre=100.001;sigma=.003
    controls=rng.normal(0,.015,(32,101))
    # Correlated noise must be propagated as spectra, not independent channels.
    controls[:,1:]+=.5*controls[:,:-1]
    controls=np.vstack([controls,-controls])
    spectrum=(1 if sign=='Emission' else -1)*np.exp(-.5*((f-centre)/sigma)**2)
    spectrum+=.2+.3*(f-100)
    if descending:f=f[::-1];spectrum=spectrum[::-1];controls=controls[:,::-1]
    cube=SimpleNamespace(frequency=f)
    ex=dict(spectrum=spectrum,controls=controls)
    snr,models,variation=assess(ex,35,65)
    q=dict(sign=sign,lo=35,hi=65,frequency_GHz=100.02,snr=snr,flags=[],variation=variation)
    original=feature_fingerprint('cube',q)
    q['aperture_measurement']=measure_feature(cube,ex,q)
    m=q['aperture_measurement']
    assert m['reliable'] and abs(m['frequency_GHz']-centre)<.0001
    assert m['sigma_kms']==pytest.approx(C*sigma/centre,rel=.02)
    assert 0<m['centroid_error_kms']<2
    assert 0<m['sigma_error_kms']<2
    assert matching_frequency(q)==m['frequency_GHz']
    assert feature_fingerprint('cube',q)==original
    # The nuclear measurement, not the deliberately displaced 3D centroid, wins.
    q['id']='F1'
    matches=match_features([q],[dict(species='aperture',rest_GHz=centre),dict(species='3D',rest_GHz=100.02)],0,10)
    assert [r['species'] for r in matches['F1']]==['aperture']


def test_unmeasurable_feature_retains_explicit_detector_fallback():
    cube=SimpleNamespace(frequency=np.linspace(100,101,30))
    ex=dict(spectrum=np.zeros(30),controls=np.zeros((10,30)))
    q=dict(sign='Emission',lo=10,hi=14,frequency_GHz=100.4,snr=1,flags=[])
    q['aperture_measurement']=measure_feature(cube,ex,q)
    assert not q['aperture_measurement']['reliable']
    assert q['aperture_measurement']['flags'] and matching_frequency(q)==100.4
