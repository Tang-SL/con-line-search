from types import SimpleNamespace
import numpy as np
from conlines.espada_backend import match_features,aperture_features
from conlines.core import C


def test_matches_use_centroids_and_keep_detection_unchanged():
    features=[dict(id='F1',frequency_GHz=100.,significant=True,snr=12.)]
    rows=[dict(species='near',rest_GHz=110.01),dict(species='far',rest_GHz=111.),dict(species='exact',rest_GHz=110.)]
    result=match_features(features,rows,.1,100.)
    assert [r['species'] for r in result['F1']]==['exact','near']
    assert np.isclose(result['F1'][1]['velocity_offset_kms'],C*(110.01/1.1-100)/(110.01/1.1))
    assert features==[dict(id='F1',frequency_GHz=100.,significant=True,snr=12.)]


def test_mask_intersection_and_absorption_sign_screen():
    rng=np.random.default_rng(3)
    controls=rng.normal(0,1,(64,60));spectrum=np.zeros(60);spectrum[28:33]=-20
    ex=dict(mask=np.array([[False,True],[False,False]]),controls=controls,spectrum=spectrum,noise=np.ones(60))
    cube=SimpleNamespace(frequency=np.linspace(100,101,60))
    components=[dict(sign=sign,lo=28,hi=32,frequency_GHz=100.5,footprint=fp) for sign,fp in [('Absorption',[1]),('Emission',[1]),('Absorption',[3])]]
    found=aperture_features(cube,ex,dict(components=components))
    assert len(found)==2
    assert found[0]['significant']
    assert not found[1]['significant'] and 'aperture sign mismatch' in found[1]['flags']
