from types import SimpleNamespace
import numpy as np
from conlines.evidence import add_evidence


def feature(ident,f,lo=49,hi=51,sign='Emission',significant=True):
    return dict(id=ident,frequency_GHz=f,lo=lo,hi=hi,sign=sign,significant=significant,
                aperture_measurement=dict(reliable=True,frequency_GHz=f))


def setup_evidence():
    loaded={k:dict(name=k,cube=SimpleNamespace(frequency=np.linspace(f-1,f+1,101),frame='LSRK')) for k,f in [('s0',100.),('s1',110.)]}
    row=lambda f,q:dict(species='Molecule v=0',species_id=42,rest_GHz=f,qns=q)
    a=row(100.,'a');b=row(110.,'b');blend=row(110.001,'c');missing=row(110.5,'d')
    lookup={'a':('s0',feature('a',100.)), 'b':('s1',feature('b',110.))}
    matches={'a':[dict(a,velocity_offset_kms=0.)]}
    catalogues={'s0':[a], 's1':[b,blend,missing]}
    return matches,lookup,catalogues,loaded


def test_companions_count_features_not_blended_rows_and_missing_not_penalized():
    matches,lookup,catalogues,loaded=setup_evidence()
    result=add_evidence(matches,lookup,catalogues,loaded,0)
    ev=result['a'][0]['evidence']
    assert ev['other_significant_features']==['b']
    assert ev['other_provisional_features']==[]
    assert any(r['Status']=='No ESPADA feature' for r in ev['transitions'])
    assert 'evidence' not in matches['a'][0]
    assert 'score' not in result['a'][0]
    lookup['b'][1]['significant']=False
    ev=add_evidence(matches,lookup,catalogues,loaded,0)['a'][0]['evidence']
    assert not ev['other_significant_features'] and ev['other_provisional_features']==['b']


def test_opposite_sign_frames_unlinked_apertures_and_missing_detection_do_not_support():
    matches,lookup,catalogues,loaded=setup_evidence()
    lookup['b'][1]['sign']='Absorption'
    ev=add_evidence(matches,lookup,catalogues,loaded,0)['a'][0]['evidence']
    assert not ev['other_significant_features']
    assert any(r['Status']=='Opposite-sign feature' for r in ev['transitions'])
    lookup['b'][1]['sign']='Emission'
    for kwargs,status in [({'linked_apertures':False},'Apertures are not linked'),({'analysed_spws':['s0']},'Detection unavailable')]:
        ev=add_evidence(matches,lookup,catalogues,loaded,0,**kwargs)['a'][0]['evidence']
        assert not ev['other_significant_features'] and any(r['Status']==status for r in ev['transitions'])
    loaded['s1']['cube'].frame='BARYCENT'
    ev=add_evidence(matches,lookup,catalogues,loaded,0)['a'][0]['evidence']
    assert not ev['other_significant_features']
    assert any(r['Status']=='Frame comparison unavailable' for r in ev['transitions'])
    del catalogues['s1']
    ev=add_evidence(matches,lookup,catalogues,loaded,0)['a'][0]['evidence']
    assert ev['unavailable_spws']==['s1']


def test_companions_follow_candidate_velocity_and_do_not_merge_isotopologues():
    matches,lookup,catalogues,loaded=setup_evidence()
    for key,(spw,q) in lookup.items():
        q['aperture_measurement']['frequency_GHz']*=1-120/299792.458
    # Widen the detector intervals to contain the shifted centroids.
    lookup['a'][1].update(lo=45,hi=51);lookup['b'][1].update(lo=45,hi=51)
    ev=add_evidence(matches,lookup,catalogues,loaded,0,velocity_tolerance=10)['a'][0]['evidence']
    assert ev['other_significant_features']==['b']
    for r in catalogues['s1']:r['species_id']=43;r['species']='13Molecule v=0'
    ev=add_evidence(matches,lookup,catalogues,loaded,0)['a'][0]['evidence']
    assert not ev['other_significant_features']
