import copy,json
import pytest
from conlines.solutions import suggest_assignments,make_solution,restore_solution,cube_fingerprint
from test_core import synthetic


def test_nearest_suggestions_and_manual_unidentified_survive_refresh():
    a=dict(species='A',rest_GHz=100.,qns='1',velocity_offset_kms=30.)
    b=dict(species='B',rest_GHz=100.1,qns='2',velocity_offset_kms=-5.)
    suggestions=suggest_assignments(['f1','f2'],{'f1':[a,b]})
    assert suggestions['f1']['transition']==b and not suggestions['f1']['reviewed']
    suggestions['f1']=dict(transition=None,reviewed=True,origin='manual')
    assert suggest_assignments(['f1'],{'f1':[a,b]},suggestions)['f1']['transition'] is None
    assert suggestions['f2']['transition'] is None
    suggestions['f1']=dict(transition=b,reviewed=True,origin='manual')
    assert suggest_assignments(['f1'],{'f1':[]},suggestions)['f1']['transition']==b


def test_solution_roundtrip_remaps_display_ids_and_checks_context():
    ctx={'redshift':.01,'data':'abc','aperture':.12}
    assignments=suggest_assignments(['old-id'],{})
    assignments['old-id']['reviewed']=True
    payload=make_solution('test',ctx,{'old-id':'stable-key'},assignments,{'old-id':(3,7)})
    saved=json.loads(json.dumps(payload))
    a,i=restore_solution(saved,ctx,{'new-id':'stable-key'},{'new-id':30})
    assert a['new-id']==assignments['old-id'] and i['new-id']==(3,7)
    for key,value in [('redshift',.02),('data','other'),('aperture',.15)]:
        changed=dict(ctx);changed[key]=value
        with pytest.raises(ValueError,match='does not match'):restore_solution(saved,changed,{'new-id':'stable-key'},{'new-id':30})
    with pytest.raises(ValueError,match='differ'):restore_solution(saved,ctx,{'new-id':'different-feature'},{'new-id':30})
    saved['features'][0]['hi']=100
    with pytest.raises(ValueError,match='interval'):restore_solution(saved,ctx,{'new-id':'stable-key'},{'new-id':30})


def test_cube_fingerprint_changes_with_pixels_or_wcs():
    cube=synthetic();copy_cube=copy.deepcopy(cube)
    assert cube_fingerprint(cube)==cube_fingerprint(copy_cube)
    copy_cube.data[0,0,0]+=1
    assert cube_fingerprint(cube)!=cube_fingerprint(copy_cube)
    copy_cube=copy.deepcopy(cube);copy_cube.header['CRVAL1']+=.1
    assert cube_fingerprint(cube)!=cube_fingerprint(copy_cube)
