from types import SimpleNamespace
from pathlib import Path
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


def test_blank_border_is_cropped_and_positions_map_back(monkeypatch):
    """pbcor cubes are blank outside the primary beam: search only the data box,
    let SoFiA see the remaining NaNs, measure noise locally, and report positions
    and footprints in the full image."""
    import subprocess, yaml
    from astropy.io import fits
    from astropy.io.votable import from_table, writeto
    from astropy.table import Table
    from conlines import espada_backend
    from test_core import synthetic
    cube=synthetic();yy,xx=np.indices(cube.data.shape[1:])
    cube.data[:,np.hypot(yy-40,xx-40)>30]=np.nan          # blank outside a circle, like a pbcor cube
    seen={}
    def fake_run(cmd,cwd,**kw):
        cfg=yaml.safe_load(open(cmd[-1]))
        data=fits.getdata(cfg['input_data_set']['cube'][0]);seen['shape']=data.shape;seen['nan']=bool(np.isnan(data).any())
        seen['par']=open(cfg['sofia_emi_file']).read()
        out=Path(cfg['output_dir'])/'espada_cube';out.mkdir(parents=True)
        mask=np.zeros(data.shape,'int32');mask[29:34,29:32,29:32]=1   # cropped coordinates: full-image centre (40,40) is (30,30)
        fits.writeto(out/'emission_cube_mask.fits',mask)
        writeto(from_table(Table(dict(id=[1],z_min=[29],z_max=[33],z=[31.],x=[30.],y=[30.]))),str(out/'emission_cube_cat.xml'))
        return subprocess.CompletedProcess(cmd,0)
    monkeypatch.setattr(espada_backend,'executables',lambda:('espada','sofia'))
    monkeypatch.setattr(espada_backend.subprocess,'run',fake_run)
    result=espada_backend.run_espada(cube,mode='Emission')
    assert seen['shape']==(65,61,61) and seen['nan']                     # cropped to the 61×61 box with data
    assert 'scaleNoise.mode = local' in seen['par'] and 'scaleNoise.windowXY = 31' in seen['par']
    c=result['components'][0]
    assert (c['x'],c['y'])==(40.,40.)
    assert set(c['footprint'])=={y*81+x for y in range(39,42) for x in range(39,42)}
    assert result['searched_box']==dict(x=[10,71],y=[10,71])


def test_all_blank_cube_is_refused(monkeypatch):
    from conlines import espada_backend
    from test_core import synthetic
    import pytest
    cube=synthetic();cube.data[:]=np.nan
    monkeypatch.setattr(espada_backend,'executables',lambda:('espada','sofia'))
    with pytest.raises(ValueError,match='no data'):espada_backend.run_espada(cube)
