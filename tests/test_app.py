from pathlib import Path
import json
import numpy as np
import pytest
from astropy.io import fits
from streamlit.testing.v1 import AppTest
from test_core import synthetic

ROOT=Path(__file__).resolve().parents[1]


def by_label(elements,label):return next(e for e in elements if e.label==label)


@pytest.fixture
def workbench(tmp_path,monkeypatch):
    import conlines.espada_backend
    c=synthetic();yy,xx=np.indices((81,81));c.data[39:44]+=.8*np.exp(-((xx-40)**2+(yy-40)**2)/8)
    path=tmp_path/'cube.fits';fits.writeto(path,c.data,c.header)
    components=[dict(sofia_id=1,sign=sign,lo=lo,hi=hi,x=40.,y=40.,frequency_GHz=float(c.frequency[(lo+hi)//2]),footprint=[40*81+40]) for sign,lo,hi in [('Absorption',29,33),('Emission',39,43)]]
    monkeypatch.setattr(conlines.espada_backend,'run_espada',lambda *a,**kw:dict(components=components,status={},parameters={},console_log='UI fixture'))
    monkeypatch.setenv('CONLINES_SOLUTION_DIR',str(tmp_path/'solutions'))
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=60).run()
    app.radio[0].set_value('Local file paths').run()
    app.text_area[0].set_value(str(path)).run()
    by_label(app.button,'Run detection across all loaded SPWs').click().run()
    assert not app.exception and not app.error,[e.value for e in app.error]
    return app,tmp_path


def test_empty_app():
    app=AppTest.from_file(str(ROOT/'app.py')).run()
    assert not app.exception
    assert 'Add a cube' in app.info[0].value


def test_detection_without_catalogue_and_unidentified_export(workbench):
    app,_=workbench
    analysis=app.session_state['analysis'];first=next(iter(analysis['results'].values()))[0]['id']
    by_label(app.multiselect,'Adopt features for mapping').set_value([first]).run()
    by_label(app.button,'Generate maps & download bundle').click().run()
    assert not app.exception and not app.error,[e.value for e in app.error]
    assert any('Generated 1 signed maps' in s.value for s in app.success)
    original_signature=app.session_state['analysis']['signature']
    by_label(app.slider,'Aperture preview zoom').set_value(5.).run()
    assert app.session_state['analysis']['signature']==original_signature
    by_label(app.number_input,'Radius / half-width (arcsec)').set_value(.04).run()
    assert any('Inputs or detection settings changed' in s.value for s in app.info)
    assert not any(b.label=='Generate maps & download bundle' for b in app.button)


def test_catalogue_failure_preserves_detections(workbench,monkeypatch):
    import conlines.core
    app,_=workbench
    monkeypatch.setattr(conlines.core,'query_catalogue',lambda *a: (_ for _ in ()).throw(RuntimeError('offline test')))
    original=app.session_state['analysis']
    by_label(app.multiselect,'Astronomical filters').set_value(['Hot Cores']).run()
    assert app.session_state['analysis']==original
    by_label(app.button,'Match features with live Splatalogue').click().run()
    assert not app.exception and app.session_state['analysis']==original
    assert 'offline test' in app.session_state['identifications']['status'][0]['Status']
    assert app.session_state['identifications']['matches']=={}


def test_batch_suggestions_sequential_review_save_load_and_maps(workbench,monkeypatch):
    import conlines.core
    app,tmp=workbench;features=next(iter(app.session_state['analysis']['results'].values()))
    z=by_label(app.number_input,'Redshift z').value
    rows=[dict(species=f'Molecule {i}',rest_GHz=f['frequency_GHz']*(1+z),qns='J=1-0',Eu_K=10.,log10_A_s=-4.,filters=['extragalactic']) for i,f in enumerate(features)]
    monkeypatch.setattr(conlines.core,'query_catalogue',lambda *a:rows)
    inspect=by_label(app.selectbox,'Inspect feature').value
    by_label(app.button,'Match features with live Splatalogue').click().run()
    assert not app.exception
    by_label(app.button,'Review identifications one by one').click().run()
    assert by_label(app.selectbox,'Identification choice').value==0
    assert any('Evidence for this identification' in m.value for m in app.markdown)
    assert any('Other compatible features' in m.value for m in app.markdown)
    by_label(app.selectbox,'Identification choice').set_value(-1).run()
    by_label(app.button,'Confirm & next').click().run()
    assert by_label(app.selectbox,'Inspect feature').value==inspect
    assert by_label(app.selectbox,'Identification choice').value==0
    by_label(app.button,'Confirm & finish').click().run()
    assert any('Review complete' in m.value for m in app.success)
    by_label(app.button,'Save solution locally').click().run()
    saved=list((tmp/'solutions').glob('*.json'));assert len(saved)==1
    payload=json.loads(saved[0].read_text());assert all(q['assignment']['reviewed'] for q in payload['features'])
    assert payload['measurement_method']=='signed aperture moments v1'
    assert payload['features'][0]['aperture_measurement']['reliable']
    assert payload['evidence_snapshot']['queried_utc']
    assert payload['features'][0]['assignment']['transition'] is None
    assert payload['features'][1]['assignment']['transition']['species']=='Molecule 1'
    payload['features'][0]['lo']=28;payload['features'][0]['hi']=34
    saved[0].write_text(json.dumps(payload))
    # Mutate a reviewed assignment, then restore the saved version.
    by_label(app.button,'Review identifications one by one').click().run()
    by_label(app.selectbox,'Identification choice').set_value(0).run()
    by_label(app.button,'Confirm & next').click().run()
    by_label(app.button,'Load solution').click().run()
    assert not app.exception
    by_label(app.button,'Review identifications one by one').click().run()
    assert by_label(app.selectbox,'Identification choice').value==-1
    headers=[];original_fits=conlines.core.fits_bytes
    def capture(image,header):
        headers.append(header.copy());return original_fits(image,header)
    monkeypatch.setattr(conlines.core,'fits_bytes',capture)
    by_label(app.multiselect,'Adopt features for mapping').set_value([f['id'] for f in features]).run()
    by_label(app.button,'Generate maps & download bundle').click().run()
    assert not app.exception and any('Generated 2 signed maps' in m.value for m in app.success)

    assert headers[0]['CHANLO']==28 and headers[0]['CHANHI']==34
    assert 'RESTFRQ' not in headers[0] and 'VREFGHZ' in headers[0]
    assert all(h['IDSTAT']=='REVIEWED' for h in headers)
    assert headers[1]['LINEID']=='Molecule 1'


def test_full_coverage_query_includes_spw_without_features(workbench,monkeypatch):
    import conlines.core
    import conlines.espada_backend
    app,tmp=workbench
    original=conlines.espada_backend.run_espada
    monkeypatch.setattr(conlines.espada_backend,'run_espada',lambda cube,*a,**kw:dict(components=[],status={},parameters={},console_log='empty SPW') if cube.frequency.max()>171 else original(cube,*a,**kw))
    c=synthetic();c.header['CRVAL3']=172e9
    second=tmp/'second.fits';fits.writeto(second,c.data,c.header)
    app.text_area[0].set_value(str(tmp/'cube.fits')+'\n'+str(second)).run()
    by_label(app.button,'Run detection across all loaded SPWs').click().run()
    assert not app.exception
    z=by_label(app.number_input,'Redshift z').value
    calls=[]
    def query(low,high,filters):
        calls.append((low/(1+z),high/(1+z)))
        return [dict(species='M',species_id=1,rest_GHz=f*(1+z),qns=str(i),Eu_K=10.,log10_A_s=-4.,filters=['extragalactic']) for i,f in enumerate([169.69,171.7])]
    monkeypatch.setattr(conlines.core,'query_catalogue',query)
    by_label(app.button,'Match features with live Splatalogue').click().run()
    assert not app.exception
    assert len(calls)==2
    assert calls[0][0]<169.36 and calls[0][1]>170
    assert calls[1][0]<171.36 and calls[1][1]>172
    options=next(v for v in app.session_state['identifications']['matches'].values() if v)
    assert any(e['SPW']=='second.fits' and e['Status']=='No ESPADA feature' for e in options[0]['evidence']['transitions'])
