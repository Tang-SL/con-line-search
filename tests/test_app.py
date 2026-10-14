from pathlib import Path
import pytest
from streamlit.testing.v1 import AppTest

ROOT=Path(__file__).resolve().parents[1]

def test_empty_app():
    app=AppTest.from_file(str(ROOT/'app.py')).run()
    assert not app.exception
    assert 'Add a cube' in app.info[0].value

@pytest.mark.skipif(not (ROOT/'targets/ic860/alma/espada_b5_all/inputs/spw0_beam060.fits').exists(),reason='Local ALMA data are not bundled')
def test_example_detection_adoption_export_and_invalidation(monkeypatch):
    import conlines.core
    monkeypatch.setattr(conlines.core,"query_catalogue",lambda *a: (_ for _ in ()).throw(RuntimeError("offline test")))
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=60).run()
    app.radio[0].set_value('IC 860 prepared example').run()
    app.button[0].click().run()
    assert not app.exception and not app.error, [e.value for e in app.error]
    analysis=app.session_state['analysis'];assert analysis['results']['example0']
    first=analysis['results']['example0'][0]['id']
    next(m for m in app.multiselect if m.label=='Adopt features for mapping').set_value([first]).run()
    assert not app.exception
    next(b for b in app.button if b.label=='Generate maps & download bundle').click().run()
    assert not app.exception and not app.error, [e.value for e in app.error]
    assert any('Generated 1 signed maps' in s.value for s in app.success)
    next(n for n in app.number_input if n.label=='Radius / half-width (arcsec)').set_value(.04).run()
    assert any('Inputs or detection settings changed' in s.value for s in app.info)
    assert not any(b.label=='Generate maps & download bundle' for b in app.button)

@pytest.mark.skipif(not (ROOT/'targets/ic860/alma/espada_b5_all/inputs/spw0_beam060.fits').exists(),reason='Local ALMA data are not bundled')
def test_catalogue_changes_and_failure_preserve_detections(monkeypatch):
    import conlines.core
    monkeypatch.setattr(conlines.core,'query_catalogue',lambda *a: (_ for _ in ()).throw(RuntimeError('offline test')))
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=120).run()
    app.radio[0].set_value('IC 860 prepared example').run()
    app.button[0].click().run()
    original=app.session_state['analysis']
    next(m for m in app.multiselect if m.label=='Astronomical filters').set_value(['Hot Cores']).run()
    assert app.session_state['analysis']==original
    next(c for c in app.checkbox if c.label=='Also match provisional features').check().run()
    next(b for b in app.button if b.label=='Match features with live Splatalogue').click().run()
    assert not app.exception
    assert app.session_state['analysis']==original
    assert 'offline test' in app.session_state['identifications']['status'][0]['Status']
    assert app.session_state['identifications']['matches']=={}
