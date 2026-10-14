from types import SimpleNamespace
import numpy as np
from conlines.plots import chart_image,spectrum_plot


def test_map_zoom_centres_on_aperture_and_preserves_full_image():
    im=np.zeros((101,201));fig=chart_image(im,'test',center=(70,60),zoom=3.)
    assert np.isclose(np.diff(fig.layout.xaxis.range)[0],201/3)
    assert np.isclose(np.mean(fig.layout.xaxis.range),70)
    assert np.isclose(np.mean(fig.layout.yaxis.range),60)
    assert np.asarray(fig.data[0].z).shape==im.shape
    edge=chart_image(im,'edge',center=(1,1),zoom=3.)
    assert edge.layout.xaxis.range[0]>=-.5


def test_all_candidates_labelled_and_significance_distinguished():
    cube=SimpleNamespace(frequency=np.linspace(166,168,20),frame='LSRK',unit='Jy/beam')
    ex=dict(spectrum=np.zeros(20),noise=np.ones(20))
    rows=[dict(id=f'spw-{i:03}',lo=i+2,hi=i+4,frequency_GHz=cube.frequency[i+3],snr=10 if i else 4,significant=bool(i),matches=[{'species':name}]) for i,name in enumerate(['H2S','CH3OH'])]
    fig=spectrum_plot(cube,ex,rows)
    assert len(fig.layout.annotations)==2
    assert 'spw-000' in fig.layout.annotations[0].text and 'spw-001' in fig.layout.annotations[1].text
    assert 'H2S' not in fig.layout.annotations[0].text
    assert fig.layout.annotations[0].font.color!=fig.layout.annotations[1].font.color
    markers=[q for q in fig.data if q.mode=='markers']
    assert {q.marker.symbol for q in markers}=={'diamond','circle-open'}
