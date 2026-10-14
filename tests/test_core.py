import io
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from conlines.core import read_cube,extract,detect,moment,fits_bytes,C


def synthetic(descending=True,velocity=False):
    rng=np.random.default_rng(32);a=rng.normal(0,.01,(65,81,81)).astype('float32')
    yy,xx=np.indices((81,81));source=np.exp(-((xx-40)**2+(yy-40)**2)/8)
    a[29:34]-=.8*source
    h=fits.Header();h['CTYPE1']='RA---TAN';h['CTYPE2']='DEC--TAN';h['CTYPE3']='FREQ';h['CRPIX1']=41;h['CRPIX2']=41;h['CRPIX3']=1;h['CRVAL1']=198.;h['CRVAL2']=24.;h['CRVAL3']=170e9;h['CDELT1']=-.01/3600;h['CDELT2']=.01/3600;h['CDELT3']=-1e7 if descending else 1e7;h['CUNIT1']='deg';h['CUNIT2']='deg';h['CUNIT3']='Hz';h['BUNIT']='Jy/beam';h['BMAJ']=.02/3600;h['BMIN']=.02/3600;h['RESTFRQ']=170e9;h['SPECSYS']='LSRK'
    if velocity:h['CTYPE3']='VRAD';h['CUNIT3']='m/s';h['CRVAL3']=0;h['CDELT3']=C*1000*1e7/170e9
    stream=io.BytesIO();fits.PrimaryHDU(a,h).writeto(stream);stream.seek(0);return read_cube(stream)


def test_frequency_velocity_equivalence():
    assert np.allclose(synthetic().frequency,synthetic(velocity=True).frequency)


def test_injected_absorption_and_signed_map():
    c=synthetic();e=extract(c,40,40,'Circle',.025)
    prior=[dict(species='test',rest_GHz=c.frequency[31],observed_GHz=c.frequency[31])]
    result=detect(c,e,prior,halfwidth=70,threshold=8)
    assert any(q['significant'] and q['sign']=='Absorption' and q['lo']<=31<=q['hi'] for q in result)
    image,h=moment(c,29,33,c.frequency[31],0,'linear')
    assert image[40,40]<0
    with fits.open(io.BytesIO(fits_bytes(image,h)),checksum=True) as hd:
        assert hd[0].data.ndim==2 and hd[0].verify_checksum()==1
        assert WCS(hd[0].header).pixel_n_dim==2
        assert 'km/s' in hd[0].header['BUNIT']


def test_width_sign_independent_of_axis_direction():
    c=synthetic();image,_=moment(c,29,33,170,0,'none')
    c.frequency=c.frequency[::-1].copy();c.data=c.data[::-1].copy()
    reversed_image,_=moment(c,31,35,170,0,'none')
    assert np.allclose(image,reversed_image)


def test_invalid_aperture_and_map_limits():
    c=synthetic()
    with pytest.raises(ValueError):extract(c,40,40,'Circle',1.)
    with pytest.raises(ValueError):moment(c,-1,8,170,0)
    c.data[:,40,40]=np.nan
    with pytest.raises(ValueError):extract(c,40,40,'Single pixel',.01)


def test_fractional_pixel_aperture():
    c=synthetic();e=extract(c,40.3,39.7,'Ellipse',.025,.7,35)
    assert len(e['controls'])>=8


def test_linear_continuum_removed():
    c=synthetic();c.data=np.broadcast_to(np.arange(65)[:,None,None]*.1+2,c.data.shape).copy()
    image,_=moment(c,29,33,170,0,'linear')
    assert np.max(abs(image))<1e-10

def test_both_signs_in_same_catalogue_window():
    c=synthetic();yy,xx=np.indices((81,81));source=np.exp(-((xx-40)**2+(yy-40)**2)/8)
    c.data[39:44]+=source
    ex=extract(c,40,40,'Circle',.025)
    prior=[dict(species='test',rest_GHz=c.frequency[36],observed_GHz=c.frequency[36])]
    found=detect(c,ex,prior,halfwidth=220,threshold=8)
    assert {q['sign'] for q in found}=={'Emission','Absorption'}
