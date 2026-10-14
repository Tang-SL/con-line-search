from dataclasses import dataclass
from pathlib import Path
import io, json, html, re
import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from radio_beam import Beam
from scipy.signal import fftconvolve

C = 299792.458

@dataclass
class Cube:
    data: np.ndarray
    header: fits.Header
    frequency: np.ndarray
    pixel_arcsec: tuple
    beam_arcsec: float
    unit: str
    frame: str
    notes: list


def read_cube(source):
    """Accept conventional RA/DEC/spectral FITS with optional singleton Stokes."""
    with fits.open(source, memmap=False) as hd:
        h = hd[0].header.copy()
        if not h.get('CTYPE1','').startswith('RA') or not h.get('CTYPE2','').startswith('DEC'):
            raise ValueError('Expected FITS axes RA, DEC, spectral (in that order). Reorder other cubes first.')
        a = np.asarray(hd[0].data)
        while a.ndim > 3 and a.shape[0] == 1: a = a[0]
        if a.ndim != 3 or min(a.shape) < 2:
            raise ValueError('Use a 3D cube, optionally with a singleton Stokes axis; multiple polarizations are not combined.')
        w = WCS(h).sub([1,2,3]); spectral = w.sub(['spectral'])
        values = np.asarray(spectral.pixel_to_world_values(np.arange(len(a))),float)
        kind = h.get('CTYPE3','').split('-')[0]
        unit = spectral.wcs.cunit[0]
        rest = h.get('RESTFRQ',h.get('RESTFREQ',0))
        if kind == 'FREQ': f = (values*u.Unit(unit)).to_value(u.GHz)
        elif kind in ['VRAD','VOPT']:
            if not rest: raise ValueError('Velocity cubes need RESTFRQ to recover observed frequencies.')
            velocity = (values*u.Unit(unit)).to_value(u.km/u.s)
            f = rest/1e9 * (1-velocity/C) if kind=='VRAD' else rest/1e9/(1+velocity/C)
        else: raise ValueError(f'Unsupported spectral axis {kind}; use FREQ, VRAD or VOPT.')
        if not np.all(np.isfinite(f)) or not (np.all(np.diff(f)>0) or np.all(np.diff(f)<0)):
            raise ValueError('Frequency axis must be finite and strictly monotonic.')
        unit_flux = h.get('BUNIT','').strip()
        if unit_flux.lower().replace(' ','') not in ['jy/beam','mjy/beam','k','jy/pixel']:
            raise ValueError(f'Unsupported brightness unit {unit_flux!r}; use Jy/beam, mJy/beam, K or Jy/pixel.')
        a = a.astype('float32',copy=True)
        if unit_flux.lower().replace(' ','')=='mjy/beam': a /= 1000; unit_flux='Jy/beam'
        scales = tuple(abs(float(q.to_value(u.arcsec))) for q in w.celestial.proj_plane_pixel_scales())
        notes=[]
        if 'BEAMS' in hd:
            table=hd['BEAMS']; rows=table.data
            if len(rows)!=len(a): raise ValueError('BEAMS table does not have exactly one beam per channel.')
            def beam_values(key):
                unit_b=table.columns[key].unit or ('deg' if key=='BPA' else 'arcsec')
                return (np.array(rows[key])*u.Unit(unit_b)).to_value(u.deg)
            majors,minors,pas=beam_values('BMAJ'),beam_values('BMIN'),beam_values('BPA')
        elif h.get('BMAJ') and h.get('BMIN'):
            majors=np.full(len(a),h['BMAJ']);minors=np.full(len(a),h['BMIN']);pas=np.full(len(a),h.get('BPA',0))
        else:
            if '/beam' in unit_flux.lower(): raise ValueError('Jy/beam cubes require BMAJ/BMIN or a BEAMS table.')
            majors=minors=np.full(len(a),max(scales)/3600);pas=np.zeros(len(a))
            notes.append('No beam supplied; spatial controls use pixel-scale spacing.')
        variable=not (np.allclose(majors,majors[0],rtol=.001) and np.allclose(minors,minors[0],rtol=.001) and np.allclose(pas,pas[0],atol=.01))
        if variable:
            if not np.isclose(scales[0],scales[1],rtol=.001): raise ValueError('Variable-beam convolution currently requires square pixels.')
            if not np.all(np.isfinite(a)): raise ValueError('Variable-beam cubes with blank pixels must be homogenized externally first.')
            target=Beam(major=float(max(majors)*1.001)*u.deg,minor=float(max(majors)*1.001)*u.deg,pa=0*u.deg)
            for i in range(len(a)):
                beam=Beam(major=majors[i]*u.deg,minor=minors[i]*u.deg,pa=pas[i]*u.deg)
                kernel=target.deconvolve(beam).as_kernel(scales[0]*u.arcsec).array
                area=(target.sr/beam.sr).value if '/beam' in unit_flux.lower() else 1
                a[i]=fftconvolve(a[i],kernel/kernel.sum(),mode='same')*area
            h['BMAJ']=target.major.to_value(u.deg);h['BMIN']=target.minor.to_value(u.deg);h['BPA']=0
            notes.append('Channel beams convolved to a common circular beam; edge apertures should be avoided.')
        else:
            h['BMAJ']=float(majors[0]);h['BMIN']=float(minors[0]);h['BPA']=float(pas[0])
        h['BUNIT']=unit_flux
        return Cube(a,h,f,scales,h['BMAJ']*3600,unit_flux,h.get('SPECSYS','UNKNOWN'),notes)


def aperture_mask(shape, x, y, kind, radius, ratio=1., angle=0., scales=(1.,1.)):
    yy,xx=np.indices(shape);dx=(xx-x)*scales[0];dy=(yy-y)*scales[1]
    t=np.deg2rad(angle);xp=dx*np.cos(t)+dy*np.sin(t);yp=-dx*np.sin(t)+dy*np.cos(t)
    if kind=='Single pixel':
        mask=np.zeros(shape,bool)
        if 0<=round(y)<shape[0] and 0<=round(x)<shape[1]: mask[round(y),round(x)]=True
        return mask
    if kind=='Rectangle': return (abs(xp)<=radius)&(abs(yp)<=radius*ratio)
    return (xp/radius)**2+(yp/(radius*ratio))**2<=1


def extract(cube, x, y, kind, radius, ratio=1., angle=0.):
    a=cube.data;ny,nx=a.shape[1:]
    extent=radius*(2**.5 if kind=='Rectangle' else 1) if kind!='Single pixel' else 0
    if min(x*cube.pixel_arcsec[0],(nx-1-x)*cube.pixel_arcsec[0],y*cube.pixel_arcsec[1],(ny-1-y)*cube.pixel_arcsec[1])<extent:
        raise ValueError('The aperture extends outside the image. Move it inward or reduce its size.')
    mask=aperture_mask((ny,nx),x,y,kind,radius,ratio,angle,cube.pixel_arcsec)
    if not mask.any(): raise ValueError('The aperture contains no pixels.')
    def spectrum(m):
        # A fixed spatial mask is used at all frequencies: no changing NaN weights.
        values=a[:,m]
        return np.mean(values,axis=1)
    spec=spectrum(mask)
    if not np.all(np.isfinite(spec)):raise ValueError('Aperture intersects blank pixels; choose a smaller or different aperture.')
    controls=[];centers=[]
    size=max(radius,cube.beam_arcsec/2,max(cube.pixel_arcsec))
    spacing=2.2*size
    step_x=max(1,int(np.ceil(spacing/cube.pixel_arcsec[0])));step_y=max(1,int(np.ceil(spacing/cube.pixel_arcsec[1])))
    for cy in np.arange(y % step_y + step_y,ny-step_y,step_y):
        for cx in np.arange(x % step_x + step_x,nx-step_x,step_x):
            distance=np.hypot((cx-x)*cube.pixel_arcsec[0],(cy-y)*cube.pixel_arcsec[1])
            if distance<3*size or distance>8*size or len(controls)>=64:continue
            m=aperture_mask((ny,nx),cx,cy,kind,radius,ratio,angle,cube.pixel_arcsec)
            if m.sum()!=mask.sum():continue
            ss=spectrum(m)
            if np.all(np.isfinite(ss)):controls.append(ss);centers.append([cx,cy])
    if len(controls)<8:raise ValueError('Fewer than 8 non-overlapping control apertures fit outside the source. Reduce aperture size or upload a wider field.')
    controls=np.asarray(controls)
    noise=np.maximum(np.std(controls,axis=0,ddof=1),1.4826*np.median(abs(controls-np.median(controls,axis=0)),axis=0))
    return dict(spectrum=spec,controls=controls,noise=noise,mask=mask,control_centers=centers,pixels=int(mask.sum()))


def query_catalogue(low, high, filters):
    rows={}
    for filt in filters:
        payload=dict(searchSpeciesAutocomplete='',quickPickerMolecules=[],telescopeBands=['any'],energyFrom='',energyTo='',energyValue='el-1',frequenciesFrom=[str(low)],frequenciesTo=[str(high)],frequenciesUnit='GHz',frequencyRedshift=0,astronomicalFilters=filt,exportType='current',exportDelimiter='tab',exportLimit='allRecords',exportStart=1,exportStop=250)
        response=requests.post('https://splatalogue.online/splata-slap/basicsearch/false',json={'body':json.dumps(payload)},timeout=90)
        response.raise_for_status();data=response.json()
        if not isinstance(data,list):raise ValueError('Unexpected Splatalogue response.')
        if len(data)==1 and data[0].get('searchErrorMessage') and not data[0].get('species_id'):raise ValueError(data[0]['searchErrorMessage'])
        for r in data:
            if not r.get('species_id'):continue
            clean=lambda s:html.unescape(re.sub('<[^>]*>','',s or '')).strip()
            freq=float(r['orderedFreq'].split(',')[0]);key=(r['species_id'],freq,clean(r['resolved_QNs']))
            if key in rows:rows[key]['filters'].append(filt);continue
            rows[key]=dict(species=clean(r['species_name']),rest_GHz=freq,qns=clean(r['resolved_QNs']),Eu_K=float(r['euk']) if r.get('euk') else None,log10_A_s=float(r['aij']) if r.get('aij') else None,filters=[filt])
    return list(rows.values())


def select_priors(rows, low, high, z, halfwidth, energy, filters, exclude=True):
    out=[]
    for r in rows:
        if not set(r.get('filters',filters))&set(filters):continue
        if r.get('Eu_K') is None or r['Eu_K']>energy:continue
        if exclude and re.search(r'D(?![a-z])|\btrans\b|^t-',r['species']):continue
        fc=r['rest_GHz']/(1+z)
        if fc*(1+halfwidth/C)<low or fc*(1-halfwidth/C)>high:continue
        out.append(dict(r,observed_GHz=fc))
    return out


def sidebands(n,lo,hi,width=10):
    left=np.arange(max(0,lo-2-width),max(0,lo-2));right=np.arange(min(n,hi+3),min(n,hi+3+width))
    return np.r_[left,right],len(left)==width and len(right)==width


def residual(values,lo,hi,model):
    sides,complete=sidebands(values.shape[-1],lo,hi,5 if model=='constant' else 10)
    if sum(sides<lo)<3 or sum(sides>hi)<3:raise ValueError('Insufficient baseline channels on both sides.')
    chan=np.arange(lo,hi+1);mid=(lo+hi)/2
    if model=='constant':base=np.median(values[...,sides],axis=-1)[...,None]
    else:
        x=sides-mid;xm=x.mean();v=values[...,sides];av=v.mean(axis=-1)
        slope=np.sum((x-xm)*v,axis=-1)/np.sum((x-xm)**2)
        base=av[...,None]+slope[...,None]*(chan-mid-xm)
    return values[...,chan]-base,sides,complete


def assess(extraction,lo,hi):
    vals=np.vstack([extraction['spectrum'],extraction['controls']]);models=[]
    for model in ['constant','linear']:
        r,sides,complete=residual(vals,lo,hi,model);summed=r.sum(axis=-1);ctrl=summed[1:]
        noise=max(np.std(ctrl,ddof=1),1.4826*np.median(abs(ctrl-np.median(ctrl))))
        signal=summed[0]-np.median(ctrl)
        chnoise=np.maximum(np.std(r[1:],axis=0,ddof=1),1.4826*np.median(abs(r[1:]-np.median(r[1:],axis=0)),axis=0))
        if noise<=0 or np.any(chnoise<=0):raise ValueError('Cannot estimate nonzero noise.')
        ch=(r[0]-np.median(r[1:],axis=0))/chnoise
        support=bool(np.any((np.sign(signal)*ch[:-1]>=3)&(np.sign(signal)*ch[1:]>=3)))
        models.append(dict(snr=float(signal/noise),signal=float(signal),complete=complete,support=support,sidebands=sides.tolist()))
    same=models[0]['signal']*models[1]['signal']>0
    score=min(abs(m['snr']) for m in models) if same else 0.
    change=abs(models[0]['signal']-models[1]['signal'])/max(abs(m['signal']) for m in models)
    return score,models,change


def detect(cube,extraction,priors,halfwidth=200,threshold=8,widths=(3,5,9),mode='Both'):
    """Targeted 1D boxcar search; local empirical S/N, not ESPADA/SoFiA."""
    if mode=='Both':
        return sorted(detect(cube,extraction,priors,halfwidth,threshold,widths,'Emission')+detect(cube,extraction,priors,halfwidth,threshold,widths,'Absorption'),key=lambda q:q['snr'],reverse=True)
    features=[];evaluations={};f=cube.frequency;n=len(f)
    for p in priors:
        centers=np.where(abs(C*(p['observed_GHz']-f)/p['observed_GHz'])<=halfwidth)[0];best=None
        for center in centers:
            for width in widths:
                lo=int(center-width//2);hi=lo+width-1
                if lo<0 or hi>=n:continue
                key=(lo,hi)
                if key not in evaluations:
                    try:evaluations[key]=assess(extraction,lo,hi)
                    except ValueError:evaluations[key]=None
                result=evaluations[key]
                if result is None:continue
                score,models,variation=result;sign='Emission' if models[0]['signal']>0 else 'Absorption'
                if mode!='Both' and sign!=mode:continue
                if best is None or score>best['snr']:
                    flags=[]
                    if not all(m['complete'] for m in models):flags.append('baseline coverage')
                    if not all(m['support'] for m in models):flags.append('channel support')
                    if variation>.3:flags.append('baseline sensitive')
                    if np.any(extraction['noise'][lo:hi+1]>3*np.median(extraction['noise'])):flags.append('high channel noise')
                    best=dict(species=p['species'],rest_GHz=p['rest_GHz'],expected_GHz=p['observed_GHz'],frequency_GHz=float(f[center]),lo=lo,hi=hi,snr=score,sign=sign,flags=flags,models=models,variation=variation,significant=score>=threshold and not flags,matches=[p])
        if best is None:continue
        # Report provisional >=3 features too, while keeping significance explicit.
        if best['snr']<3:continue
        other=next((q for q in features if q['sign']==best['sign'] and max(q['lo'],best['lo'])<=min(q['hi'],best['hi'])),None)
        if other:
            matches=other['matches']+[p]
            if best['snr']>other['snr']:other.update(best)
            other['matches']=matches
        else:features.append(best)
    return sorted(features,key=lambda q:q['snr'],reverse=True)


def moment(cube,lo,hi,rest,z,baseline='linear'):
    if not 0<=lo<=hi<len(cube.data):raise ValueError('Integration channels are outside this cube.')
    if rest<=0 or z<=-1:raise ValueError('Invalid rest frequency or redshift.')
    values=np.moveaxis(cube.data,0,-1)
    if baseline=='none':r=values[...,lo:hi+1];sides=[]
    else:r,sides,_=residual(values,lo,hi,baseline)
    freq=cube.frequency;edges=np.r_[freq[0]-(freq[1]-freq[0])/2,(freq[1:]+freq[:-1])/2,freq[-1]+(freq[-1]-freq[-2])/2]
    dv=C*np.abs(np.diff(edges))/(rest/(1+z))
    image=np.sum(r*dv[lo:hi+1],axis=-1)
    h=WCS(cube.header).celestial.to_header()
    for k in ['BMAJ','BMIN','BPA','OBJECT','SPECSYS']:
        if k in cube.header:h[k]=cube.header[k]
    h['BUNIT']=cube.unit+' km/s';h['RESTFRQ']=float(rest*1e9);h['REDSHIFT']=float(z);h['CHANLO']=lo;h['CHANHI']=hi;h['BASEMOD']=baseline
    h.add_history('Signed moment 0; no clipping or S/N mask. Channel indices are zero-based.')
    h.add_history('Velocity widths: c * abs(dnu) / (rest frequency / (1+z)).')
    h.add_history('Baseline channels: '+','.join(map(str,sides)))
    return image,h


def fits_bytes(image,header):
    out=io.BytesIO();fits.PrimaryHDU(image.astype('float32'),header).writeto(out,checksum=True);return out.getvalue()
