import numpy as np
from types import SimpleNamespace
from conlines import blends
from conlines.core import C

Z=0.04
def channels(rest,n=200,dv=5.):
    e=rest/(1+Z);df=e*dv/C
    return e+df*(np.arange(n)-n/2)


def spectrum(freq,lines,sigma=60.,v0=10.,noise=.002,seed=1):
    y=np.zeros(len(freq))
    for rest,a in lines:
        e=rest/(1+Z);v=C*(e-freq)/e;y+=a*np.exp(-.5*((v-v0)/sigma)**2)
    return y+np.random.default_rng(seed).normal(0,noise,len(freq))


def comp(key,rest,region=0,sign='Emission',tie=None):
    return dict(key=key,row=dict(species=key,rest_GHz=rest,qns=''),expected_GHz=rest/(1+Z),region=region,sign=sign,tie=tie)


def test_two_blended_lines_are_separated_and_a_spurious_one_is_rejected():
    hcn,vib,other=265.886434,265.852709,265.70
    f=channels(hcn,300);y=spectrum(f,[(hcn,.05),(vib,.02)])
    regions=[dict(id='F1',spw='a',freq=f,spectrum=y,noise=np.full(len(f),.002),lo=0,hi=len(f)-1)]
    fit=blends.decompose(regions,[comp('HCN',hcn),comp('vib',vib),comp('X',other)],'gauss',Z)
    q={c['key']:c for c in fit['components']}
    # 40 km/s apart under a 140 km/s FWHM: the two lines trade flux through the shared
    # velocity. The fit must say so instead of reporting a confident split.
    assert set(q['HCN']['degenerate_with'])=={'vib'} and q['HCN']['verdict']=='Needed, but not separable'
    assert q['vib']['verdict']=='Needed, but not separable' and q['X']['verdict']=='Not needed'
    assert q['HCN']['amplitude_error']>.01
    assert abs(q['HCN']['amplitude']+q['vib']['amplitude']-.07)<.006        # the total is well determined
    a=blends.assignments_from_fit(fit)['F1']
    assert {a['transition']['species'],*[r['species'] for r in a['blend']]}=={'HCN','vib'}


def test_tie_across_two_windows_constrains_the_blend():
    hcn,e_,f_,hco=265.886434,265.852709,267.199283,267.557626
    fa=channels(hcn,300);fb=channels(hco,300)
    ya=spectrum(fa,[(hcn,.05),(e_,.02)],seed=2);yb=spectrum(fb,[(f_,.02),(hco,.06)],seed=3)
    regions=[dict(id='A',spw='a',freq=fa,spectrum=ya,noise=np.full(300,.002),lo=0,hi=299),
             dict(id='B',spw='b',freq=fb,spectrum=yb,noise=np.full(300,.002),lo=0,hi=299)]
    comps=[comp('HCN',hcn,0),comp('e',e_,0,tie=('f',1.)),comp('f',f_,1),comp('HCO+',hco,1)]
    fit=blends.decompose(regions,comps,'gauss',Z)
    q={c['key']:c for c in fit['components']}
    assert abs(q['f']['amplitude']-.02)<.004 and q['e']['amplitude']==q['f']['amplitude']
    assert abs(q['HCN']['amplitude']-.05)<.006 and q['HCN']['amplitude_error']<.004   # the tie breaks the degeneracy
    assert q['f']['verdict']=='Needed'


def test_self_absorbed_profile_with_two_gaussians():
    rest=267.557626;f=channels(rest,300)
    e=rest/(1+Z);v=C*(e-f)/e
    y=.06*np.exp(-.5*(v/120)**2)-.04*np.exp(-.5*((v-20)/30)**2)+np.random.default_rng(4).normal(0,.002,300)
    regions=[dict(id='F',spw='a',freq=f,spectrum=y,noise=np.full(300,.002),lo=0,hi=299)]
    one=blends.decompose(regions,[comp('HCO+',rest)],'gauss',Z)
    two=blends.decompose(regions,[comp('HCO+',rest)],'gauss2',Z)
    assert two['reduced_chi2']<1.5<one['reduced_chi2'] and two['params'][4]<0


def test_absorption_amplitudes_are_negative_and_inventory_uses_range():
    rest=265.886434;f=channels(rest,100);y=-spectrum(f,[(rest,.03)],noise=.001)
    regions=[dict(id='F',spw='a',freq=f,spectrum=y,noise=np.full(100,.001),lo=0,hi=99)]
    fit=blends.decompose(regions,[comp('HCN',rest,sign='Absorption')],'gauss',Z)
    assert fit['components'][0]['amplitude']<0 and fit['components'][0]['flux']<0
    feature=dict(lo=40,hi=60,frequency_GHz=float(f[50]),sign='Absorption')
    rows=[dict(species='in',rest_GHz=rest,qns=''),dict(species='edge',rest_GHz=f[62]*(1+Z),qns=''),dict(species='out',rest_GHz=f[95]*(1+Z),qns='')]
    names=[r['species'] for r in blends.inventory(feature,f,rows,Z,margin_kms=30)]
    assert names==sorted(names,key=lambda s:['edge','in'].index(s)) or set(names)=={'in','edge'}
    assert 'out' not in names


def test_component_maps_recover_per_pixel_amplitudes_and_ties():
    from astropy.io import fits
    hcn,vib=265.886434,265.70
    f=channels(hcn,120)
    yy,xx=np.indices((9,9));amp1=.05*np.exp(-((xx-4)**2+(yy-4)**2)/6);amp2=.02*(xx/8)
    e1,e2=hcn/(1+Z),vib/(1+Z)
    prof=lambda e:np.exp(-.5*(((C*(e-f)/e)-0)/40)**2)
    data=(prof(e1)[:,None,None]*amp1+prof(e2)[:,None,None]*amp2+np.random.default_rng(5).normal(0,1e-4,(120,9,9))).astype('float32')
    h=fits.Header();h['CTYPE1']='RA---TAN';h['CTYPE2']='DEC--TAN';h['CRPIX1']=5;h['CRPIX2']=5;h['CRVAL1']=10;h['CRVAL2']=5;h['CDELT1']=-1e-5;h['CDELT2']=1e-5
    cube=SimpleNamespace(data=data,frequency=f,header=h)
    regions=[dict(id='F',spw='a',freq=f,spectrum=data[:,4,4],noise=np.full(120,1e-3),lo=0,hi=119)]
    fit=blends.decompose(regions,[comp('A',hcn),comp('B',vib)],'gauss',Z)
    integral=fit['integral_kms']
    maps,diag=blends.component_maps(fit,{0:cube},{0:np.full(120,1e-3)},local=False)
    assert np.allclose(maps['A'][0]/integral,amp1,atol=1e-3) and np.allclose(maps['B'][0]/integral,amp2,atol=1e-3)
    maps,diag=blends.component_maps(fit,{0:cube},{0:np.full(120,1e-3)})
    sig=np.isfinite(diag[0]['velocity']);assert sig.sum()>=20
    assert np.allclose(maps['A'][0][sig]/integral,amp1[sig],rtol=.06,atol=2e-4) and np.allclose(maps['B'][0][sig]/integral,amp2[sig],rtol=.06,atol=2e-4)


def test_local_maps_follow_a_velocity_gradient_without_leaking_into_the_blend():
    from astropy.io import fits
    a_,b_=265.886434,265.70            # B is 210 km/s from A
    f=channels(a_,160)
    yy,xx=np.indices((11,11));amp=.05*np.exp(-((xx-5)**2+(yy-5)**2)/8);vfield=(xx-5)*18.   # rotation +-90 km/s, within half the 210 km/s separation
    ea=a_/(1+Z);va=C*(ea-f)/ea
    data=(np.exp(-.5*((va[:,None,None]-vfield)/40)**2)*amp+np.random.default_rng(6).normal(0,1e-4,(160,11,11))).astype('float32')    # only line A is present
    h=fits.Header();h['CTYPE1']='RA---TAN';h['CTYPE2']='DEC--TAN';h['CRPIX1']=6;h['CRPIX2']=6;h['CRVAL1']=10;h['CRVAL2']=5;h['CDELT1']=-1e-5;h['CDELT2']=1e-5
    cube=SimpleNamespace(data=data,frequency=f,header=h)
    aperture=data[:,4:7,4:7].mean(axis=(1,2))
    regions=[dict(id='F',spw='a',freq=f,spectrum=aperture,noise=np.full(160,1e-3),lo=0,hi=159)]
    fit=blends.decompose(regions,[comp('A',a_),comp('B',b_)],'gauss',Z)
    fixed,_=blends.component_maps(fit,{0:cube},{0:np.full(160,1e-3)},local=False)
    local,diag=blends.component_maps(fit,{0:cube},{0:np.full(160,1e-3)})
    leak_fixed=np.nanmax(np.abs(fixed['B'][0]));leak_local=np.nanmax(np.abs(local['B'][0]))
    assert leak_local<.25*leak_fixed                                  # B stays (nearly) empty when the shape follows the gradient
    assert abs(diag[0]['velocity'][5,8]-54)<8 and abs(diag[0]['velocity'][5,2]+54)<8
    assert abs(np.nansum(local['A'][0])/np.nansum(amp*np.sqrt(2*np.pi)*40)-1)<.1   # A's flux is kept
    assert not diag[0]['limited']
    # A velocity field beyond half the separation is reported, not silently mis-split.
    steep=(np.exp(-.5*((va[:,None,None]-(xx-5)*45.)/40)**2)*amp+np.random.default_rng(7).normal(0,1e-4,(160,11,11))).astype('float32')
    _,diag2=blends.component_maps(fit,{0:SimpleNamespace(data=steep,frequency=f,header=h)},{0:np.full(160,1e-3)})
    assert diag2[0]['limited']
