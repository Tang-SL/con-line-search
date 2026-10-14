"""Assess independently rerun targeted ESPADA components and catalogue ambiguity."""
from pathlib import Path
import json,ast
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'targets/ic860/alma/catalogue_spw0';OLD=P.parent/'espada_b5_all';C=299792.458
cfg=json.loads((P/'config/settings.json').read_text());priors=json.loads((P/'catalogue/priors.json').read_text());catalog=json.loads((P/'catalogue/normalized.json').read_text())
a,h=fits.getdata(OLD/'inputs/spw0_beam060.fits',header=True)
f=h['RESTFRQ']/1e9*(1-(h['CRVAL3']+(np.arange(len(a))+1-h['CRPIX3'])*h['CDELT3'])/1000/C)
# Reuse the previously reviewed empirical-noise/baseline measurement functions only.
tree=ast.parse((ROOT/'scripts/ic860_espada_b5_rank.py').read_text())
ns={'np':np};exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'baseline_helpers','exec'),ns)
raw=[];runstatus=[]
for g in priors:
 for sign,mode in [(-1,'absorption'),(1,'emission')]:
  folder=P/'run'/f"espada_{g['id']}_spw0";path=folder/f"{mode}_{g['id']}_spw0_cat.xml"
  rec=dict(prior=g['id'],mode=mode,catalogue_exists=path.exists(),accepted=0)
  if path.exists():
   for r in Table.read(path):
    zz=float(r['z'])+g['channel_start'];fc=float(np.interp(zz,np.arange(len(f)),f));dv=C*(g['observed_GHz']-fc)/g['observed_GHz']
    if abs(dv)>cfg['velocity_halfwidth_kms']:continue
    lo=int(r['z_min'])+g['channel_start'];hi=int(r['z_max'])+g['channel_start'];rec['accepted']+=1
    raw.append(dict(sign=sign,mode=mode,x=float(r['x']),y=float(r['y']),z=zz,frequency_GHz=fc,lo=lo,hi=hi,trigger=g['id'],margin=min(lo-g['channel_start'],g['channel_stop']-hi)))
  runstatus.append(rec)
# Overlapping queries can recover the same source repeatedly. Prefer an untruncated
# realization with most padding; grouping tolerance is 0.024 arcsec and 3 channels.
features=[]
for r in sorted(raw,key=lambda r:r['margin'],reverse=True):
 match=next((q for q in features if q['sign']==r['sign'] and np.hypot(q['x']-r['x'],q['y']-r['y'])<=6 and abs(q['z']-r['z'])<=3),None)
 if match:match['triggered_by'].append(r['trigger'])
 else:features.append(dict(r,triggered_by=[r['trigger']]))
yy,xx=np.mgrid[-150:151,-150:151]*.004;rr=np.hypot(xx,yy)
for i,r in enumerate(sorted(features,key=lambda r:r['frequency_GHz']),1):
 r['id']=f'F{i:02}';x=int(round(r['x']));y=int(round(r['y']));dist=np.hypot(xx-(x-150)*.004,yy-(y-150)*.004)
 masks=[(rr>.25)&(rr<.4)&(dist>.09),(rr>.4)&(rr<.55)&(dist>.09)]
 excluded=np.unique([k for q in features if np.hypot(q['x']-x,q['y']-y)*.004<=.09 for k in range(q['lo']-1,q['hi']+2)])
 models,_=ns['evaluate'](a,r['lo'],r['hi'],x,y,r['sign'],masks,abs(h['CDELT3'])/1000,excluded)
 r['models']=models;valid=all(m['snr'] is not None for m in models)
 r['worst_snr']=float(min(m['snr'] for m in models)) if valid else None
 variation=abs(models[0]['integral']-models[1]['integral'])/max(abs(m['integral']) for m in models) if valid else None
 r['baseline_variation']=variation;r['flags']=[]
 if not all(m['complete'] for m in models):r['flags'].append('baseline-coverage')
 if not valid or min(m['consecutive3'] for m in models)<2:r['flags'].append('channel-support')
 if valid and variation>.3:r['flags'].append('baseline-sensitive')
 if r['margin']<2:r['flags'].append('cutout-edge')
 channel_noise=np.load(OLD/'qa/spectra.npz')['spw0'][5]
 if np.max(channel_noise[r['lo']:r['hi']+1])>3*np.median(channel_noise):r['flags'].append('noisy-channel')
 if r['lo']<2 or r['hi']>=len(a)-2:r['flags'].append('spectral-edge')
 if min(x,y,300-x,300-y)<15:r['flags'].append('spatial-edge')
 r['robust']=bool(valid and r['worst_snr']>=8 and not r['flags'])
 r['alternatives']=[]
 for q in catalog:
  dv=C*(q['observed_GHz']-r['frequency_GHz'])/q['observed_GHz']
  if not q['excluded'] and abs(dv)<=200:r['alternatives'].append(dict(q,velocity_offset_kms=dv))
 r['prior_matches']=[g['id'] for g in priors if abs(C*(g['observed_GHz']-r['frequency_GHz'])/g['observed_GHz'])<=200]
features.sort(key=lambda r:r['id'])
(P/'qa/features.json').write_text(json.dumps(features,indent=2));(P/'qa/run_status.json').write_text(json.dumps(runstatus,indent=2))
lines=['# IC 860 SPW0: catalogue-conditioned search','', '| Prior | Species | Rest GHz | Expected observed GHz | Eu (K) | Candidate features within ±200 km/s |','|---|---|---:|---:|---:|---|']
for g in priors:
 r=g['members'][0];found=[q['id']+(' (robust feature)' if q['robust'] else ' (provisional)') for q in features if g['id'] in q['prior_matches']]
 lines.append(f"| {g['id']} | {r['species']} | {r['rest_GHz']:.7f} | {g['observed_GHz']:.7f} | {r['Eu_K']:.1f} | {', '.join(found) or 'No accepted component'} |")
lines+=['','A feature can match several priors. These are candidate associations, not molecule detections.','', '| Feature | Sign | Observed GHz | Lower integrated S/N | Robust feature? | Flags |','|---|---|---:|---:|---|---|']
for r in features:lines.append(f"| {r['id']} | {r['mode']} | {r['frequency_GHz']:.7f} | {r['worst_snr']:.2f} | {r['robust']} | {', '.join(r['flags'])} |" if r['worst_snr'] is not None else f"| {r['id']} | {r['mode']} | {r['frequency_GHz']:.7f} | N/A | False | {', '.join(r['flags'])} |")
(P/'qa/results.md').write_text('\n'.join(lines)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(15,5));ax.step(f,a[:,150,150]*1000,where='mid',color='#273747',label='Nucleus, 0.060 arcsec beam')
for i,g in enumerate(priors):
 fc=g['observed_GHz']
 if fc<f.min() or fc>f.max():continue
 ax.axvline(fc,color='#368d95',alpha=.5,lw=.8)
 ax.text(fc,2.0+(i%3)*1.6,g['id']+' '+g['members'][0]['species'].split(' v')[0],rotation=75,ha='center',va='bottom',fontsize=8,color='#23727a')
for r in features:
 if r['robust']:
  ax.annotate(r['id']+' absorption',xy=(r['frequency_GHz'],a[int(round(r['z'])),150,150]*1000),xytext=(r['frequency_GHz']-.09,-7.5),arrowprops={'arrowstyle':'->'},fontsize=10)
ax.set(xlabel='Observed LSRK frequency (GHz)',ylabel='Flux density (mJy/beam)',title='IC 860 Band 5 SPW0 — catalogue positions and nuclear spectrum',ylim=(-9,9),xlim=(f.min(),f.max()));ax.axhline(0,color='gray',lw=.6);ax.legend(loc='lower left');fig.tight_layout();fig.savefig(P/'qa/catalogue_overview.png',dpi=160);plt.close(fig)
fig,axes=plt.subplots(6,3,figsize=(16,20))
for ax,g in zip(axes.flat,priors):
 r=g['members'][0];v=C*(g['observed_GHz']-f)/g['observed_GHz'];sel=abs(v)<450
 ax.step(v[sel],a[sel,150,150]*1000,where='mid',color='#273747',label='Nucleus')
 matched=[q for q in features if g['id'] in q['prior_matches']]
 for q in matched:
  dv=C*(g['observed_GHz']-q['frequency_GHz'])/g['observed_GHz'];ax.axvline(dv,color='#b35d30' if q['robust'] else '#aaa',ls=':',lw=1)
  ax.text(dv,.93,q['id'],transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
 ax.axvspan(-200,200,color='#53a6ab',alpha=.10);ax.axvline(0,color='#23727a',lw=.8);ax.axhline(0,color='#aaa',lw=.6)
 ax.set_title(f"{g['id']}  {r['species']}\n{r['rest_GHz']:.6f} GHz; Eu={r['Eu_K']:.1f} K",fontsize=10)
 ax.set(xlim=(-450,450),xlabel='Velocity offset from predicted position (km/s)',ylabel='mJy/beam')
fig.suptitle('Catalogue-conditioned windows: nuclear spectra\nShading: ±200 km/s search; F labels mark spatial components (not necessarily nuclear)',fontsize=14);fig.tight_layout(rect=(0,0,1,.96));fig.savefig(P/'qa/targeted_windows.png',dpi=130);plt.close(fig)
print('Accepted window-components:',len(raw),'unique:',len(features),'robust:',sum(r['robust'] for r in features))
for r in features:print(r['id'],r['mode'],r['frequency_GHz'],r['worst_snr'],r['robust'],r['prior_matches'])
