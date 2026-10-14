"""Catalogue-conditioned ESPADA pilot; settings and raw catalogue are saved."""
from pathlib import Path
import json,re,html,hashlib,datetime
import numpy as np
from astropy.io import fits
import yaml
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'targets/ic860/alma/catalogue_spw0'; OLD=P.parent/'espada_b5_all'
for d in ['inputs','config','qa','run']: (P/d).mkdir(exist_ok=True)
C=299792.458; z=3880/C
cfg=dict(z=z,cz_optical_kms=3880,systemic_uncertainty_kms=20,velocity_halfwidth_kms=200,upper_energy_max_K=300,padding_channels=5,primary_filter='extragalactic',alternative_filter='ism-hot-core',source='https://arxiv.org/abs/1905.07275',frame_note='Published cz adopted as working redshift; paper does not explicitly identify reference frame in systemic-velocity paragraph. Cube is LSRK; no unverified frame correction applied.',exclusions=['deuterated species (uppercase D atom)','trans conformers'],catalogue_query_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
if (P/'config/settings.json').exists():
 cfg.update(json.loads((P/'config/settings.json').read_text()))
z=cfg['z']
(P/'config/settings.json').write_text(json.dumps(cfg,indent=2))
clean=lambda x:html.unescape(re.sub('<[^>]*>','',x or '')).strip()
def number(v):
 try:return float(v)
 except (ValueError,TypeError):return None
rows={}
for source in ['ism-hot-core','extragalactic']:
 for r in json.loads((P/'catalogue'/f'{source}_raw.json').read_text()):
  if not r['species_id']:continue
  name=clean(r['species_name']);freq=float(r['orderedFreq'].split(',')[0]);qn=clean(r['resolved_QNs'])
  key=(r['species_id'],freq,qn)
  if key in rows:rows[key]['filters'].append(source);continue
  reasons=[]
  if re.search(r'D(?![a-z])',name):reasons.append('deuterated')
  if re.search(r'\btrans\b|^t-',name,re.I):reasons.append('trans conformer')
  rows[key]=dict(species=name,species_id=r['species_id'],rest_GHz=freq,observed_GHz=freq/(1+z),qns=qn,Eu_K=number(r['euk']),log10_A_s=number(r['aij']),catalogue=r['linelists'],filters=[source],excluded=reasons)
allrows=sorted(rows.values(),key=lambda r:r['rest_GHz'])
(P/'catalogue/normalized.json').write_text(json.dumps(allrows,indent=2))
a,h=fits.getdata(OLD/'inputs/spw0_beam060.fits',header=True)
f=h['RESTFRQ']/1e9*(1-(h['CRVAL3']+(np.arange(len(a))+1-h['CRPIX3'])*h['CDELT3'])/1000/C)
priors=[r for r in allrows if not r['excluded'] and 'extragalactic' in r['filters'] and r['Eu_K'] is not None and r['Eu_K']<=cfg['upper_energy_max_K']]
# Combine unresolved entries of same species within one channel; preserve members.
groups=[]
for r in priors:
 if groups and r['species_id']==groups[-1]['members'][0]['species_id'] and abs(r['observed_GHz']-groups[-1]['members'][0]['observed_GHz'])<abs(f[1]-f[0]):groups[-1]['members'].append(r)
 else:groups.append(dict(members=[r]))
base=yaml.safe_load((OLD/'config/pilot.template.yaml').read_text());base['input_data_set']={};base['output_dir']=str(P/'run');base['make_report']=False
coverage=np.zeros(len(f),bool)
for i,g in enumerate(groups,1):
 g['id']=f'L{i:02}';g['observed_GHz']=float(np.mean([r['observed_GHz'] for r in g['members']]))
 vel=C*(g['observed_GHz']-f)/g['observed_GHz'];idx=np.where(abs(vel)<=cfg['velocity_halfwidth_kms'])[0]
 if not len(idx):continue
 coverage[idx]=True;lo=max(0,int(idx.min())-cfg['padding_channels']);hi=min(len(a)-1,int(idx.max())+cfg['padding_channels'])
 g.update(channel_start=lo,channel_stop=hi,search_channels=idx.tolist(),edge_truncated=bool(idx.min()==0 or idx.max()==len(a)-1))
 hh=h.copy();hh['CRPIX3']-=lo;hh['HISTORY']='Catalogue-prior cutout; detection is active, no forced source mask.'
 dest=P/'inputs'/f"{g['id']}_spw0.fits";fits.writeto(dest,a[lo:hi+1],hh,overwrite=True,checksum=True)
 base['input_data_set'][g['id']]=[str(dest)]
for sign in ['abs','emi']:
 dest=P/'config'/f'{sign}.par';dest.write_text((OLD/'config'/f'{sign}.par').read_text());base[f'sofia_{sign}_file']=str(dest)
(P/'config/espada.yaml').write_text(yaml.safe_dump(base,sort_keys=False))
(P/'catalogue/priors.json').write_text(json.dumps(groups,indent=2))
summary=dict(observed_GHz=[float(f.min()),float(f.max())],normalized_rows=len(allrows),excluded=sum(bool(r['excluded']) for r in allrows),primary_entries=len(priors),groups=len(groups),target_channels=int(coverage.sum()),total_channels=len(f),coverage_fraction=float(coverage.mean()))
(P/'qa/preparation.json').write_text(json.dumps(summary,indent=2));print(summary)
