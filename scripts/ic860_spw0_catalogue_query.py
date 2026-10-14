import requests,json
from pathlib import Path
from astropy.io import fits
import numpy as np
out=Path(__file__).resolve().parents[1]/'targets/ic860/alma/catalogue_spw0'
h=fits.getheader(out.parent/'espada_b5_all/inputs/spw0_beam060.fits')
settings=json.loads((out/'config/settings.json').read_text())
c=299792.458;z=settings['z'];width=settings['velocity_halfwidth_kms']
v=(h['CRVAL3']+(np.arange(h['NAXIS3'])+1-h['CRPIX3'])*h['CDELT3'])/1000
f=h['RESTFRQ']/1e9*(1-v/c)
lo=f.min()*(1+z)*(1-width/c);hi=f.max()*(1+z)*(1+width/c)
print('Observed',f.min(),f.max(),'Rest query',lo,hi,flush=True)
for filt in ['ism-hot-core','extragalactic']:
 p=dict(searchSpeciesAutocomplete='',quickPickerMolecules=[],telescopeBands=['any'],energyFrom='',energyTo='',energyValue='el-1',frequenciesFrom=[str(lo)],frequenciesTo=[str(hi)],frequenciesUnit='GHz',frequencyRedshift=0,astronomicalFilters=filt,exportType='current',exportDelimiter='tab',exportLimit='allRecords',exportStart=1,exportStop=250)
 (out/'catalogue'/f'{filt}_request.json').write_text(json.dumps(p,indent=2))
 r=requests.post('https://splatalogue.online/splata-slap/basicsearch/false',json={'body':json.dumps(p)},timeout=90);r.raise_for_status()
 (out/'catalogue'/f'{filt}_raw.json').write_text(r.text)
 data=r.json();print(filt,len(data),str(data[:1])[:1800],str(data[-1:])[:600],flush=True)
