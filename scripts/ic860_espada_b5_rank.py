"""Rank all four-SPW ESPADA candidates using predeclared selection_rules.md."""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'targets/ic860/alma/espada_b5_all'
Q=P/'qa'
ss=np.load(Q/'spectra.npz')
yy,xx=np.mgrid[-150:151,-150:151]*.004
rr=np.hypot(xx,yy)
rows=[]; display={}


def stats(image, masks):
    """Noise estimates use untrimmed real background and no independence claim."""
    scales=[]
    for mask in masks:
        v=image[...,mask]
        med=np.median(v,axis=-1)
        scales.append(1.4826*np.median(abs(v-med[...,None]),axis=-1))
        scales.append(np.std(v,axis=-1,ddof=1))
    return np.median(image[...,masks[0]|masks[1]],axis=-1),np.maximum.reduce(scales)


def longest_run(bools):
    best=now=0
    for b in bools:
        now=now+1 if b else 0; best=max(best,now)
    return best


def evaluate(a, lo, hi, x, y, sign, masks, dv, excluded):
    mid=(lo+hi)/2
    sel=np.arange(lo,hi+1)
    outcomes=[]; plots=[]
    for model,width in [('A',5),('B',10)]:
        # Nearest available clean channels, searched at most 30 channels out.
        # This avoids turning a line-rich spectrum into an arbitrary edge case
        # solely because another detected line occupies a fixed sideband.
        left=np.arange(lo-3,max(-1,lo-31),-1)
        right=np.arange(hi+3,min(len(a),hi+31))
        left=left[~np.isin(left,excluded)][:width]
        right=right[~np.isin(right,excluded)][:width]
        sides=np.r_[left,right]
        complete=bool(len(left)==width and len(right)==width)
        enough=bool(np.sum(sides<lo)>=3 and np.sum(sides>hi)>=3)
        if not enough:
            outcomes.append({'model':model,'complete':False,'snr':None,'integral':None,
                             'noise':None,'consecutive3':0,'sideband_channels':sides.tolist()})
            plots.append(None); continue
        if model=='A':
            intercept=np.median(a[sides],axis=0); slope=np.zeros_like(intercept)
        else:
            z=sides-mid; meanz=z.mean(); av=a[sides].mean(axis=0)
            slope=np.einsum('i,ijk->jk',z-meanz,a[sides])/np.sum((z-meanz)**2)
            intercept=av-meanz*slope
        residual=a[sel]-intercept-(sel-mid)[:,None,None]*slope
        integ=residual.sum(axis=0)*dv
        background,noise=stats(integ,masks)
        sig=float(integ[y,x]-background)
        bgch,nch=stats(residual,masks)
        consecutive=longest_run(sign*(residual[:,y,x]-bgch)/nch>=3)
        outcomes.append(dict(model=model,complete=complete,snr=sign*sig/float(noise),
            integral=sig,noise=float(noise),consecutive3=int(consecutive),sideband_channels=sides.tolist()))
        plots.append(dict(map=integ-background,baseline=intercept[y,x]+(np.arange(len(a))-mid)*slope[y,x],
                          noise=float(noise)))
    return outcomes,plots


for spw in range(4):
    a=fits.getdata(P/'inputs'/f'spw{spw}_beam060.fits')
    paired=fits.getdata(P/'inputs'/f'spw{spw}_pluscont_beam060.fits')
    h=fits.getheader(P/'inputs'/f'spw{spw}_beam060.fits'); dv=abs(h['CDELT3'])/1000
    freq=ss[f'spw{spw}'][1]; channel_noise=ss[f'spw{spw}'][5]
    neighbours=[]
    for cat in (P/'run'/f'espada_spw{spw}_beam060').glob('*cat.xml'):
        neighbours.extend(list(Table.read(cat)))
    for mode in ['absorption','emission']:
        p=P/'run'/f'espada_spw{spw}_beam060'/f'{mode}_spw{spw}_beam060_cat.xml'
        if not p.exists(): continue
        for r in Table.read(p):
            ident=f"S{spw}-{'A' if mode=='absorption' else 'E'}{int(r['id'])}"
            x=int(round(float(r['x']))); y=int(round(float(r['y'])))
            lo=int(r['z_min']); hi=int(r['z_max']); sign=-1 if mode=='absorption' else 1
            distance=np.hypot(xx-(x-150)*.004,yy-(y-150)*.004)
            masks=[(rr>.25)&(rr<.40)&(distance>.09),(rr>.40)&(rr<.55)&(distance>.09)]
            # Apply to every candidate: exclude detected neighbours at <=1.5 beams
            # from the baseline windows, with a one-channel guard on each side.
            excluded=[]
            for other in neighbours:
                if np.hypot(float(other['x'])-x,float(other['y'])-y)*.004<=.09:
                    excluded.extend(range(int(other['z_min'])-1,int(other['z_max'])+2))
            excluded=np.unique(excluded)
            models,plots=evaluate(a,lo,hi,x,y,sign,masks,dv,excluded)
            pairmodels,_=evaluate(paired,lo,hi,x,y,sign,masks,dv,excluded)
            valid=all(m['snr'] is not None for m in models)
            score=min(m['snr'] for m in models) if valid else -999.
            integrals=[m['integral'] for m in models]
            variation=abs(integrals[0]-integrals[1])/max(abs(v) for v in integrals) if valid and any(integrals) else 1.
            flags=[]
            if not valid: flags.append('insufficient-clean-sidebands')
            if not all(m['complete'] for m in models): flags.append('baseline-coverage')
            if min(x,y,300-x,300-y)<15: flags.append('spatial-edge')
            if lo<2 or hi>=len(a)-2: flags.append('spectral-edge')
            if np.max(channel_noise[lo:hi+1])>3*np.median(channel_noise): flags.append('noisy-channel')
            if valid and min(m['consecutive3'] for m in models)<2: flags.append('channel-support')
            if valid and variation>.30: flags.append('baseline-sensitive')
            tier='strict' if score>=8 and not flags else ('follow-up' if score>=5 and not flags else 'flagged/weak')
            z=float(r['z']); vel=(h['CRVAL3']+(z+1-h['CRPIX3'])*h['CDELT3'])/1000
            fc=h['RESTFRQ']*(1-vel/299792.458)/1e9
            row=dict(key=ident,spw=spw,mode=mode,catalog_id=int(r['id']),
              frequency_GHz=float(fc),x=x,y=y,east_arcsec=-(x-150)*.004,north_arcsec=(y-150)*.004,
              channel_min=lo,channel_max=hi,channels=hi-lo+1,worst_snr=float(score) if valid else None,
              baseline_variation=float(variation) if valid else None,tier=tier,flags=flags,models=models,
              paired_models=pairmodels,excluded_baseline_channels=excluded.tolist())
            rows.append(row)
            if valid:
                worse=int(np.argmin([m['snr'] for m in models]))
                display[ident]=dict(freq=freq,profile=a[:,y,x]*1000,noise=channel_noise*1000,
                    plots=plots,selected=plots[worse],lo=lo,hi=hi)
    print(f'SPW{spw} assessed',flush=True)

rows.sort(key=lambda r:r['worst_snr'] if r['worst_snr'] is not None else -np.inf,reverse=True)
for i,r in enumerate(rows,1):r['rank']=i
(Q/'ranked_candidates.json').write_text(json.dumps(rows,indent=2)+'\n')

text=['# Band 5: all four SPWs — ranked candidate components','',
 'Ranked by the smaller baseline-corrected integrated S/N at each candidate’s',
 'fixed catalogue centroid. This is an operational robustness ranking, not a',
 'formal detection probability or a molecular-line identification.','',
 '| Rank | ID | SPW | Sign | LSRK GHz | S/N A | S/N B | Worst S/N | Baseline change | Class | Flags |',
 '|---:|---|---:|---|---:|---:|---:|---:|---:|---|---|']
for r in rows:
    nums=[f"{m['snr']:.1f}" if m['snr'] is not None else 'N/A' for m in r['models']]
    score_text=f"{r['worst_snr']:.1f}" if r['worst_snr'] is not None else 'N/A'
    change_text=f"{100*r['baseline_variation']:.0f}%" if r['baseline_variation'] is not None else 'N/A'
    text.append(f"| {r['rank']} | {r['key']} | {r['spw']} | {r['mode']} | {r['frequency_GHz']:.6f} | {nums[0]} | {nums[1]} | {score_text} | {change_text} | {r['tier']} | {', '.join(r['flags']) or '—'} |")
(Q/'ranked_candidates.md').write_text('\n'.join(text)+'\n')

colors={'strict':'#176c66','follow-up':'#ca8626','flagged/weak':'#a6a9af'}
fig,ax=plt.subplots(figsize=(12,max(6,len(rows)*.28)),layout='constrained')
positions=np.arange(len(rows))
ax.barh(positions,[r['worst_snr'] if r['worst_snr'] is not None else 0 for r in rows],color=[colors[r['tier']] for r in rows])
for i,r in enumerate(rows):
    if r['worst_snr'] is None: ax.text(.2,i,'Not assessable: insufficient clean sidebands',va='center',fontsize=8)
ax.set_yticks(positions,[f"{r['key']} · {r['frequency_GHz']:.4f} GHz" for r in rows],fontsize=9)
ax.invert_yaxis(); ax.axvline(8,color='#176c66',ls='--',label='Strict S/N threshold: 8')
ax.axvline(5,color='#ca8626',ls=':',label='Follow-up S/N threshold: 5')
ax.set_xlabel('Lower integrated S/N across two baseline models')
ax.set_title('IC 860 — all Band 5 candidates ranked by robustness\nGreen: passes all strict rules · orange: follow-up · grey: weak or quality-flagged')
ax.legend(loc='lower right'); ax.grid(axis='x',alpha=.15)
fig.savefig(Q/'ranking.png',dpi=160); plt.close(fig)

# Show every strict candidate; if none pass, show the six best-scoring candidates.
chosen=[r for r in rows if r['tier']=='strict'] or rows[:6]
for old in Q.glob('shortlist_[0-9]*.png'):
    old.unlink()  # Remove only this script's stale derived shortlist pages.
for start in range(0,len(chosen),4):
    subset=chosen[start:start+4]; n=len(subset)
    fig,axes=plt.subplots(2,n,figsize=(max(8,4.6*n),8),squeeze=False,layout='constrained')
    for j,r in enumerate(subset):
        d=display[r['key']]; f=d['freq'];lo=d['lo'];hi=d['hi']
        low,high=sorted([f[lo],f[hi]])
        ax=axes[0,j]
        for k,model in enumerate(d['plots']):
            ax.plot(f,model['baseline']*1000,color=['#ad538a','#d68620'][k],ls=['--',':'][k],label=f'Baseline {"AB"[k]}')
        ax.step(f,d['profile'],where='mid',color='#19628b',label='At catalogue centroid')
        ax.axvspan(low,high,color='#ddba55',alpha=.2)
        ax.axhline(0,color='.5',lw=.7)
        ax.set(xlim=(low-.11,high+.11),xlabel='LSRK frequency (GHz)',ylabel='mJy / 0.06″ beam',
          title=f"{r['key']} · {r['mode']}\nWorst S/N {r['worst_snr']:.1f} · baseline change {100*r['baseline_variation']:.0f}%")
        window=(f>low-.11)&(f<high+.11); mn=float(d['profile'][window].min());mx=float(d['profile'][window].max())
        ax.set_ylim(mn-(mx-mn)*.15,mx+(mx-mn)*.6)
        ax.ticklabel_format(axis='x',useOffset=False);ax.legend(fontsize=8,loc='upper left')
        ax=axes[1,j];m=d['selected']['map'];noise=d['selected']['noise']
        region=rr<.3;lim=np.max(abs(m[region]))
        im=ax.imshow(m,origin='lower',extent=[.602,-.602,-.602,.602],cmap='RdBu_r',vmin=-lim,vmax=lim)
        ax.plot(0,0,'+',color='#26b295',ms=10)
        ax.plot(r['east_arcsec'],r['north_arcsec'],'x',color='black',ms=7)
        ax.add_patch(Circle((.22,-.22),.03,fc='white',ec='black'))
        ax.set(xlim=(.3,-.3),ylim=(-.3,.3),xlabel='East offset (″)',ylabel='North offset (″)')
        cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.05,pad=.1);cb.set_label('Jy beam⁻¹ km s⁻¹')
    fig.suptitle('IC 860 — strict Band 5 shortlist\nMaps use the lower-S/N baseline model; + nucleus, × measured position',fontsize=14)
    fig.savefig(Q/f'shortlist_{start//4+1}.png',dpi=160);plt.close(fig)
print(json.dumps([{k:r[k] for k in ['key','frequency_GHz','worst_snr','baseline_variation','tier','flags']} for r in rows],indent=2))
