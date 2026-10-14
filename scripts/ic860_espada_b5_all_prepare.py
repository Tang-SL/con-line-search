"""Prepare a reproducible ESPADA search on all four Sabine Band 5 SPWs.

Original FITS files remain untouched. Detection pilot, not PB-corrected photometry.
Dependencies: astropy, numpy, scipy, radio-beam, matplotlib, espada-astro.
"""
from pathlib import Path
import json
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy import units as u
from radio_beam import Beam
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml
import adplib

warnings.filterwarnings('ignore', category=FITSFixedWarning)
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'targets/ic860/alma/Sabine/b5'
OUT = ROOT/'targets/ic860/alma/espada_b5_all'
for sub in ['inputs', 'config', 'qa']: (OUT/sub).mkdir(parents=True,exist_ok=True)
TARGET = Beam(major=.06*u.arcsec, minor=.06*u.arcsec, pa=0*u.deg)
PIX = .004
HALF = 150  # 1.204 arcsec wide; padded during convolution
PAD = 40
with fits.open(SOURCE/'IC860_cont_natural_new_spw4cont.image.fits') as hd:
    cont = hd[0].data.squeeze()
    cy,cx = np.unravel_index(np.nanargmax(cont),cont.shape)
    center = WCS(hd[0].header).celestial.pixel_to_world(cx,cy)


def prepared(path):
    with fits.open(path,memmap=True) as hd:
        h = hd[0].header
        x,y = WCS(h).celestial.world_to_pixel(center)
        ix,iy = int(round(float(x))),int(round(float(y)))
        n = HALF+PAD
        a = hd[0].data.squeeze()[:,iy-n:iy+n+1,ix-n:ix+n+1]
        assert np.isfinite(a).all(), f'Invalid input pixels: {path}'
        beams=hd['BEAMS'].data
        out=[]
        for plane,b in zip(a,beams):
            native=Beam(b['BMAJ']*u.arcsec,b['BMIN']*u.arcsec,b['BPA']*u.deg)
            kernel=TARGET.deconvolve(native).as_kernel(PIX*u.arcsec).array
            # Normalized convolution, then Jy/native-beam -> Jy/target-beam.
            sm=fftconvolve(plane,kernel/kernel.sum(),mode='same')*(TARGET.sr/native.sr).value
            out.append(sm[PAD:-PAD,PAD:-PAD])
        out=np.asarray(out,dtype='float32')
        w=WCS(h).sub([1,2,3])
        w.wcs.crpix[:2]-=[ix-HALF,iy-HALF]
        oh=w.to_header()
        for k in ['OBJECT','TELESCOP','DATE-OBS','OBSERVER','ORIGIN']: oh[k]=h[k]
        oh['BUNIT']='Jy/beam'
        oh['BMAJ']=.06/3600; oh['BMIN']=.06/3600; oh['BPA']=0.
        oh['PBSTAT']='UNKNOWN'
        oh.add_history('Derived from Sabine/'+path.name)
        oh.add_history('Convolved channel-by-channel to 0.060 arcsec circular beam.')
        oh.add_history('Beam-area scaling applied; central 1.204 arcsec cutout.')
        oh.add_history('Primary beam status unknown; no additional PB correction applied.')
        oh.add_history('No additional continuum or background subtraction applied.')
        vel=(h['CRVAL3']+(np.arange(len(out))+1-h['CRPIX3'])*h['CDELT3'])/1000
        freq=h['RESTFRQ']*(1-vel/299792.458)/1e9
        return out,oh,vel,freq


summary={'center':center.to_string('hmsdms'),'beam_arcsec':.06,
 'cutout_arcsec':(2*HALF+1)*PIX,'PB_status':'unknown; no new correction',
 'noise':'1.4826*MAD in annulus 0.35–0.55 arcsec, each channel after smoothing',
 'continuum_subtraction':'no additional subtraction; compare paired input products',
 'datasets':{}}
fig,axs=plt.subplots(3,4,figsize=(20,10),layout='constrained')
yy,xx=np.mgrid[-HALF:HALF+1,-HALF:HALF+1]*PIX
r=np.hypot(xx,yy); off=(r>.35)&(r<.55); aperture=r<=.05
arrays={}
for j,(spw,label) in enumerate([(0,'H2S'),(1,''),(2,''),(3,'H2O')]):
    stem=f'IC860-spw{spw}'+(f'_{label}' if label else '')
    a,h,v,f=prepared(SOURCE/f'{stem}_20kms_natural.image.fits')
    b,hplus,_,_=prepared(SOURCE/f'{stem}+cont_20kms_natural.image.fits')
    fits.PrimaryHDU(b,hplus).writeto(OUT/'inputs'/f'spw{spw}_pluscont_beam060.fits',overwrite=True,checksum=True)
    fits.PrimaryHDU(a,h).writeto(OUT/'inputs'/f'spw{spw}_beam060.fits',overwrite=True,checksum=True)
    med=np.median(a[:,off],axis=1)
    noise=1.4826*np.median(np.abs(a[:,off]-med[:,None]),axis=1)
    noise_b=1.4826*np.median(np.abs(b[:,off]-np.median(b[:,off],axis=1)[:,None]),axis=1)
    core=a[:,HALF,HALF]; core_b=b[:,HALF,HALF]
    mean=a[:,aperture].mean(axis=1)
    delta=(b-a)[:,aperture].mean(axis=1)
    # This pilot retains all channels. Noise is scaled spectrally by SoFiA.
    summary['datasets'][f'spw{spw}']={
       'source_file':f'{stem}_20kms_natural.image.fits','nchan':len(v),
       'noise_mJy_beam_min_median_max':(1000*np.array([noise.min(),np.median(noise),noise.max()])).tolist(),
       'aperture_pair_difference_mJy_beam_p16_median_p84':(1000*np.percentile(delta,[16,50,84])).tolist(),
       'core_SNR_min_max':[float((core/noise).min()),float((core/noise).max())]}
    arrays[f'spw{spw}']=np.array([v,f,core,core_b,mean,noise,noise_b,med,delta])
    axs[0,j].step(f,core_b*1000,where='mid',color='#c77c20',lw=1,label='+cont file')
    axs[0,j].step(f,core*1000,where='mid',color='#156693',lw=1,label='without +cont')
    axs[0,j].fill_between(f,-3*noise*1000,3*noise*1000,color='.7',alpha=.25,label='±3 × off-source noise')
    axs[0,j].set_title(f'Band 5 SPW{spw} · {label} window')
    axs[0,j].set_ylabel('At continuum peak (mJy / beam)')
    axs[0,j].legend(fontsize=8)
    axs[1,j].plot(f,noise*1000,label='without +cont',color='#156693')
    axs[1,j].plot(f,noise_b*1000,label='+cont',color='#c77c20',ls='--')
    axs[1,j].set_ylabel('Channel noise (mJy / beam)')
    axs[1,j].legend(fontsize=8)
    axs[2,j].step(f,core/noise,where='mid',color='#156693')
    axs[2,j].axhspan(-3,3,color='.7',alpha=.25)
    axs[2,j].set_ylabel('Peak-pixel brightness / channel noise')
    for ax in axs[:,j]:
        ax.axhline(0,color='.4',lw=.7); ax.set_xlabel('LSRK frequency (GHz)')
        ax.set_xlim(f.min(),f.max()); ax.ticklabel_format(axis='x',useOffset=False)
fig.suptitle('IC 860 Band 5 — pre-search checks at a common 0.06″ beam\nNo new continuum subtraction; noise from a 0.35–0.55″ annulus',fontsize=15)
fig.savefig(OUT/'qa'/'preflight.png',dpi=170)
np.savez_compressed(OUT/'qa'/'spectra.npz',**arrays)
(OUT/'qa'/'preflight.json').write_text(json.dumps(summary,indent=2)+'\n')

# Preserve full parameter files and explicit overrides for reproducibility.
package=Path(adplib.__file__).parent
overrides={
 'scfind.kernelsXY':'0, 15', 'scfind.kernelsZ':'0, 3, 5',
 'scfind.threshold':'5.0','scfind.fluxRange':'both',
 'scaleNoise.enable':'true','scaleNoise.mode':'spectral','scaleNoise.fluxRange':'both',
 'linker.radiusXY':'3','linker.radiusZ':'1','linker.minSizeXY':'5',
 'linker.minSizeZ':'2','linker.minPixels':'50','linker.positivity':'true',
 'reliability.enable':'false','reliability.plot':'false',
 'filter.minSNR':'0','filter.minPixels':'0','filter.discardNegative':'true',
 'dilation.enable':'false','output.overwrite':'true','output.writePV':'true',
 'output.writeNoise':'true','output.writeDiagnosticPlot':'false'}
# Reliability based on positive/negative symmetry is inappropriate with real
# absorption present and few sources. Outputs are candidates, not validated IDs.
for mode in ['abs','emi']:
    lines=[]
    for line in (package/'sofia'/f'sofia_{mode}_default.par').read_text().splitlines():
        key=line.split('=')[0].strip()
        lines.append(f'{key} = {overrides[key]}' if key in overrides else line)
    (OUT/'config'/f'{mode}.par').write_text('\n'.join(lines)+'\n')
cfg=yaml.safe_load((package/'config.yaml').read_text())
cfg.update(make_report=True, verbose=True,num_cores=1,clear_logs=False,
 output_dir=str(OUT/'run'),input_data_set={
 'spw0':[str(OUT/'inputs'/'spw0_beam060.fits')],
 'spw3':[str(OUT/'inputs'/'spw3_beam060.fits')],
 'spw1':[str(OUT/'inputs'/'spw1_beam060.fits')],
 'spw2':[str(OUT/'inputs'/'spw2_beam060.fits')]},
 enable_tap_service=False,enable_sofia=True,run_mode='both',use_mask=False,
 abs_flag_cube=False,auto_setup=False,enable_sip=False,enable_group=False,
 sofia_abs_file=str(OUT/'config'/'abs.par'),sofia_emi_file=str(OUT/'config'/'emi.par'))
# SIP disabled for this mixed-line pilot: avoid assigning one molecule to every
# candidate. ESPADA still creates masks/moments/catalogues and its HTML report.
(OUT/'config'/'pilot.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
print(json.dumps(summary,indent=2))
