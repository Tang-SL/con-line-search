"""Prepare full-field IC 860 examples without changing the original FITS files."""
import argparse
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from conlines.core import read_cube


def prepare(source,output,spws):
    output.mkdir(parents=True,exist_ok=True)
    with fits.open(source/'IC860_cont_natural_new_spw4cont.image.fits') as hd:
        image=hd[0].data.squeeze();y,x=np.unravel_index(np.nanargmax(image),image.shape)
        centre=WCS(hd[0].header).celestial.pixel_to_world(x,y)
    names={0:'IC860-spw0_H2S_20kms_natural.image.fits',1:'IC860-spw1_20kms_natural.image.fits',2:'IC860-spw2_20kms_natural.image.fits',3:'IC860-spw3_H2O_20kms_natural.image.fits'}
    for spw in spws:
        target=output/f'spw{spw}_full_beam060.fits'
        if target.exists():raise ValueError(f'Refusing to overwrite {target}')
        cube=read_cube(source/names[spw],target_beam_arcsec=.06)
        h=WCS(cube.header).sub([1,2,3]).to_header()
        for key in ['BUNIT','BMAJ','BMIN','BPA','OBJECT','SPECSYS','RESTFRQ','RESTFREQ']:
            if key in cube.header:h[key]=cube.header[key]
        cx,cy=WCS(h).celestial.world_to_pixel(centre)
        h['NUCX']=(float(cx),'Continuum peak x pixel (zero-based)')
        h['NUCY']=(float(cy),'Continuum peak y pixel (zero-based)')
        h['ORIGFILE']=names[spw]
        h.add_history('Full original spatial field retained; no spatial cutout.')
        h.add_history('Each channel convolved to circular 0.060 arcsec beam with beam-area scaling.')
        h.add_history('Convolution zero-pads at image boundary; avoid edge apertures.')
        h.add_history('No additional primary-beam or continuum correction applied.')
        fits.writeto(target,cube.data,h,checksum=True)
        print(f'SPW{spw}: {cube.data.shape} -> {target}',flush=True)
        del cube

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path,help='Original Sabine/b5 directory')
    parser.add_argument('--output',type=Path,default=ROOT/'targets/ic860/alma/fullfield/inputs')
    parser.add_argument('--spws',nargs='+',type=int,choices=range(4),default=list(range(4)))
    args=parser.parse_args();prepare(args.source,args.output,args.spws)
