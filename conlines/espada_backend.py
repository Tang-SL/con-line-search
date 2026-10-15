"""Run the real ESPADA/SoFiA 3D finder, independently of line catalogues."""
from pathlib import Path
import os, shutil, subprocess, tempfile, json, re
import numpy as np
import yaml
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from .core import assess, C

ROOT=Path(__file__).resolve().parents[1]


def executables():
    import sys
    espada=os.environ.get('ESPADA_EXECUTABLE') or shutil.which('espada')
    if not espada:
        candidate=Path(sys.executable).parent/'espada'
        if candidate.is_file():espada=str(candidate)
    sofia=os.environ.get('SOFIA_EXECUTABLE') or shutil.which('sofia')
    if not espada or not sofia:
        raise ValueError('ESPADA and SoFiA are required. Install espada-astro and put sofia on PATH (or set SOFIA_EXECUTABLE), then restart the app.')
    return espada,sofia


def run_espada(cube, widths=(3,5,9), mode='Both', clip=5., continuum_subtracted=True):
    """Search the full spatial/spectral cube; cacheable without aperture/catalogue."""
    espada,sofia=executables()
    if not np.all(np.isfinite(cube.data)):
        raise ValueError('This ESPADA integration requires finite pixels. Crop blank borders or prepare a finite cube before detection.')
    modes=['absorption','emission'] if mode=='Both' else [mode.lower()]
    with tempfile.TemporaryDirectory(prefix='conlines-espada-') as folder:
        work=Path(folder);h=WCS(cube.header).sub([1,2,3]).to_header()
        for k in ['BUNIT','BMAJ','BMIN','BPA','RESTFRQ','RESTFREQ','OBJECT','SPECSYS']:
            if k in cube.header:h[k]=cube.header[k]
        source=work/'cube.fits';fits.writeto(source,cube.data,h,checksum=True)
        cfg=yaml.safe_load((ROOT/'targets/ic860/alma/espada_b5_all/config/pilot.template.yaml').read_text())
        cfg.update(make_report=False,output_dir=str(work/'run'),input_data_set={'cube':[str(source)]},
                   run_mode='both' if mode=='Both' else mode.lower(),num_cores=1)
        # One extra beam-wide spatial smoothing, in pixels; preserve all channels.
        spatial=max(1,round(cube.beam_arcsec/max(cube.pixel_arcsec)))
        parameters={}
        for label in ['abs','emi']:
            overrides={'scfind.threshold':str(clip),'scfind.kernelsXY':f'0, {spatial}',
                       'scfind.kernelsZ':', '.join(map(str,[0,*widths])),
                       'contsub.enable':'false' if continuum_subtracted else 'true',
                       'contsub.order':'1','output.writeNoise':'false'}
            text=(ROOT/f'targets/ic860/alma/espada_b5_all/config/{label}.par').read_text()
            lines=[f'{line.split("=")[0].strip()} = {overrides[line.split("=")[0].strip()]}' if line.split('=')[0].strip() in overrides else line for line in text.splitlines()]
            path=work/f'{label}.par';path.write_text('\n'.join(lines)+'\n');cfg[f'sofia_{label}_file']=str(path)
            parameters[label]=path.read_text()
        config=work/'espada.yaml';config.write_text(yaml.safe_dump(cfg))
        env=dict(os.environ,PATH=str(Path(sofia).resolve().parent)+os.pathsep+os.environ.get('PATH',''),MPLBACKEND='Agg',MPLCONFIGDIR=str(work/'mpl'),XDG_CACHE_HOME=str(work/'cache'))
        logpath=work/'console.log'
        with logpath.open('w') as stream:
            try:proc=subprocess.run([espada,'-c',str(config)],cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1200,check=False)
            except subprocess.TimeoutExpired as e:raise ValueError('ESPADA exceeded the 20-minute run limit. Use a smaller spatial cube.') from e
        log=logpath.read_text(errors='replace')
        if proc.returncode:raise ValueError('ESPADA failed: '+log[-1800:])
        failures=re.findall(r'failed with exit code (\d+)',log)
        if any(code!='8' for code in failures):raise ValueError('SoFiA failed: '+log[-1800:])
        components=[];status={}
        for sign in modes:
            directory=work/'run/espada_cube';cat=directory/f'{sign}_cube_cat.xml';maskpath=directory/f'{sign}_cube_mask.fits'
            if not cat.exists():
                # ESPADA can exit zero even if a child fails: explicitly verify no-source status.
                pattern=rf'failed with exit code 8:[^\n]*Mode: {sign}'
                if not re.search(pattern,log):raise ValueError(f'ESPADA produced no {sign} catalogue and did not report a verified no-source outcome.')
                status[sign]='No sources (SoFiA exit 8)';continue
            if not maskpath.exists():raise ValueError(f'Missing {sign} source mask.')
            mask=fits.getdata(maskpath);table=Table.read(cat);status[sign]=f'{len(table)} components'
            for r in table:
                lo,hi=int(r['z_min']),int(r['z_max']);source_id=int(r['id'])
                footprint=np.flatnonzero(np.any(mask[lo:hi+1]==source_id,axis=0)).tolist()
                fc=float(np.interp(float(r['z']),np.arange(len(cube.frequency)),cube.frequency))
                components.append(dict(sofia_id=source_id,sign=sign.capitalize(),lo=lo,hi=hi,
                                       x=float(r['x']),y=float(r['y']),frequency_GHz=fc,footprint=footprint))
        # No spectral catalogue, redshift or molecular information enters this result.
        return dict(components=components,status=status,parameters=parameters,console_log=log,
                    engine='ESPADA/SoFiA 3D',settings=dict(clip=clip,widths=list(widths),mode=mode,continuum_subtracted=continuum_subtracted))


def aperture_features(cube,extraction,result,threshold=8):
    from .measurements import measure_feature
    features=[];selected=set(np.flatnonzero(extraction['mask']).tolist())
    for component in result['components']:
        if not selected.intersection(component['footprint']):continue
        q={k:v for k,v in component.items() if k!='footprint'}
        flags=[]
        try:
            score,models,variation=assess(extraction,q['lo'],q['hi'])
            sign=1 if q['sign']=='Emission' else -1
            if any(sign*m['signal']<=0 for m in models):score=0.;flags.append('aperture sign mismatch')
            if not all(m['complete'] for m in models):flags.append('baseline coverage')
            if not all(m['support'] for m in models):flags.append('channel support')
            if variation>.3:flags.append('baseline sensitive')
        except ValueError:
            score=None;models=[];variation=None;flags.append('insufficient baseline/noise')
        if q['lo']<2 or q['hi']>=len(cube.frequency)-2:flags.append('spectral edge')
        if np.any(extraction['noise'][q['lo']:q['hi']+1]>3*np.median(extraction['noise'])):flags.append('high channel noise')
        q.update(snr=score,models=models,variation=variation,flags=flags,
                 significant=score is not None and score>=threshold and not flags,matches=[])
        q['aperture_measurement']=measure_feature(cube,extraction,q)
        features.append(q)
    return sorted(features,key=lambda q:q['frequency_GHz'])


def match_features(features,rows,z,tolerance):
    """Match measured centroids, never overlapping detection search windows."""
    from .measurements import matching_frequency
    matches={}
    for feature in features:
        found=[]
        for row in rows:
            expected=row['rest_GHz']/(1+z)
            offset=C*(expected-matching_frequency(feature))/expected
            if abs(offset)<=tolerance:
                found.append(dict(row,expected_GHz=expected,velocity_offset_kms=float(offset)))
        matches[feature['id']]=sorted(found,key=lambda r:abs(r['velocity_offset_kms']))
    return matches
