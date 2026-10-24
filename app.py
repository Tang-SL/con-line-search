"""CON Line Search: a guided, step-by-step line search. Run with streamlit run app.py.

Steps: 1 Load data → 2 Find lines → 3 Choose region → 4 Identify lines → 5 Make maps.
Only the current step is open; finished steps collapse to a one-line summary and
can be reopened. Settings appear in the step that uses them.
"""
from pathlib import Path
import hashlib, io, json, zipfile, os, sys, subprocess
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from conlines.espada_backend import run_espada, aperture_features, match_features
from conlines.plots import chart_image, spectrum_plot
from conlines.solutions import cube_fingerprint, feature_fingerprint, digest
from conlines.review_ui import identification_review
from conlines.measurements import matching_frequency
from conlines.evidence import add_evidence
from conlines import admit_lineid, amase_bridge
from conlines.reliability import channel_noise,component_properties,detection_reliability,neighbour_counts,identification_probabilities,chance_match_probability,anchor_statistics
from astropy.wcs import WCS
from conlines.core import read_cube,extract,query_catalogue,select_priors,moment,fits_bytes,C

ROOT=Path(__file__).resolve().parent
IC860_Z=3880/C
st.set_page_config(page_title='CON Line Search',page_icon='🔭',layout='wide',initial_sidebar_state='collapsed')
st.markdown('''<style>.block-container{max-width:1300px;padding-top:2rem}h1{letter-spacing:-.03em}
.stepper{display:flex;gap:6px;margin:.2rem 0 1.4rem 0;flex-wrap:wrap}
.stepper span{flex:1;min-width:120px;padding:8px 10px;border-radius:8px;font-size:.92rem;text-align:center;background:#eef2f6;color:#7a8794}
.stepper span.done{background:#e3f1ec;color:#0b6b4f}
.stepper span.now{background:#138c87;color:#fff;font-weight:600}
.todo{background:#f3f7fa;border-left:4px solid #138c87;padding:.6rem .9rem;border-radius:6px;margin:.2rem 0 .8rem 0}</style>''',unsafe_allow_html=True)


# ---------------------------------------------------------------- caches
@st.cache_resource(max_entries=8)
def load_path(path,mtime):return read_cube(path)
@st.cache_resource(max_entries=8)
def load_upload(digest,_file):
    _file.seek(0)
    return read_cube(_file)
@st.cache_data(max_entries=8,show_spinner=False)
def search_cube(key,widths,mode,clip,continuum_subtracted,_cube):
    return run_espada(_cube,widths,mode,clip,continuum_subtracted)
@st.cache_data(max_entries=16,show_spinner=False)
def identify_cube(key,_cube):return cube_fingerprint(_cube)
@st.cache_data(max_entries=32,show_spinner=False)
def measure(key,config,_cube):return extract(_cube,**config)
@st.cache_data(max_entries=64,show_spinner=False)
def screen_features(key,spw,_cube,_extraction,_run,threshold):
    # key covers the cubes, aperture and detection settings; aperture-only work, never a new ESPADA run.
    return aperture_features(_cube,_extraction,_run,threshold)
@st.cache_data(max_entries=8,show_spinner=False)
def preview_image(key,_cube):return np.nanmedian(_cube.data,axis=0)
@st.cache_data(max_entries=8,show_spinner=False)
def cube_noise(key,_cube):return channel_noise(_cube)
@st.cache_data(max_entries=16,show_spinner=False)
def strengths(key,_cube,_components,_noise):return component_properties(_cube,_components,_noise)
@st.cache_data(max_entries=16,show_spinner=False)
def region_areas(key,_cube,cx,cy,region_radius,near_factor=3.):
    yy,xx=np.indices(_cube.data.shape[1:]);sx,sy=_cube.pixel_arcsec
    r=np.hypot((xx-cx)*sx,(yy-cy)*sy)
    return dict(inside=int((r<=region_radius).sum()),near=int(((r>region_radius)&(r<=near_factor*region_radius)).sum()),far=int((r>near_factor*region_radius).sum()))


# ---------------------------------------------------------------- page helpers
STEPS=['Load data','Find lines','Choose region','Identify lines','Make maps']
st.session_state.setdefault('stage',1)

def go_to(n):st.session_state['stage']=n
def clamp(n):
    """An earlier step became incomplete (e.g. new cubes): make it the current step."""
    if st.session_state['stage']>n:st.session_state['stage']=n

def keep(key):st.session_state['_keep_'+key]=st.session_state[key]
def setting(widget,label,key,default,**kw):
    """A widget whose value survives while its step is hidden (Streamlit forgets hidden widgets)."""
    store='_keep_'+key
    value=st.session_state.get(store,default)
    # New data can make a remembered value invalid (an SPW no longer loaded, a pixel outside a smaller cube).
    if kw.get('options') is not None:
        options=list(kw['options'])
        if isinstance(value,list):value=[v for v in value if v in options]
        elif value not in options:value=default if default in options else options[0]
    if isinstance(value,(int,float)) and not isinstance(value,bool):
        if kw.get('min_value') is not None and value<kw['min_value']:value=kw['min_value']
        if kw.get('max_value') is not None and value>kw['max_value']:value=kw['max_value']
    st.session_state[store]=value
    if st.session_state.get(key,None)!=value:st.session_state[key]=value
    return widget(label,key=key,on_change=keep,args=(key,),**kw)
def remembered(key,default):return st.session_state.get('_keep_'+key,default)

st.title('CON Line Search')
st.caption('Find, identify and map spectral lines in ALMA cubes, one step at a time.')
top=st.columns([6,1])
stepper=top[0].empty()
if top[1].button('Start over',help='Clear all steps and start again with new data.'):
    st.session_state.clear();st.rerun()

def draw_stepper(current):
    parts=[]
    for i,name in enumerate(STEPS,1):
        cls='done' if i<current else 'now' if i==current else ''
        mark='✓ ' if i<current else ''
        parts.append(f'<span class="{cls}">{mark}{i} · {name}</span>')
    stepper.markdown('<div class="stepper">'+''.join(parts)+'</div>',unsafe_allow_html=True)

def finish(current):
    draw_stepper(current);st.stop()

def step(n):
    """The current step is open with a heading; a finished step folds into one line with its summary."""
    done=st.session_state['stage']>n
    if done:
        summary=st.session_state.get(f'summary_{n}','')
        return None,st.expander(f"{n} · {STEPS[n-1]}"+(f"  —  {summary}" if summary else ''),key=f'step_{n}',icon=':material/check_circle:'),done
    return st.empty(),st.container(),done

def header(head,n,done,summary='',todo=''):
    previous=st.session_state.get(f'summary_{n}')
    st.session_state[f'summary_{n}']=summary
    if done:
        # The folded label is drawn before the step runs; refresh it once if the step's result changed.
        if previous!=summary:st.rerun()
    else:head.markdown(f"### {n} · {STEPS[n-1]}"+(f"\n<div class='todo'><b>What to do:</b> {todo}</div>" if todo else ''),unsafe_allow_html=True)


# ================================================================ 1 · Load data
head,body,done=step(1)
loaded={};errors=[]
with body:
    input_mode=st.radio('Input source',['Upload FITS cubes','Local file paths','IC 860 full-field example'],horizontal=True,key='input_mode',
                        help='Give one FITS cube per spectral window. For large cubes, file paths avoid copying the data through the browser.')
    if input_mode=='Upload FITS cubes':
        uploads=st.file_uploader('Choose one or more cubes',type=['fits','fit'],accept_multiple_files=True,key='uploads',help='Files stay on this computer. Up to 4 GB per file.')
        for i,file in enumerate(uploads or []):
            try:
                key=hashlib.sha256(file.getbuffer()).hexdigest()
                loaded[key]=dict(name=f'{i+1}: {file.name}',cube=load_upload(key,file))
            except Exception as e:errors.append(f'{file.name}: {e}')
    elif input_mode=='Local file paths':
        paths=st.text_area('Absolute FITS paths, one per line',key='paths',placeholder='/Users/me/data/spw0.fits').splitlines()
        for path in filter(None,(p.strip() for p in paths)):
            try:
                p=Path(path).expanduser().resolve();key=hashlib.sha256((str(p)+str(p.stat().st_mtime_ns)).encode()).hexdigest();loaded[key]=dict(name=p.name,cube=load_path(str(p),p.stat().st_mtime_ns))
            except Exception as e:errors.append(f'{path}: {e}')
    else:
        base=ROOT/'targets/ic860/alma/fullfield/inputs'
        chosen=st.multiselect('Example spectral windows',[0,1,2,3],default=[0],key='example_spws',format_func=lambda s:f'SPW{s}')
        for spw in chosen:
            p=base/f'spw{spw}_full_beam060.fits'
            try:loaded[f'example{spw}']=dict(name=f'IC 860 · SPW{spw}',cube=load_path(str(p),p.stat().st_mtime_ns),cache_key=str(p)+str(p.stat().st_mtime_ns))
            except Exception as e:errors.append(f'SPW{spw}: Full-field example unavailable. Run scripts/ic860_fullfield_prepare.py on the original Sabine/b5 directory; see README. ({e})')
    for error in errors:st.error(error)
    if loaded:
        for k,v in loaded.items():v['fingerprint']=identify_cube(v.get('cache_key',k),v['cube'])
        st.dataframe(pd.DataFrame([dict(Cube=v['name'],Channels=len(v['cube'].frequency),Frequencies_GHz=f"{v['cube'].frequency.min():.3f}–{v['cube'].frequency.max():.3f}",Beam_arcsec=round(v['cube'].beam_arcsec,4),Size=' × '.join(map(str,v['cube'].data.shape[1:])),Frame=v['cube'].frame) for v in loaded.values()]),hide_index=True,width='stretch')
        beams=[v['cube'].beam_arcsec for v in loaded.values()]
        if not np.allclose(beams,beams[0],rtol=.01):st.warning('These cubes have different beam sizes, so their spectra are not resolution-matched.')
        frames={v['cube'].frame for v in loaded.values()}
        if len(frames)>1 or 'UNKNOWN' in frames:st.warning('Spectral reference frames differ or are unknown. No frame conversion is done; check frames before comparing velocities.')
        if not done:st.button('Continue →',type='primary',on_click=go_to,args=(2,),key='continue_1')
header(head,1,done,f"{len(loaded)} cube{'s' if len(loaded)!=1 else ''}: "+', '.join(v['name'] for v in loaded.values()),
       'choose the spectral cubes to analyse, one FITS cube per spectral window, then press <b>Continue</b>.')
if not loaded:
    clamp(1);st.info('Add a cube to begin. You can also open the IC 860 example without uploading files.');finish(1)
if st.session_state['stage']<2:finish(1)


# ================================================================ 2 · Find lines
head,body,done=step(2)
with body:
    with st.expander('Detection settings (advanced)'):
        a,b=st.columns(2)
        clip=setting(a.number_input,'ESPADA voxel threshold (σ)','clip',5.,min_value=3.,max_value=10.,step=.5,help='How far above the noise a voxel must be to seed a detection.')
        widths=setting(b.multiselect,'Trial widths (channels)','widths',[3,5,9],options=[3,5,7,9,13,17],help='Spectral smoothing widths tried by the finder.')
        mode=setting(a.selectbox,'Search for','mode','Both',options=['Both','Emission','Absorption'])
        continuum_subtracted=setting(b.checkbox,'Input cubes are continuum-subtracted','contsub',True,help='Untick to let SoFiA remove a linear continuum.')
    finder_settings=dict(engine='espada-v2-fullfield',widths=sorted(widths),mode=mode,clip=clip,continuum_subtracted=continuum_subtracted)
    detection_signature=digest(dict(datasets=[loaded[k]['fingerprint'] for k in loaded],bindings=list(loaded),settings=finder_settings))
    def component_prefix(k,i):return 'SPW'+k.removeprefix('example') if k.startswith('example') else f'Cube{i+1}'
    detection=st.session_state.get('detection')
    valid=bool(detection) and detection['signature']==detection_signature
    if st.button('Find lines in all cubes' if not valid else 'Find lines again',type='primary' if not valid else 'secondary',disabled=not widths,key='find_lines'):
        runs={};status=[]
        for i,(k,v) in enumerate(loaded.items()):
            try:
                with st.spinner(f"Searching {v['name']}… (can take a few minutes)"):
                    run=search_cube(v['fingerprint'],tuple(widths),mode,clip,continuum_subtracted,v['cube'])
                # IDs belong to the cube's detections, so they stay fixed when the aperture moves.
                ordered=sorted(run['components'],key=lambda c:(c['frequency_GHz'],c['sign']))
                runs[k]=dict(run,components=[dict(c,id=f'{component_prefix(k,i)}-F{n+1:03}') for n,c in enumerate(ordered)])
                status.append(dict(Cube=v['name'],Features_found=len(ordered),Status='Done'))
            except Exception as e:status.append(dict(Cube=v['name'],Features_found=None,Status=str(e)))
        st.session_state['detection']=dict(signature=detection_signature,runs=runs,status=status)
        st.session_state['stage']=max(st.session_state['stage'],3);st.rerun()
    if valid:
        st.dataframe(pd.DataFrame(detection['status']),hide_index=True,width='stretch')
        if not done:st.button('Continue →',type='primary',on_click=go_to,args=(3,),key='continue_2')
        if not done:
            for n,(k,run) in enumerate(detection['runs'].items()):
                components=run['components']
                with st.expander(f"Where were features found in {loaded[k]['name']}? ({len(components)})",expanded=n==0):
                    if not components:st.info('No features found in this cube.');continue
                    fig=chart_image(preview_image(loaded[k]['fingerprint'],loaded[k]['cube']),'Features on the median-channel image')
                    for sign,color,symbol in [('Absorption','#00c2ff','circle-open'),('Emission','#ff4d6d','diamond-open')]:
                        points=[c for c in components if c['sign']==sign]
                        if points:fig.add_trace(go.Scatter(x=[c['x'] for c in points],y=[c['y'] for c in points],mode='markers+text',text=[c['id'].rsplit('-',1)[-1] for c in points],textposition='top center',marker=dict(size=14,color=color,symbol=symbol,line=dict(width=2.5)),textfont=dict(color=color,size=11),name=sign,hovertemplate='%{text}<br>x %{x:.1f}, y %{y:.1f}<extra>'+sign+'</extra>'))
                    fig.update_layout(showlegend=True,legend=dict(orientation='h',y=-.15,x=0))
                    col1,col2=st.columns([1,1.2])
                    col1.plotly_chart(fig,width='stretch',key='detection_map_'+k)
                    col2.dataframe(pd.DataFrame([dict(ID=c['id'],Type=c['sign'],GHz=round(c['frequency_GHz'],6),x=round(c['x'],1),y=round(c['y'],1),Channels=f"{c['lo']}–{c['hi']}") for c in components]),hide_index=True,width='stretch',height=min(430,38+35*len(components)))
    elif detection:st.info('The cubes or detection settings changed since the last search. Press **Find lines in all cubes** again.')
total=sum(len(r['components']) for r in detection['runs'].values()) if valid else 0
header(head,2,done,f"{total} feature{'s' if total!=1 else ''} found across {len(loaded)} cube{'s' if len(loaded)!=1 else ''}",
       'press <b>Find lines in all cubes</b>. Every cube is searched for emission and absorption; this can take a few minutes per cube and only needs repeating if you change the cubes or the settings.')
if not valid:
    clamp(2);finish(2)
if st.session_state['stage']<3:finish(2)
all_components=[(k,c) for k,run in detection['runs'].items() for c in run['components']]


# ================================================================ 3 · Choose region
head,body,done=step(3)
with body:
    keys=list(loaded)
    if len(keys)>1:
        reference=setting(st.selectbox,'Reference SPW','reference',keys[0],options=keys,format_func=lambda k:loaded[k]['name'],help='The cube whose image is used to place the aperture.')
        if reference not in keys:reference=keys[0]
    else:reference=keys[0]
    ref=loaded[reference]['cube']
    centre_options=(['Continuum peak (example)'] if 'NUCX' in ref.header else [])+['Image centre','Brightness peak']+(['Detected component'] if all_components else [])+['Custom pixel position']
    center_mode=setting(st.radio,'Centre','centre',centre_options[0],options=centre_options,horizontal=True,help='Where to extract the spectrum, usually the nucleus.')
    preview=preview_image(loaded[reference]['fingerprint'],ref)
    if center_mode=='Continuum peak (example)':y0,x0=ref.header['NUCY'],ref.header['NUCX']
    elif center_mode=='Brightness peak':y0,x0=np.unravel_index(np.nanargmax(abs(preview)),preview.shape)
    elif center_mode=='Detected component':
        pick=setting(st.selectbox,'Centre on component','pick',0,options=range(len(all_components)),format_func=lambda i:f"{all_components[i][1]['id']} · {all_components[i][1]['sign']} · {all_components[i][1]['frequency_GHz']:.6f} GHz")
        pick=min(pick,len(all_components)-1)
        source_key,component=all_components[pick]
        world=WCS(loaded[source_key]['cube'].header).celestial.pixel_to_world(component['x'],component['y'])
        x0,y0=(float(v) for v in WCS(ref.header).celestial.world_to_pixel(world))
    else:y0,x0=(np.array(ref.data.shape[1:])-1)/2
    if center_mode=='Custom pixel position':
        a,b=st.columns(2)
        x0=setting(a.number_input,'Centre x (zero-based)','cx',float(x0),min_value=0.,max_value=float(ref.data.shape[2]-1))
        y0=setting(b.number_input,'Centre y (zero-based)','cy',float(y0),min_value=0.,max_value=float(ref.data.shape[1]-1))
    a,b,c,d=st.columns(4)
    shape=setting(a.selectbox,'Aperture','shape','Circle',options=['Circle','Ellipse','Rectangle','Single pixel'])
    radius=setting(b.number_input,'Radius / half-width (arcsec)','radius',float(ref.beam_arcsec),min_value=.0001,format='%.4f',help='Default: one beam.')
    if shape in ['Ellipse','Rectangle']:
        ratio=setting(c.number_input,'Minor / major ratio','ratio',1.,min_value=.1,max_value=1.,step=.1)
        angle=setting(d.number_input,'Angle from image x axis (deg)','angle',0.,min_value=-180.,max_value=180.,step=5.)
    else:ratio=1.;angle=0.
    with st.expander('More options: spectral windows, preview, significance and reliability (advanced)'):
        a,b=st.columns(2)
        linked=setting(a.checkbox,'Use the same sky position and angular aperture for every SPW','linked',True)
        aperture_zoom=setting(b.slider,'Aperture preview zoom','zoom',3.,min_value=1.,max_value=8.,step=.5,help='Display only; analysis always uses the full field.')
        threshold=setting(a.number_input,'Minimum integrated S/N','threshold',8.,min_value=3.,max_value=50.,step=1.,help='A feature is “significant” if its integrated S/N in the aperture reaches this under two baselines, with ≥2 adjacent channels at 3σ and ≤30% baseline sensitivity.')
        region_radius=setting(b.number_input,'Source region radius (arcsec)','region_radius',float(5*ref.beam_arcsec),min_value=float(ref.beam_arcsec),format='%.3f',help='Features outside this circle around the aperture centre are used to estimate how often artefacts occur. Default: five beams.')
        min_reliability=setting(a.slider,'Minimum detection reliability','min_reliability',.5,min_value=0.,max_value=1.,step=.05,help='Features below this are set aside as likely artefacts.')
        st.caption('Reliability compares each feature with similar features outside the source region, where real lines are rare, so it estimates how likely a feature is an artefact. It uses no line catalogue and is a cautious estimate, not a false-alarm probability.')
    sky=WCS(ref.header).celestial.pixel_to_world(x0,y0)
    extractions={};configs={}
    for k,v in loaded.items():
        cube=v['cube']
        x,y=WCS(cube.header).celestial.world_to_pixel(sky) if linked else ((cube.data.shape[2]-1)/2,(cube.data.shape[1]-1)/2)
        if not linked:
            with st.expander(f"Centre for {v['name']}"):
                x=st.number_input('x pixel',0.,float(cube.data.shape[2]-1),float(x),key=k+'x');y=st.number_input('y pixel',0.,float(cube.data.shape[1]-1),float(y),key=k+'y')
        cfg=dict(x=float(x),y=float(y),kind=shape,radius=radius,ratio=ratio,angle=angle);configs[k]=cfg
        try:
            ex=measure(v['fingerprint'],cfg,cube);extractions[k]=ex
            if not done:
                with st.expander(f"Spectrum of {v['name']} in this aperture",expanded=k==reference):
                    for note in cube.notes:st.info(note)
                    col1,col2=st.columns([1,1.4]);col1.plotly_chart(chart_image(preview_image(v['fingerprint'],cube),'Aperture (orange) and noise-control positions (white)',ex['mask'],ex['control_centers'],center=(cfg['x'],cfg['y']),zoom=aperture_zoom),width='stretch',key=k+'image');col2.plotly_chart(spectrum_plot(cube,ex),width='stretch',key=k+'spectrum')
                    st.download_button('Download this spectrum (CSV)',pd.DataFrame(dict(frequency_GHz=cube.frequency,aperture_mean=ex['spectrum'],control_noise=ex['noise'])).to_csv(index=False),file_name=f'spectrum_{k[:12]}.csv',key=k+'csv')
        except ValueError as e:st.error(f"{v['name']}: {e}")
    detection_settings=dict(finder_settings,threshold=threshold)
    detection_settings={k_:detection_settings[k_] for k_ in ['engine','threshold','widths','mode','clip','continuum_subtracted']}
    detection_context=dict(datasets=sorted([dict(fingerprint=v['fingerprint'],aperture={name:round(value,8) if isinstance(value,float) else value for name,value in configs[k].items()}) for k,v in loaded.items()],key=lambda v:v['fingerprint']),settings=detection_settings)
    signature=digest(dict(context=detection_context,bindings=list(loaded),measurement_version=1))
    results={};status=[]
    for k,ex in extractions.items():
        run=detection['runs'].get(k)
        if run is None:continue
        found=screen_features(signature,k,loaded[k]['cube'],ex,run,threshold)
        results[k]=found
        status.append(dict(SPW=loaded[k]['name'],Components=len(run['components']),Aperture_features=len(found),Significant=sum(q['significant'] for q in found)))
    analysis=dict(signature=signature,detection_signature=detection_signature,results=results,status=status,
                  runs={k:{name:value for name,value in run.items() if name!='components'} for k,run in detection['runs'].items()})
    st.session_state['analysis']=analysis
    # Detection reliability for every component (catalogue-free).
    scores={};positions=[];position_ids=[]
    for k,run in detection['runs'].items():
        cube=loaded[k]['cube'];components=run['components']
        if not components:continue
        wcs=WCS(cube.header).celestial
        props=strengths(detection_signature+k,cube,components,cube_noise(loaded[k]['fingerprint'],cube))
        world=wcs.pixel_to_world([c['x'] for c in components],[c['y'] for c in components])
        dra,ddec=sky.spherical_offsets_to(world)
        offsets=np.c_[dra.to_value('arcsec'),ddec.to_value('arcsec')];distance=np.hypot(offsets[:,0],offsets[:,1])
        cx,cy=(float(v) for v in wcs.world_to_pixel(sky))
        areas=region_areas(loaded[k]['fingerprint'],cube,round(cx,3),round(cy,3),round(region_radius,6))
        for c,p,r,rel,xy in zip(components,props,distance,detection_reliability(props,distance,region_radius,areas),offsets):
            scores[c['id']]=dict(p,**rel,distance_arcsec=float(r),spw=k);positions.append(xy);position_ids.append(c['id'])
    for ident,n in zip(position_ids,neighbour_counts(positions,max(v['cube'].beam_arcsec for v in loaded.values()))):scores[ident]['neighbours']=n
    def reliable(ident):
        rel=scores.get(ident,{}).get('reliability')
        return rel is not None and rel>=min_reliability
    lookup={q['id']:(k,q) for k,features in results.items() for q in features}
    n_feat=len(lookup);n_sig=sum(q['significant'] for k,q in lookup.values());n_rel=sum(reliable(i) for i in lookup)
    st.markdown('**Features at this position**')
    if not lookup:st.info('No features at this position. Move the aperture, for example to the nucleus or a detected feature.')
    else:
        m=st.columns(4)
        m[0].metric('Features here',n_feat);m[1].metric('Significant',n_sig,help='Pass the integrated S/N checks in this aperture.')
        m[2].metric('Reliable',n_rel,help=f'Detection reliability ≥ {min_reliability:.2f}: unlikely to be artefacts.');m[3].metric('Set aside',n_feat-n_rel,help='Likely artefacts; they are not identified unless you include them in Step 4.')
        def status_label(q):
            if not reliable(q['id']):return 'Likely artefact'
            return 'Reliable, significant' if q['significant'] else 'Reliable, weak'
        simple=[dict(ID=i,Window=loaded[k]['name'],Type=q['sign'],GHz=round(matching_frequency(q),6),SN=round(q['snr'],1) if q['snr'] is not None else None,Reliability=round(scores[i]['reliability'],2) if scores.get(i,{}).get('reliability') is not None else None,Status=status_label(q)) for i,(k,q) in lookup.items()]
        st.dataframe(pd.DataFrame(simple),hide_index=True,width='stretch',column_config={'SN':st.column_config.NumberColumn('S/N'),'GHz':st.column_config.NumberColumn('Observed GHz',format='%.6f')})
        with st.expander('All measurements and reliability details (advanced)'):
            full=[]
            for i,(k,q) in lookup.items():
                m_=q.get('aperture_measurement',{});s_=scores.get(i,{})
                full.append(dict(ID=i,Window=loaded[k]['name'],Sign=q['sign'],SNR=round(q['snr'],2) if q['snr'] is not None else None,Significant=q['significant'],Observed_GHz=matching_frequency(q),Centroid='Aperture' if m_.get('reliable') else 'ESPADA fallback',Centroid_error_kms=m_.get('centroid_error_kms') if m_.get('reliable') else None,Sigma_kms=m_.get('sigma_kms') if m_.get('reliable') else None,Channels=f"{q['lo']}–{q['hi']}",Reliability=s_.get('reliability'),Best_estimate=s_.get('reliability_best'),Peak_SNR=round(s_['peak_snr'],1) if s_ else None,Neighbours=s_.get('neighbours'),Stronger_inside=s_.get('observed'),Stronger_near=s_.get('control_near'),Stronger_far=s_.get('control_far'),Flags=', '.join(q['flags'])))
            st.dataframe(pd.DataFrame(full),hide_index=True,width='stretch')
            st.caption('Observed frequency uses the signed aperture centroid when measurable, otherwise the ESPADA centroid. Sigma is the second-moment width within the detected channels. Errors include spatial-control noise and baseline sensitivity, not redshift uncertainty.')
        if st.session_state.get('inspect') not in lookup:st.session_state.pop('inspect',None)
        inspect=st.selectbox('Inspect feature',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['sign']} · {status_label(lookup[ident][1])}",key='inspect')
        if not done:
            k,q=lookup[inspect];st.plotly_chart(spectrum_plot(loaded[k]['cube'],extractions[k],results[k],selected_id=inspect),width='stretch',key='inspect_'+inspect)
        if not done:st.button('Continue to line identification →',type='primary',on_click=go_to,args=(4,),key='continue_3')
header(head,3,done,f"{shape} aperture, radius {radius:.3f}″, at {center_mode.lower()}: {n_feat} feature{'s' if n_feat!=1 else ''}, {n_rel} reliable",
       'choose where to extract a spectrum (usually the nucleus). Features at that position are measured and checked; likely artefacts are set aside. Then press <b>Continue to line identification</b>.')
if not lookup:
    clamp(3);finish(3)
if st.session_state['stage']<4:finish(3)
feature_keys={ident:feature_fingerprint(loaded[k]['fingerprint'],q) for ident,(k,q) in lookup.items()}


# ================================================================ 4 · Identify lines
head,body,done=step(4)
methods={'Velocity likelihood (built-in)':'builtin','ADMIT LineID rules':'admit','Astro AMASE':'amase'}
explain={'builtin':'Looks up catalogue lines near each feature (Splatalogue) and ranks them by how well the velocity matches.',
         'admit':'The rule-based identification of the ALMA Data Mining Toolkit: a list of common bright lines first, then the simplest catalogue match (lowest energy, lighter molecule, strongest line).',
         'amase':'Models each candidate molecule’s full spectrum and checks its other lines are present, using chemical relationships between molecules (Fried et al.). Runs on the spectrum from Step 3 and takes several minutes.'}
with body:
    is_ic860=any(k.startswith('example') for k in loaded) or any('IC860' in str(v['cube'].header.get('OBJECT','')).replace(' ','').upper() for v in loaded.values())
    a,b=st.columns([1,2])
    z=setting(a.number_input,'Redshift z','z',IC860_Z if is_ic860 else 0.,min_value=-.99,format='%.7f',help='Source redshift, used to compare observed and rest frequencies.')
    b.caption(('IC 860 working value: cz = 3880 km/s (its reference frame is still to be verified). ' if is_ic860 else '')+'Used for all methods and for the map velocity axes.')
    if z==0:b.warning('Enter the source redshift first: with z = 0 the catalogue search looks for lines at rest frequencies.')
    method_name=setting(st.radio,'Identification method','id_method','Velocity likelihood (built-in)',options=list(methods),horizontal=True,help='All methods only suggest; you confirm each line below.')
    method=methods[method_name]
    st.caption(explain[method])
    with st.expander('Which features to identify (advanced)'):
        a,b,c=st.columns(3)
        include_provisional=setting(a.checkbox,'Also match provisional features','include_provisional',True,help='Include features that did not pass the S/N checks.')
        include_low=setting(b.checkbox,'Also match and review low-reliability features','include_low',False,help='Include features set aside as likely artefacts.')
        companion_tolerance=setting(c.number_input,'Companion velocity agreement (km/s)','companion_tolerance',50.,min_value=5.,max_value=500.,step=5.,help='For the evidence panel: other lines of the same molecule must sit at the same velocity within this tolerance.')
    def passes_reliability(ident):return include_low or reliable(ident)
    targets={k:[q for q in features if (q['significant'] or include_provisional) and passes_reliability(q['id'])] for k,features in results.items()}
    n_targets=sum(map(len,targets.values()))
    if not n_targets:st.warning('No features to identify with the current choices. Include provisional or low-reliability features (advanced), or change the region in Step 3.')

    def catalogue_controls(full):
        """Splatalogue settings; the built-in method uses all of them, ADMIT and AMASE only the filters."""
        a,b=st.columns([2,1])
        filters=setting(a.multiselect,'Astronomical filters','filters',['Extragalactic'],options=['Extragalactic','Hot Cores'],help='Which Splatalogue line lists to search.' + (' ADMIT used the full catalogue: choosing both is closest.' if not full else ''))
        if full:
            halfwidth=setting(b.number_input,'Match window ± (km/s)','halfwidth',200.,min_value=20.,max_value=2000.,step=20.,help='Catalogue lines within this velocity of a feature are candidates.')
            with st.expander('Catalogue filters (advanced)'):
                a,b=st.columns(2)
                energy=setting(a.number_input,'Maximum upper-state energy (K)','energy',300.,min_value=1.,max_value=10000.,step=50.)
                strength_filter=setting(b.checkbox,'Filter very weak transitions','strength_filter',True,help='Screens by Einstein A; not a flux prediction.')
                min_log_a=setting(b.number_input,'Minimum log₁₀ Einstein A (s⁻¹)','min_log_a',-5.,min_value=-15.,max_value=0.,step=.5) if strength_filter else None
                exclude=setting(a.checkbox,'Exclude deuterated / trans species','exclude',True)
                exclude_multiple=setting(a.checkbox,'Exclude multiple rare-isotope substitutions','exclude_multiple',True,help='Single substitutions (13C, 18O, 34S, 30Si…) are kept.')
                exclude_excited=setting(b.checkbox,'Exclude explicitly excited vibrational/torsional states','exclude_excited',False,help='Off by default: excited lines matter in CONs.')
        else:
            halfwidth=remembered('halfwidth',200.);energy=remembered('energy',300.);exclude=remembered('exclude',True);exclude_multiple=remembered('exclude_multiple',True);exclude_excited=remembered('exclude_excited',False)
            min_log_a=remembered('min_log_a',-5.) if remembered('strength_filter',True) else None
        return tuple('extragalactic' if v=='Extragalactic' else 'ism-hot-core' for v in filters),halfwidth,energy,exclude,min_log_a,exclude_multiple,exclude_excited

    def search_splatalogue(selected_filters,halfwidth,energy,exclude,min_log_a,exclude_multiple,exclude_excited,catalogue_settings):
        audits={};catalogues={};selected_catalogues={};query_status=[]
        for k,item in loaded.items():
            low=float(min(item['cube'].frequency));high=float(max(item['cube'].frequency))
            try:
                with st.spinner(f"Searching Splatalogue for {loaded[k]['name']}…"):
                    raw=query_catalogue(low*(1+z)/(1+halfwidth/C),high*(1+z)/(1-halfwidth/C),selected_filters)
                    selected,audit=select_priors(raw,low,high,z,halfwidth,energy,selected_filters,exclude,min_log_a,exclude_multiple,exclude_excited,return_audit=True)
                    audits[k]=audit;catalogues[k]=raw;selected_catalogues[k]=selected
                    query_status.append(dict(Cube=loaded[k]['name'],Status='Search completed',Transitions=len(selected)))
            except Exception as e:query_status.append(dict(Cube=loaded[k]['name'],Status=f'Search failed: {e}'))
        st.session_state['identifications']=dict(signature=digest(catalogue_settings),audits=audits,catalogues=catalogues,selected=selected_catalogues,status=query_status,queried_utc=datetime.now(timezone.utc).isoformat(),settings=catalogue_settings)

    if method=='amase':
        with st.expander('Optional: also search Splatalogue, for the companion-line evidence panel'):
            cat=catalogue_controls(False)
            st.caption('Not needed for Astro AMASE itself, which uses its own CDMS/JPL database.')
            amase_catalogue_button=st.button('Search Splatalogue',disabled=not cat[0],key='search_splatalogue_amase')
    else:cat=catalogue_controls(method=='builtin')
    selected_filters,halfwidth,energy,exclude,min_log_a,exclude_multiple,exclude_excited=cat
    catalogue_settings=dict(datasets=sorted(v['fingerprint'] for v in loaded.values()),z=z,tolerance=halfwidth,filters=selected_filters,energy=energy,exclude=exclude,min_log_a=min_log_a,exclude_multiple=exclude_multiple,exclude_excited=exclude_excited)
    match_settings=dict(catalogue_settings,detection=signature,include_provisional=include_provisional,companion_tolerance=companion_tolerance,linked_apertures=linked,method=method)
    identifications=st.session_state.get('identifications')
    if identifications and identifications['signature']!=digest(catalogue_settings):identifications=None
    matches={};method_suggestions=None;method_notes={};revision=identifications['queried_utc']+'|builtin' if identifications and method=='builtin' else None
    def channel_bounds(k,q):
        f=loaded[k]['cube'].frequency;half=abs(f[1]-f[0])/2
        return float(min(f[q['lo']],f[q['hi']])-half),float(max(f[q['lo']],f[q['hi']])+half)

    if method in ('builtin','admit'):
        if method=='admit':
            with st.expander('ADMIT settings (advanced)'):
                a,b,c=st.columns(3)
                admit_width=setting(a.number_input,'Tier 1 window (km/s; 0 = ADMIT default)','admit_width',0.,min_value=0.,max_value=1000.,step=10.,help='ADMIT default: 40 km/s, or 300 km/s when |v_sys| > 150 km/s.')
                admit_exotics=setting(b.checkbox,'Allow exotic atoms (Al, Cl, Mg, Mn, F, Li, Na, K, Ti)','admit_exotics',False)
                admit_recomb=setting(c.selectbox,'Recombination lines','admit_recomb','shallow',options=['shallow','off','deep'],help='shallow: H and He α/β only; deep: all.')
            match_settings.update(admit_tier1_width_kms=admit_width,admit_allow_exotics=admit_exotics,admit_recombination=admit_recomb)
        if st.button('Search Splatalogue' if not identifications else 'Search Splatalogue again',type='primary' if not identifications else 'secondary',disabled=not selected_filters or not n_targets,key='search_splatalogue'):
            search_splatalogue(*cat,catalogue_settings);st.rerun()
        if not identifications and st.session_state.get('identifications'):st.info('The redshift or catalogue settings changed. Press **Search Splatalogue** again.')
        if method=='builtin' and identifications:
            for k,features in targets.items():
                if k in identifications['selected']:matches.update(match_features(features,identifications['selected'][k],z,halfwidth))
        elif method=='admit' and identifications:
            method_suggestions={}
            for k,features in targets.items():
                if k not in identifications['catalogues'] or not features:continue
                f=loaded[k]['cube'].frequency
                feats=[dict(id=q['id'],spw=k,sign=q['sign'],centre_GHz=matching_frequency(q),lo_GHz=channel_bounds(k,q)[0],hi_GHz=channel_bounds(k,q)[1]) for q in features]
                for ident,r in admit_lineid.identify(feats,identifications['catalogues'][k],z,{k:(float(min(f)),float(max(f)))},None,admit_width or None,admit_exotics,admit_recomb).items():
                    rows=[dict(c,method_label=('ADMIT choice: '+r['rule']) if i==0 and r['winner'] is not None else ('ADMIT blend' if c in r['blends'] else 'ADMIT candidate')) for i,c in enumerate(r['candidates'])]
                    matches[ident]=rows;method_suggestions[ident]=rows[0] if r['winner'] is not None else None
                    method_notes[ident]=('ADMIT: '+'; '.join(r['notes'][:3])) if r['winner'] is not None else 'ADMIT: U (no catalogue line inside the feature)'
            revision=identifications['queried_utc']+'|admit|'+digest([admit_width,admit_exotics,admit_recomb])
    else:
        if amase_catalogue_button:
            search_splatalogue(*cat,catalogue_settings);st.rerun()
        if not st.session_state.get('amase_python'):st.session_state['amase_python']=amase_bridge.default_python()
        if not st.session_state.get('amase_data'):st.session_state['amase_data']=amase_bridge.default_data_dir()
        problems=amase_bridge.check_installation(st.session_state['amase_python'],st.session_state['amase_data'])
        if problems:
            st.warning('Astro AMASE is not installed yet. '+' '.join(problems))
            a,b=st.columns([1,2])
            if sys.platform=='darwin':
                if a.button('Set up Astro AMASE…',help='Opens Terminal and installs Astro AMASE, a private Python 3.11 and its 5 GB database into .local-tools/.'):
                    subprocess.Popen(['open','-a','Terminal',str(ROOT/'Set up Astro AMASE.command')])
                    st.info('Setup is running in Terminal. When it says “Done”, press “Use installed Astro AMASE”.')
            else:a.code(f'"{ROOT/"scripts"/"setup_amase.sh"}"',language='sh')
            def reload_paths():st.session_state.update(amase_python=amase_bridge.default_python(),amase_data=amase_bridge.default_data_dir())
            b.button('Use installed Astro AMASE',on_click=reload_paths,help='Reads the paths recorded by the setup.')
        a,b=st.columns(2)
        amase_temperature=setting(a.number_input,'Excitation temperature (K)','amase_temperature',150.,min_value=3.,max_value=2000.,step=10.,help='AMASE uses one temperature for all molecules; ~150 K suits hot nuclei.')
        amase_signs=setting(b.selectbox,'Analyse','amase_signs','Emission and absorption',options=['Emission and absorption','Emission only','Absorption only'],help='Absorption is analysed in a separate run on the sign-flipped spectrum.')
        with st.expander('Astro AMASE settings (advanced)'):
            a,b=st.columns(2)
            amase_python=a.text_input('Astro AMASE Python executable',key='amase_python')
            amase_data=b.text_input('Astro AMASE database folder',key='amase_data')
            a,b,c=st.columns(3)
            amase_refine=setting(a.checkbox,'Let AMASE refine v_sys and T_ex','amase_refine',False,help='Off: v_sys from the redshift, T_ex fixed. On: AMASE searches ±250 km/s and ±100 K.')
            amase_sigma=setting(b.number_input,'AMASE peak threshold (σ)','amase_sigma',5.,min_value=3.,max_value=20.,step=.5)
            amase_timeout=setting(c.number_input,'Time limit per run (min)','amase_timeout',120,min_value=1,max_value=600)
            a,b=st.columns([1.2,1])
            amase_atoms=setting(a.multiselect,'Allowed elements','amase_atoms',['C','H','O','N','S'],options=['C','H','O','N','S','Si','P','Cl','Na','Mg','Al','K','Fe'],help='AMASE default is C, H, O, N, S. Add Si for SiO and its isotopologues.')
            amase_size=setting(b.number_input,'Source size (arcsec; 0 = fills the beam)','amase_size',0.,min_value=0.,max_value=100.,step=.05)
            amase_seeds=setting(st.text_input,'Known molecules to seed chemical relevance (SMILES, comma-separated)','amase_seeds','',help="AMASE's hot-core suggestion: [C-]#[O+], CO, CCO, C=O, COC=O, CC#N, CC(C)=O")
        v_sys=amase_bridge.radio_velocity(z)
        ref_cube=loaded[reference]['cube']
        amase_params=dict(temperature=float(amase_temperature),sigma_threshold=float(amase_sigma),observation_type='interferometric',
                          beam_major_axis=float(ref_cube.header['BMAJ']*3600),beam_minor_axis=float(ref_cube.header['BMIN']*3600),
                          source_size=float(amase_size) if amase_size>0 else 1e20,valid_atoms=list(amase_atoms),overwrite_old_files=True)
        if amase_refine:amase_params.update(vlsr_range=[v_sys-250,v_sys+250],temperature_is_exact=False)
        else:amase_params.update(vlsr=v_sys)
        seeds=[x.strip() for x in amase_seeds.split(',') if x.strip()]
        if seeds:amase_params['known_molecules']=seeds
        signs={'Emission and absorption':['Emission','Absorption'],'Emission only':['Emission'],'Absorption only':['Absorption']}[amase_signs]
        amase_signature=amase_bridge.signature(signature,amase_params,signs,sorted(extractions))
        match_settings.update(amase_params=amase_params,amase_signs=signs)
        amase_run=st.session_state.get('amase')
        if amase_run and amase_run['signature']!=amase_signature:amase_run=None
        if st.button('Run Astro AMASE' if not amase_run else 'Run Astro AMASE again',type='primary' if not amase_run else 'secondary',disabled=not extractions,key='run_amase'):
            summaries={};errors_=[]
            amase_bridge.remember_paths(amase_python,amase_data)
            parts=[dict(frequency_GHz=loaded[k]['cube'].frequency,spectrum=ex['spectrum'],noise=ex['noise'],unit=loaded[k]['cube'].unit,bmaj_arcsec=loaded[k]['cube'].header['BMAJ']*3600,bmin_arcsec=loaded[k]['cube'].header['BMIN']*3600) for k,ex in extractions.items()]
            for sign in signs:
                try:
                    freq,intensity,rms=amase_bridge.build_spectrum(parts,1. if sign=='Emission' else -1.)
                    with st.spinner(f'Running Astro AMASE on the {sign.lower()} spectrum… this can take several minutes.'):
                        summaries[sign]=amase_bridge.run(amase_python,amase_data,freq,intensity,dict(amase_params,rms_noise=rms,save_name='conlines_'+sign.lower()+'_'+amase_signature[:8]),timeout=60*int(amase_timeout))
                except Exception as e:errors_.append(f'{sign}: {e}')
            st.session_state['amase']=dict(signature=amase_signature,summaries=summaries,errors=errors_,run_utc=datetime.now(timezone.utc).isoformat())
            st.rerun()
        if st.session_state.get('amase') and not amase_run:st.info('The region or AMASE settings changed. Press **Run Astro AMASE** again.')
        if amase_run:
            for e in amase_run['errors']:st.error(e)
            if amase_run['summaries']:
                st.caption('Astro AMASE run '+amase_run['run_utc'][:16].replace('T',' ')+' UTC · '+', '.join(f"{s}: v_sys {r['vlsr']:.1f} km/s, T_ex {r['temperature']:.0f} K, FWHM {r['linewidth']:.1f} km/s, {len(r['lines'])} lines" for s,r in amase_run['summaries'].items()))
                for s_,r in amase_run['summaries'].items():
                    files=r.get('output_files',{})
                    if files.get('interactive_plot') and Path(files['interactive_plot']).exists():
                        st.download_button(f'Download AMASE {s_.lower()} fit (interactive HTML)',Path(files['interactive_plot']).read_bytes(),f'amase_{s_.lower()}_fit.html','text/html',key='amase_html_'+s_)
            method_suggestions={}
            channel={k:float(abs(v['cube'].frequency[1]-v['cube'].frequency[0])) for k,v in loaded.items()}
            feats=[dict(id=q['id'],spw=k,sign=q['sign'],lo_GHz=channel_bounds(k,q)[0],hi_GHz=channel_bounds(k,q)[1]) for k,features in targets.items() for q in features]
            for sign,summary in amase_run['summaries'].items():
                for ident,r in amase_bridge.map_to_features(summary,feats,sign,channel).items():
                    matches[ident]=r['candidates'];method_suggestions[ident]=r['suggestion'];method_notes[ident]='AMASE: '+r['note']
            revision='amase|'+amase_run['run_utc']
    if matches:matches=add_evidence(matches,lookup,identifications['selected'] if identifications else {},loaded,z,companion_tolerance,linked,results)
    id_scores={}
    if method=='builtin':
        # Identification probabilities: which transition, given that the feature is real.
        store=st.session_state.get('assignments_by_feature',{})
        anchor_offsets=[]
        for ident,(k,q) in lookup.items():
            a_=store.get(feature_keys[ident])
            if a_ and a_.get('reviewed') and a_.get('transition') and q['significant'] and reliable(ident):
                anchor_offsets.append(C*(1-matching_frequency(q)/(a_['transition']['rest_GHz']/(1+z))))
        anchors=anchor_statistics(anchor_offsets)
        with st.expander('Scoring (advanced)'):
            def use_anchors():
                values=dict(id_sigma=round(max(anchors['sigma_kms'],5.),1),id_mean=round(anchors['mean_kms'],1))
                st.session_state.update(values);st.session_state.update({'_keep_'+k_:v_ for k_,v_ in values.items()})
            a,b,c=st.columns([1,1,1.4])
            id_sigma=setting(a.number_input,'Line velocity scatter σ (km/s)','id_sigma',30.,min_value=5.,max_value=500.,step=5.,help='Spread of true lines around the systemic velocity, before adding each centroid error.')
            id_mean=setting(b.number_input,'Mean velocity offset (km/s)','id_mean',0.,min_value=-500.,max_value=500.,step=5.,help='Corrects a systematic offset, for example from the redshift or reference frame.')
            if anchors:
                c.caption(f"{anchors['n']} reviewed, significant and reliable anchor lines: mean offset {anchors['mean_kms']:+.1f} km/s, scatter {anchors['sigma_kms']:.1f} km/s.")
                c.button('Use anchor calibration',on_click=use_anchors,disabled=anchors['n']<3,help='Needs at least three anchors.')
            else:c.caption('Once you have confirmed a few secure lines, they can calibrate σ and the mean offset.')
            st.caption(f'A candidate’s probability assumes the feature is real: Gaussian in velocity offset, against an “unidentified” alternative spread over ±{halfwidth:.0f} km/s. Combined = detection reliability × identification probability. “Chance match” is the probability that an unrelated catalogue line lies this close by chance.')
        match_settings.update(id_sigma_kms=id_sigma,id_mean_kms=id_mean)
        if identifications:
            for ident,options in matches.items():
                k,q=lookup[ident];measurement=q.get('aperture_measurement',{})
                error=measurement.get('centroid_error_kms') if measurement.get('reliable') else None
                probabilities,p_none=identification_probabilities([o['velocity_offset_kms'] for o in options],id_sigma,id_mean,halfwidth,error)
                rel=scores.get(ident,{}).get('reliability')
                for o,pr in zip(options,probabilities):o.update(id_probability=pr,combined_probability=rel*pr if rel is not None else None)
                f=loaded[k]['cube'].frequency;low,high=float(min(f)),float(max(f));middle=(low+high)/2
                coverage=(high-low)+2*halfwidth/C*middle
                id_scores[ident]=dict(p_unidentified=p_none,chance=chance_match_probability(len(identifications['selected'].get(k,[])),coverage,matching_frequency(q),id_sigma,id_mean))
    match_settings.update(region_radius_arcsec=region_radius,min_reliability=min_reliability,include_low_reliability=include_low)
    st.session_state['identification_matches']=matches
    if identifications:
        with st.expander('Catalogue search details'):
            st.caption('Splatalogue searched '+identifications['queried_utc'][:16].replace('T',' ')+' UTC.')
            st.dataframe(pd.DataFrame(identifications['status']),hide_index=True,width='stretch')
            st.json({loaded[k]['name']:a['reason_counts'] for k,a in identifications['audits'].items()})
            st.download_button('Download catalogue query and filter audit',json.dumps(identifications,indent=2),'live_catalogue_audit.json','application/json')
    solution_context=dict(detection=detection_context,redshift=z)
    solution_scope=digest(dict(context=solution_context,feature_keys=feature_keys))
    # Channel ranges for maps (edited in Step 5) travel with saved solutions.
    intervals_store=st.session_state.setdefault('map_intervals',{})
    restored=st.session_state.pop('solution_intervals_'+solution_scope,None)
    if restored:
        for ident,(lo,hi) in restored.items():intervals_store[feature_keys[ident]]=[int(lo),int(hi)]
    intervals={ident:tuple(intervals_store.get(feature_keys[ident],(q['lo'],q['hi']))) for ident,(k,q) in lookup.items()}
    ready=bool(matches) or (method=='builtin' and identifications) or method_suggestions is not None
    if ready:review_box=st.container();review_box.markdown('**Review the suggestions**')
    else:review_box=st.expander('Current identifications, and loading a saved solution',key='review_waiting')
    review_scores={ident:dict(reliability=scores.get(ident,{}).get('reliability'),region=scores.get(ident,{}).get('region'),method_note=method_notes.get(ident),**id_scores.get(ident,{})) for ident in lookup}
    provenance=dict(queried_utc=identifications['queried_utc'] if identifications else None,method=method,settings=match_settings) if (identifications or method_suggestions is not None) else None
    with review_box:assignments,solution=identification_review(lookup,matches,solution_context,feature_keys,loaded,extractions,intervals,provenance,revision,Path(os.environ.get('CONLINES_SOLUTION_DIR',ROOT/'solutions')),scores=review_scores,review_ids=[i for i in lookup if passes_reliability(i)],suggestions=method_suggestions,method=method)
    if not done:st.button('Continue to maps →',type='primary',on_click=go_to,args=(5,),key='continue_4')
reviewed=sum(a_['reviewed'] for a_ in assignments.values());named=sum(a_['transition'] is not None for a_ in assignments.values())
header(head,4,done,f"{method_name}: {named} feature{'s' if named!=1 else ''} with an identification, {reviewed} reviewed",
       {'builtin':'press <b>Search Splatalogue</b>, then review each suggestion with <b>Review identifications one by one</b>. Continue to maps when you are done.',
        'admit':'press <b>Search Splatalogue</b> (ADMIT then applies its rules), review each suggestion with <b>Review identifications one by one</b>, then continue to maps.',
        'amase':'press <b>Run Astro AMASE</b> (several minutes), review each suggestion with <b>Review identifications one by one</b>, then continue to maps.'}[method])
if st.session_state['stage']<5:finish(4)


# ================================================================ 5 · Make maps
head,body,done=step(5)
with body:
    label=lambda ident:f"{ident} · {lookup[ident][1]['sign']} · "+(assignments[ident]['transition']['species'] if assignments[ident]['transition'] else 'unidentified')
    default=[i for i in lookup if assignments[i]['reviewed'] and assignments[i]['transition']]
    adopted_ids=setting(st.multiselect,'Features to map','map_features',default,options=list(lookup),format_func=label,help='Defaults to the features you identified and reviewed.')
    adopted_ids=[i for i in adopted_ids if i in lookup]
    if not adopted_ids:st.info('Choose at least one feature to map.')
    else:
        st.caption('Each map integrates the feature over its channel range. Adjust a range if the line is wider or blended; this changes the maps only, not the detection.')
        table=pd.DataFrame([dict(ID=i,Identification=(assignments[i]['transition']['species'] if assignments[i]['transition'] else 'Unidentified'),Reviewed=assignments[i]['reviewed'],Start=intervals[i][0],End=intervals[i][1]) for i in adopted_ids])
        edited=st.data_editor(table,hide_index=True,width='stretch',key='intervals_'+digest(adopted_ids)+'_'+str(st.session_state.get('interval_revision_'+solution_scope,0)),disabled=['ID','Identification','Reviewed'],
                              column_config={'Start':st.column_config.NumberColumn('First channel',min_value=0,step=1),'End':st.column_config.NumberColumn('Last channel',min_value=0,step=1)})
        for _,r in edited.iterrows():
            st.session_state['map_intervals'][feature_keys[r.ID]]=[int(r.Start),int(r.End)];intervals[r.ID]=(int(r.Start),int(r.End))
        a,b=st.columns(2)
        baseline=setting(a.selectbox,'Map baseline subtraction','baseline','linear',options=['linear','constant','none'],help='Fitted to the channels beside each line, per pixel.')
        map_zoom=setting(b.slider,'Map preview zoom','map_zoom',3.,min_value=1.,max_value=8.,step=.5)
        if st.button('Generate maps & download bundle',type='primary'):
            output=io.BytesIO();manifest=dict(detection=analysis,identifications=dict(identifications or {},matches=matches,settings=match_settings,method=method,amase=st.session_state.get('amase') if method=='amase' else None),solution=solution,apertures=configs,map_preview_zoom=map_zoom,lines=[])
            with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
                for ident in adopted_ids:
                    k,q=lookup[ident];cube=loaded[k]['cube'];line=assignments[ident]['transition']
                    try:
                        lo,hi=intervals[ident]
                        im,h=moment(cube,lo,hi,line['rest_GHz'] if line else matching_frequency(q),z if line else 0.,baseline)
                        name=line['species'] if line else 'Unidentified '+ident
                        if line is None:
                            for keyword in ['RESTFRQ','RESTFREQ','REDSHIFT']:
                                if keyword in h:del h[keyword]
                            h['VREFGHZ']=(matching_frequency(q),'Observed reference for velocity channel widths')
                        h['IDSTAT']='REVIEWED' if assignments[ident]['reviewed'] else 'TENTATIVE'
                        h['LINEID']=name.encode('ascii','backslashreplace').decode();h['SRCFILE']=loaded[k]['name'].encode('ascii','backslashreplace').decode();h.add_history('Aperture: '+json.dumps(configs[k]))
                        archive.writestr(ident+'.fits',fits_bytes(im,h))
                        manifest['lines'].append(dict(id=ident,source=loaded[k]['name'],transition=line,lo=lo,hi=hi,baseline=baseline,frame=cube.frame,unit=h['BUNIT']))
                        fig=chart_image(im,f"{loaded[k]['name']} · {name} · channels {lo}–{hi}",center=(configs[k]['x'],configs[k]['y']),zoom=map_zoom)
                        fig.data[0].colorscale='RdBu_r';fig.data[0].zmid=0;fig.data[0].colorbar.title=h['BUNIT'];st.plotly_chart(fig,width='stretch',key='map_'+ident)
                        archive.writestr(ident+'.html',fig.to_html(include_plotlyjs=True))
                    except (ValueError,TypeError) as e:st.error(f'{ident}: {e}')
                features_csv=pd.DataFrame([dict(ID=i,Window=loaded[k]['name'],Sign=q['sign'],Observed_GHz=matching_frequency(q),SNR=q['snr'],Significant=q['significant'],Reliability=scores.get(i,{}).get('reliability'),Start=intervals[i][0],End=intervals[i][1],Identification=assignments[i]['transition']['species'] if assignments[i]['transition'] else None,Reviewed=assignments[i]['reviewed']) for i,(k,q) in lookup.items()])
                archive.writestr('identification_solution.json',json.dumps(solution,indent=2));archive.writestr('manifest.json',json.dumps(manifest,indent=2));archive.writestr('features.csv',features_csv.to_csv(index=False))
            if manifest['lines']:
                st.download_button('Download FITS maps + interactive previews + provenance',output.getvalue(),'adopted_line_maps.zip','application/zip',on_click='ignore')
                st.success(f"Generated {len(manifest['lines'])} signed maps. Absorption keeps its negative sign; nothing is clipped.")
header(head,5,done,'','choose the features to map, check their channel ranges, then press <b>Generate maps & download bundle</b>. Unidentified features can be mapped too.')
draw_stepper(5)
