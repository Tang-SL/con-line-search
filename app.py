"""Local, multi-SPW aperture spectroscopy workbench. Run with streamlit run app.py."""
from pathlib import Path
import hashlib, io, json, zipfile, os
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
st.set_page_config(page_title='CON Line Search',page_icon='🔭',layout='wide')
st.markdown('''<style>.block-container{max-width:1450px;padding-top:2rem}h1{letter-spacing:-.04em} [data-testid="stMetric"]{background:#eaf0f5;padding:16px;border-radius:12px}</style>''',unsafe_allow_html=True)
st.caption('CON LINE SEARCH  /  SPECTRAL WORKBENCH')
st.title('From a cube to an adopted line')
st.write('Detect emission and absorption in full cubes with ESPADA, screen the detections in comparable apertures across spectral windows, identify features with live Splatalogue, and export signed integrated maps.')

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
def identify_cube(key,_cube):
    fingerprint=cube_fingerprint(_cube)
    return fingerprint

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


with st.sidebar:
    st.subheader('Source & catalogue')
    preset=st.selectbox('Source preset',['IC 860','Custom source'])
    z=st.number_input('Redshift z',min_value=-.99,value=3880/C if preset=='IC 860' else 0.,format='%.7f',key='redshift_'+preset)
    st.caption('Use a redshift appropriate to the cube reference frame. The IC 860 preset is a working value; its precise frame remains to be verified.')
    filters=st.multiselect('Astronomical filters',['Extragalactic','Hot Cores'],default=['Extragalactic'])
    selected_filters=tuple('extragalactic' if v=='Extragalactic' else 'ism-hot-core' for v in filters)
    energy=st.number_input('Maximum upper-state energy (K)',1.,10000.,300.,step=50.)
    exclude=st.checkbox('Exclude deuterated / trans species',True)
    exclude_multiple=st.checkbox('Exclude multiple rare-isotope substitutions',True,help='Retains single substitutions such as 13C, 15N, 18O, 34S and 30Si. Applies to explicit isotope labels, not ordinary atom counts.')
    strength_filter=st.checkbox('Filter very weak transitions',True)
    min_log_a=st.number_input('Minimum log₁₀ Einstein A (s⁻¹)',-15.,0.,-5.,step=.5,disabled=not strength_filter) if strength_filter else None
    exclude_excited=st.checkbox('Exclude explicitly excited vibrational/torsional states',False,help='Off by default: excited lines can be important in hot cores and CONs. Mixed-state catalogue entries are retained unless the transition identifies an excited state.')
    st.caption('Default A ≥10⁻⁵ s⁻¹ is a screening choice, not a flux prediction. Weak transitions and multiple isotopic substitutions can still be real. Disable these filters for completeness checks.')
    st.caption('Catalogue source: live Splatalogue only. Each matching request fetches fresh results; no saved fallback.')
    halfwidth=st.number_input('Catalogue match ± velocity (km/s)',20.,2000.,200.,step=20.)
    companion_tolerance=st.number_input('Companion velocity agreement (km/s)',5.,500.,50.,step=5.,help='Other transitions must coincide with a detected feature at the same velocity offset within this tolerance. This is a compatibility check, not a probability.')
    st.divider();st.subheader('Detection settings')
    clip=st.number_input('ESPADA voxel threshold (σ)',3.,10.,5.,step=.5)
    continuum_subtracted=st.checkbox('Input cubes are continuum-subtracted',True)
    threshold=st.number_input('Minimum integrated S/N',3.,50.,8.,step=1.)
    widths=st.multiselect('Trial widths (channels)',[3,5,7,9,13,17],default=[3,5,9])
    mode=st.selectbox('Search for',['Both','Emission','Absorption'])
    st.caption('Two local baselines, ≥2 adjacent channels at 3σ, and ≤30% baseline sensitivity. S/N is not a search-wide false-alarm probability.')

st.subheader('1 · Load spectral windows')
input_mode=st.radio('Input source',['Upload FITS cubes','Local file paths','IC 860 full-field example'],horizontal=True)
loaded={};errors=[]
if input_mode=='Upload FITS cubes':
    uploads=st.file_uploader('Choose one or more cubes',type=['fits','fit'],accept_multiple_files=True,help='Local app: files stay on this computer. Up to 4 GB per file; available RAM still limits the batch.')
    for i,file in enumerate(uploads):
        try:
            key=hashlib.sha256(file.getbuffer()).hexdigest()
            loaded[key]=dict(name=f'{i+1}: {file.name}',cube=load_upload(key,file))
        except Exception as e:errors.append(f'{file.name}: {e}')
elif input_mode=='Local file paths':
    paths=st.text_area('Absolute FITS paths, one per line',help='Avoids browser upload copies for large cubes. This server is bound to localhost.').splitlines()
    for path in filter(None,(p.strip() for p in paths)):
        try:
            p=Path(path).expanduser().resolve();key=hashlib.sha256((str(p)+str(p.stat().st_mtime_ns)).encode()).hexdigest();loaded[key]=dict(name=p.name,cube=load_path(str(p),p.stat().st_mtime_ns))
        except Exception as e:errors.append(f'{path}: {e}')
else:
    base=ROOT/'targets/ic860/alma/fullfield/inputs'
    chosen=st.multiselect('Example SPWs',[0,1,2,3],default=[0])
    for spw in chosen:
        p=base/f'spw{spw}_full_beam060.fits'
        try:loaded[f'example{spw}']=dict(name=f'IC 860 · SPW{spw}',cube=load_path(str(p),p.stat().st_mtime_ns),cache_key=str(p)+str(p.stat().st_mtime_ns))
        except Exception as e:errors.append(f'SPW{spw}: Full-field example unavailable. Run scripts/ic860_fullfield_prepare.py on the original Sabine/b5 directory; see README. ({e})')
for error in errors:st.error(error)
if not loaded:
    st.info('Add a cube to begin. You can also open the full-field IC 860 example without uploading files.')
    st.stop()
summary=[]
for k,v in loaded.items():
    v['fingerprint']=identify_cube(v.get('cache_key',k),v['cube'])
    c=v['cube'];summary.append(dict(SPW=v['name'],Channels=len(c.frequency),Shape=str(c.data.shape[1:]),Minimum_GHz=c.frequency.min(),Maximum_GHz=c.frequency.max(),Frame=c.frame,Beam_arcsec=c.beam_arcsec))
st.dataframe(pd.DataFrame(summary),hide_index=True,width='stretch')
beams=[v['cube'].beam_arcsec for v in loaded.values()]
if not np.allclose(beams,beams[0],rtol=.01):st.warning('SPWs have different beam sizes. Their spectra are not resolution-matched across files.')
frames={v['cube'].frame for v in loaded.values()}
if len(frames)>1 or 'UNKNOWN' in frames:st.warning('Reference frames differ or are unknown. No frame conversion or spectral stacking is performed. Check frames before interpreting velocity matches.')

st.subheader('2 · Detect features with ESPADA')
st.caption('ESPADA/SoFiA searches each full 3D cube for emission and absorption. No aperture, redshift or line catalogue is involved, so detection runs once per cube and settings: choosing or moving an aperture later never reruns it.')
detection_settings=dict(engine='espada-v2-fullfield',threshold=threshold,widths=sorted(widths),mode=mode,clip=clip,continuum_subtracted=continuum_subtracted)
finder_settings={name:value for name,value in detection_settings.items() if name!='threshold'}
detection_signature=digest(dict(datasets=[loaded[k]['fingerprint'] for k in loaded],bindings=list(loaded),settings=finder_settings))
def component_prefix(k,i):return 'SPW'+k.removeprefix('example') if k.startswith('example') else f'Cube{i+1}'
if st.button('Run detection across all loaded SPWs',type='primary',disabled=not widths):
    runs={};status=[]
    for i,(k,v) in enumerate(loaded.items()):
        try:
            with st.spinner(f"Running ESPADA on {v['name']}…"):
                run=search_cube(v['fingerprint'],tuple(widths),mode,clip,continuum_subtracted,v['cube'])
            # IDs belong to the cube's detections, so they stay fixed when the aperture moves.
            ordered=sorted(run['components'],key=lambda c:(c['frequency_GHz'],c['sign']))
            runs[k]=dict(run,components=[dict(c,id=f'{component_prefix(k,i)}-F{n+1:03}') for n,c in enumerate(ordered)])
            status.append(dict(SPW=v['name'],Components=len(ordered),Status='Completed'))
        except Exception as e:status.append(dict(SPW=v['name'],Components=None,Status=str(e)))
    st.session_state['detection']=dict(signature=detection_signature,runs=runs,status=status)
detection=st.session_state.get('detection')
if not detection or detection['signature']!=detection_signature:
    st.info('Cubes or detection settings changed. Rerun detection to continue.' if detection else 'Run detection to continue. The aperture is chosen after detection.')
    st.stop()
st.dataframe(pd.DataFrame(detection['status']),hide_index=True,width='stretch')
for n,(k,run) in enumerate(detection['runs'].items()):
    components=run['components']
    with st.expander(f"{loaded[k]['name']} · {len(components)} ESPADA components in the full field",expanded=n==0):
        if not components:
            st.info('No components detected in this cube.');continue
        fig=chart_image(preview_image(loaded[k]['fingerprint'],loaded[k]['cube']),'Detected components on the median-channel image')
        for sign,color,symbol in [('Absorption','#00c2ff','circle-open'),('Emission','#ff4d6d','diamond-open')]:
            points=[c for c in components if c['sign']==sign]
            if points:fig.add_trace(go.Scatter(x=[c['x'] for c in points],y=[c['y'] for c in points],mode='markers+text',text=[c['id'].rsplit('-',1)[-1] for c in points],textposition='top center',marker=dict(size=14,color=color,symbol=symbol,line=dict(width=2.5)),textfont=dict(color=color,size=11),name=sign,hovertemplate='%{text}<br>x %{x:.1f}, y %{y:.1f}<extra>'+sign+'</extra>'))
        fig.update_layout(showlegend=True,legend=dict(orientation='h',y=-.15,x=0))
        col1,col2=st.columns([1,1.2])
        col1.plotly_chart(fig,width='stretch',key='detection_map_'+k)
        col2.dataframe(pd.DataFrame([dict(ID=c['id'],Sign=c['sign'],ESPADA_GHz=round(c['frequency_GHz'],6),x=round(c['x'],1),y=round(c['y'],1),Channels=f"{c['lo']}–{c['hi']}") for c in components]),hide_index=True,width='stretch',height=min(430,38+35*len(components)))
all_components=[(k,c) for k,run in detection['runs'].items() for c in run['components']]

st.subheader('3 · Choose an aperture and screen features')
keys=list(loaded);reference=st.selectbox('Reference SPW',keys,format_func=lambda k:loaded[k]['name']);ref=loaded[reference]['cube']
a,b,c,d=st.columns(4)
with a:shape=st.selectbox('Aperture',['Circle','Ellipse','Rectangle','Single pixel'])
with b:radius=st.number_input('Radius / half-width (arcsec)',min_value=.0001,value=float(ref.beam_arcsec),format='%.4f',key='radius_'+reference)
with c:ratio=st.number_input('Minor / major ratio',.1,1.,1.,step=.1,disabled=shape in ['Circle','Single pixel'])
with d:angle=st.number_input('Angle from image x axis (deg)',-180.,180.,0.,step=5.,disabled=shape in ['Circle','Single pixel'])
if shape in ['Circle','Single pixel']:ratio=1.;angle=0.
centre_options=(['Continuum peak (example)'] if 'NUCX' in ref.header else [])+['Image centre','Brightness peak']+(['Detected component'] if all_components else [])+['Custom pixel position']
center_mode=st.radio('Centre',centre_options,horizontal=True)
preview=preview_image(loaded[reference]['fingerprint'],ref)
if center_mode=='Continuum peak (example)':y0,x0=ref.header['NUCY'],ref.header['NUCX']
elif center_mode=='Brightness peak':
    # Useful for continuum cubes; inspect the preview for continuum-subtracted inputs.
    y0,x0=np.unravel_index(np.nanargmax(abs(preview)),preview.shape)
elif center_mode=='Detected component':
    pick=st.selectbox('Centre on component',range(len(all_components)),format_func=lambda i:f"{all_components[i][1]['id']} · {all_components[i][1]['sign']} · {all_components[i][1]['frequency_GHz']:.6f} GHz")
    source_key,component=all_components[pick]
    world=WCS(loaded[source_key]['cube'].header).celestial.pixel_to_world(component['x'],component['y'])
    x0,y0=(float(v) for v in WCS(ref.header).celestial.world_to_pixel(world))
else:y0,x0=(np.array(ref.data.shape[1:])-1)/2
if center_mode=='Custom pixel position':
    a,b=st.columns(2);x0=a.number_input('Centre x (zero-based)',0.,float(ref.data.shape[2]-1),float(x0));y0=b.number_input('Centre y (zero-based)',0.,float(ref.data.shape[1]-1),float(y0))
linked=st.checkbox('Use the same sky position and angular aperture for every SPW',True)
sky=WCS(ref.header).celestial.pixel_to_world(x0,y0)
st.caption(f'Reference centre: {sky.to_string("hmsdms")} · Spectra are aperture means, not total fluxes. White dots show spatial noise-control apertures; inspect them for source contamination.')
aperture_zoom=st.slider('Aperture preview zoom',1.,8.,3.,step=.5,help='Display only. Extraction, noise controls and detection use the full loaded cube.')
st.caption('The full spatial field is used for analysis. Set preview zoom to 1 to inspect all control positions; zoom never removes pixels.')
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
        with st.expander(f"{v['name']} · {ex['pixels']} aperture pixels · {len(ex['controls'])} controls",expanded=k==reference):
            for note in cube.notes:st.info(note)
            col1,col2=st.columns([1,1.4]);col1.plotly_chart(chart_image(preview_image(v['fingerprint'],cube),'Median-channel image',ex['mask'],ex['control_centers'],center=(cfg['x'],cfg['y']),zoom=aperture_zoom),width='stretch',key=k+'image');col2.plotly_chart(spectrum_plot(cube,ex),width='stretch',key=k+'spectrum')
            csv=pd.DataFrame(dict(frequency_GHz=cube.frequency,aperture_mean=ex['spectrum'],control_noise=ex['noise'])).to_csv(index=False)
            st.download_button('Download extracted spectrum',csv,file_name=f'spectrum_{k[:12]}.csv',key=k+'csv')
    except ValueError as e:st.error(f"{v['name']}: {e}")

st.markdown('**Features in this aperture**')
st.caption('ESPADA components whose 3D masks intersect your aperture are screened automatically using its extracted spectrum and spatial noise controls; moving the aperture updates this table without rerunning detection. Green features pass the integrated S/N and baseline checks; others are provisional and need review. Neither threshold is a search-wide false-alarm probability.')
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
st.dataframe(pd.DataFrame(status),hide_index=True,width='stretch')

st.markdown('**Detection reliability**')
st.caption('Is a component real? This uses no line catalogue. Components outside the source region are the artefact control sample: for each component inside it, the app counts same-sign components outside that are at least as strong in both peak integrated S/N and number of channels, scales that count to the source-region area, and compares it with the number observed inside. The near (≤3× radius) and far control zones are scored separately and the worse one is used. Reliability is 1 − (84% upper limit of expected artefacts) / (observed). It assumes real features outside the region are rare, so it errs on the cautious side; it is not a calibrated false-alarm probability.')
a,b,c=st.columns(3)
region_radius=a.number_input('Source region radius (arcsec)',min_value=float(ref.beam_arcsec),value=float(5*ref.beam_arcsec),format='%.3f',help='Centred on the aperture centre. Default: five beams.',key='region_radius_'+reference)
min_reliability=b.slider('Minimum detection reliability',0.,1.,.5,step=.05,help='Features below this are not matched or reviewed unless you include them.')
include_low=c.checkbox('Also match and review low-reliability features',False)
scores={};positions=[];position_ids=[]
for k,run in detection['runs'].items():
    cube=loaded[k]['cube'];components=run['components']
    if not components:continue
    wcs=WCS(cube.header).celestial
    props=strengths(detection_signature+k,cube,components,cube_noise(loaded[k]['fingerprint'],cube))
    world=wcs.pixel_to_world([c['x'] for c in components],[c['y'] for c in components])
    dra,ddec=sky.spherical_offsets_to(world)
    offsets=np.c_[dra.to_value('arcsec'),ddec.to_value('arcsec')];radius=np.hypot(offsets[:,0],offsets[:,1])
    cx,cy=(float(v) for v in wcs.world_to_pixel(sky))
    areas=region_areas(loaded[k]['fingerprint'],cube,round(cx,3),round(cy,3),round(region_radius,6))
    for c,p,r,rel,xy in zip(components,props,radius,detection_reliability(props,radius,region_radius,areas),offsets):
        scores[c['id']]=dict(p,**rel,distance_arcsec=float(r),spw=k);positions.append(xy);position_ids.append(c['id'])
for ident,n in zip(position_ids,neighbour_counts(positions,max(v['cube'].beam_arcsec for v in loaded.values()))):scores[ident]['neighbours']=n
def passes_reliability(ident):
    rel=scores.get(ident,{}).get('reliability')
    return include_low or (rel is not None and rel>=min_reliability)
in_region=[s for s in scores.values() if s['region']=='Source']
st.caption(f"{len(in_region)} of {len(scores)} components lie inside the {region_radius:.3f}″ source region; {sum(1 for s in in_region if s['reliability'] is not None and s['reliability']>=min_reliability)} reach reliability ≥ {min_reliability:.2f}. Neighbours count other components (any sign or window) within one beam.")
with st.expander('Reliability of every component in the source region'):
    if in_region:
        st.dataframe(pd.DataFrame([dict(ID=i,Sign=s['sign'],Reliability=s['reliability'],Best_estimate=s['reliability_best'],Peak_SNR=round(s['peak_snr'],1),Channels=s['channels'],Pixels=s['pixels'],Neighbours=s['neighbours'],Distance_arcsec=round(s['distance_arcsec'],3),Stronger_inside=s['observed'],Stronger_near=s['control_near'],Stronger_far=s['control_far']) for i,s in scores.items() if s['region']=='Source']).sort_values('Reliability',ascending=False),hide_index=True,width='stretch')
    else:st.info('No components inside the source region. Check the aperture centre and region radius.')
rows=[];lookup={}
for k,features in analysis['results'].items():
    for q in features:
        lookup[q['id']]=(k,q)
        measurement=q.get('aperture_measurement',{})
        rows.append(dict(ID=q['id'],SPW=loaded[k]['name'],Sign=q['sign'],SNR=round(q['snr'],2) if q['snr'] is not None else None,Significant=q['significant'],Observed_GHz=matching_frequency(q),Centroid='Aperture' if measurement.get('reliable') else 'ESPADA fallback',Centroid_error_kms=measurement.get('centroid_error_kms') if measurement.get('reliable') else None,Sigma_kms=measurement.get('sigma_kms') if measurement.get('reliable') else None,Start=q['lo'],End=q['hi'],Reliability=scores.get(q['id'],{}).get('reliability'),Region=scores.get(q['id'],{}).get('region'),Peak_SNR=round(scores[q['id']]['peak_snr'],1) if q['id'] in scores else None,Channels=q['hi']-q['lo']+1,Neighbours=scores.get(q['id'],{}).get('neighbours'),Flags=', '.join(q['flags'])))
if not rows:st.info('No ESPADA components intersect this aperture. Try another aperture or inspect detection settings.');st.stop()
if not any(r['Significant'] for r in rows):st.info('No features pass the current aperture significance screen. ESPADA components remain below for provisional review; enable “Also match provisional features” to inspect their possible identifications.')
solution_context=dict(detection=detection_context,redshift=z)
feature_keys={ident:feature_fingerprint(loaded[k]['fingerprint'],q) for ident,(k,q) in lookup.items()}
solution_scope=digest(dict(context=solution_context,feature_keys=feature_keys))
restored_intervals=st.session_state.get('solution_intervals_'+solution_scope,{})
for row in rows:
    if row['ID'] in restored_intervals:row['Start'],row['End']=restored_intervals[row['ID']]
st.write('Feature IDs describe detections, not molecule assignments. Integration channels are zero-based and inclusive. Editing them changes maps only; the reported significance belongs to the original ESPADA interval.')
st.caption('Observed frequency uses the signed aperture centroid when measurable; otherwise it explicitly falls back to ESPADA. Sigma is the second-moment width inside the original interval, not a Gaussian FWHM. Errors include spatial-control noise and baseline sensitivity; they exclude source-redshift uncertainty.')
edited=st.data_editor(pd.DataFrame(rows).style.apply(lambda row: ['background-color: #d9f0e7; color: #075e46' if row['Significant'] else '' for _ in row],axis=1),hide_index=True,width='stretch',key='review_'+signature+'_'+str(st.session_state.get('interval_revision_'+solution_scope,0)),disabled=[c for c in rows[0] if c not in ['Start','End']],column_config={'Start':st.column_config.NumberColumn(min_value=0,step=1),'End':st.column_config.NumberColumn(min_value=0,step=1)})
inspect=st.selectbox('Inspect feature',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['sign']}")
k,q=lookup[inspect];st.plotly_chart(spectrum_plot(loaded[k]['cube'],extractions[k],analysis['results'][k],selected_id=inspect),width='stretch',key='inspect_'+inspect)

st.subheader('4 · Identify features')
methods={'Velocity likelihood (built-in)':'builtin','ADMIT LineID rules':'admit','Astro AMASE':'amase'}
method_name=st.radio('Identification method',list(methods),horizontal=True,key='id_method',help='Built-in: closest transitions scored by velocity offset. ADMIT: the rule-based LineID of the ALMA Data Mining Toolkit. Astro AMASE: LTE line-pattern checks plus chemical relevance (Fried et al.).')
method=methods[method_name]
st.caption({'builtin':'Candidates within the match window, ranked by a velocity-offset probability against an “unidentified” alternative.',
            'admit':'ADMIT’s LineID rules, ported: its Tier 1 list of bright common lines first (±40 km/s, or ±300 km/s when the systemic velocity exceeds 150 km/s), then catalogue lines inside each feature narrowed to one by lowest upper-state energy, molecular mass, isotopologue and line-strength rules. Uses the live Splatalogue rows of the selected astronomical filters without the energy or Einstein A screens; Hot Cores plus Extragalactic is closest to ADMIT’s full catalogue.',
            'amase':'Astro AMASE runs on the aperture spectrum in its own Python environment, with its own CDMS/JPL catalogues: LTE simulation of each candidate, checks that its other strong lines are present and none is predicted too strong, and chemical relevance to molecules already assigned. Absorption features are analysed in a separate run on the sign-flipped spectrum.'}[method]+' All suggestions stay tentative until reviewed.')
include_provisional=st.checkbox('Also match provisional features',True)
targets={k:[q for q in features if (q['significant'] or include_provisional) and passes_reliability(q['id'])] for k,features in analysis['results'].items()}
# The live query depends only on the cubes' frequency coverage and catalogue settings.
# Matches are recomputed from that query for the current aperture features, so
# moving the aperture never needs a new catalogue request.
catalogue_settings=dict(datasets=sorted(v['fingerprint'] for v in loaded.values()),z=z,tolerance=halfwidth,filters=selected_filters,energy=energy,exclude=exclude,min_log_a=min_log_a,exclude_multiple=exclude_multiple,exclude_excited=exclude_excited)
catalogue_signature=digest(catalogue_settings)
match_settings=dict(catalogue_settings,detection=signature,include_provisional=include_provisional,companion_tolerance=companion_tolerance,linked_apertures=linked,method=method)
if method=='amase':st.caption('Splatalogue is optional with Astro AMASE: it only feeds the companion-line evidence panel.')
if st.button('Match features with live Splatalogue',disabled=not selected_filters):
    audits={};catalogues={};selected_catalogues={};query_status=[]
    for k,item in loaded.items():
        low=float(min(item['cube'].frequency));high=float(max(item['cube'].frequency))
        try:
            with st.spinner(f"Fetching live Splatalogue for {loaded[k]['name']}…"):
                raw=query_catalogue(low*(1+z)/(1+halfwidth/C),high*(1+z)/(1-halfwidth/C),selected_filters)
                selected,audit=select_priors(raw,low,high,z,halfwidth,energy,selected_filters,exclude,min_log_a,exclude_multiple,exclude_excited,return_audit=True)
                audits[k]=audit;catalogues[k]=raw;selected_catalogues[k]=selected
                query_status.append(dict(SPW=loaded[k]['name'],Status='Live query completed',Transitions=len(selected)))
        except Exception as e:query_status.append(dict(SPW=loaded[k]['name'],Status=f'Live query failed: {e}'))
    st.session_state['identifications']=dict(signature=catalogue_signature,audits=audits,catalogues=catalogues,selected=selected_catalogues,status=query_status,queried_utc=datetime.now(timezone.utc).isoformat(),settings=catalogue_settings)
identifications=st.session_state.get('identifications')
if identifications and identifications['signature']!=catalogue_signature:
    st.info('Identification settings changed. Fetch live matches again. Feature detections remain valid.');identifications=None
matches={};method_suggestions=None;method_notes={};revision=identifications['queried_utc'] if identifications else None
def channel_bounds(k,q):
    f=loaded[k]['cube'].frequency;half=abs(f[1]-f[0])/2
    return float(min(f[q['lo']],f[q['hi']])-half),float(max(f[q['lo']],f[q['hi']])+half)
if method=='builtin' and identifications:
    for k,features in targets.items():
        if k in identifications['selected']:matches.update(match_features(features,identifications['selected'][k],z,halfwidth))
elif method=='admit':
    with st.expander('ADMIT settings',expanded=False):
        a,b,c=st.columns(3)
        admit_width=a.number_input('Tier 1 window (km/s; 0 = ADMIT default)',0.,1000.,0.,step=10.,help='ADMIT default: 40 km/s, or 300 km/s when |v_sys| > 150 km/s.')
        admit_exotics=b.checkbox('Allow exotic atoms (Al, Cl, Mg, Mn, F, Li, Na, K, Ti)',False)
        admit_recomb=c.selectbox('Recombination lines',['shallow','off','deep'],help='shallow: H and He α/β only; deep: all.')
    match_settings.update(admit_tier1_width_kms=admit_width,admit_allow_exotics=admit_exotics,admit_recombination=admit_recomb)
    if identifications:
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
elif method=='amase':
    with st.expander('Astro AMASE settings',expanded=not st.session_state.get('amase')):
        st.caption('Install Astro AMASE in its own environment (see README: “Astro AMASE”) and download its database files. The paths can also come from the AMASE_PYTHON and AMASE_DATA environment variables.')
        a,b=st.columns(2)
        amase_python=a.text_input('Astro AMASE Python executable',value=amase_bridge.default_python(),key='amase_python')
        amase_data=b.text_input('Astro AMASE database folder',value=amase_bridge.default_data_dir(),key='amase_data')
        a,b,c,d=st.columns(4)
        amase_temperature=a.number_input('Excitation temperature (K)',3.,2000.,150.,step=10.,help='AMASE uses one temperature for all molecules.')
        amase_refine=b.checkbox('Let AMASE refine v_sys and T_ex',False,help='Off: v_sys from the redshift, T_ex fixed. On: AMASE searches ±250 km/s and ±100 K.')
        amase_sigma=c.number_input('AMASE peak threshold (σ)',3.,20.,5.,step=.5)
        amase_signs=d.selectbox('Analyse',['Emission and absorption','Emission only','Absorption only'])
        a,b,c=st.columns([1.2,1,1])
        amase_atoms=a.multiselect('Allowed elements',['C','H','O','N','S','Si','P','Cl','Na','Mg','Al','K','Fe'],default=['C','H','O','N','S'],help='AMASE default is C, H, O, N, S. Add Si for SiO and its isotopologues.')
        amase_size=b.number_input('Source size (arcsec; 0 = fills the beam)',0.,100.,0.,step=.05)
        amase_timeout=c.number_input('Time limit per run (min)',1,600,120)
        amase_seeds=st.text_input('Known molecules to seed chemical relevance (SMILES, comma-separated)',value='',help="AMASE's hot-core suggestion: [C-]#[O+], CO, CCO, C=O, COC=O, CC#N, CC(C)=O")
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
    st.caption(f'Systemic velocity passed to AMASE: {v_sys:.1f} km/s (radio convention, equivalent to z = {z:.7f}).')
    if st.button('Run Astro AMASE',type='primary',disabled=not extractions):
        summaries={};errors=[]
        amase_bridge.remember_paths(amase_python,amase_data)
        parts=[dict(frequency_GHz=loaded[k]['cube'].frequency,spectrum=ex['spectrum'],noise=ex['noise'],unit=loaded[k]['cube'].unit,bmaj_arcsec=loaded[k]['cube'].header['BMAJ']*3600,bmin_arcsec=loaded[k]['cube'].header['BMIN']*3600) for k,ex in extractions.items()]
        for sign in signs:
            try:
                freq,intensity,rms=amase_bridge.build_spectrum(parts,1. if sign=='Emission' else -1.)
                with st.spinner(f'Running Astro AMASE on the {sign.lower()} spectrum… this can take several minutes.'):
                    summaries[sign]=amase_bridge.run(amase_python,amase_data,freq,intensity,dict(amase_params,rms_noise=rms,save_name='conlines_'+sign.lower()+'_'+amase_signature[:8]),timeout=60*int(amase_timeout))
            except Exception as e:errors.append(f'{sign}: {e}')
        st.session_state['amase']=dict(signature=amase_signature,summaries=summaries,errors=errors,run_utc=datetime.now(timezone.utc).isoformat())
    amase_run=st.session_state.get('amase')
    if amase_run and amase_run['signature']!=amase_signature:
        st.info('Aperture, detection or AMASE settings changed. Run Astro AMASE again.');amase_run=None
    if amase_run:
        for e in amase_run['errors']:st.error(e)
        if amase_run['summaries']:
            st.caption('Astro AMASE run: '+amase_run['run_utc']+' · '+', '.join(f"{s}: v_sys {r['vlsr']:.1f} km/s, T_ex {r['temperature']:.0f} K, FWHM {r['linewidth']:.1f} km/s, {len(r['lines'])} lines" for s,r in amase_run['summaries'].items()))
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
if matches:matches=add_evidence(matches,lookup,identifications['selected'] if identifications else {},loaded,z,companion_tolerance,linked,analysis['results'])
id_scores={}
if method=='builtin':
    # Identification probabilities: which transition, given that the feature is real.
    store=st.session_state.get('assignments_by_feature',{})
    anchor_offsets=[]
    for ident,(k,q) in lookup.items():
        a=store.get(feature_keys[ident])
        if a and a.get('reviewed') and a.get('transition') and q['significant'] and (scores.get(ident,{}).get('reliability') or 0)>=min_reliability:
            anchor_offsets.append(C*(1-matching_frequency(q)/(a['transition']['rest_GHz']/(1+z))))
    anchors=anchor_statistics(anchor_offsets)
    st.session_state.setdefault('id_sigma',30.);st.session_state.setdefault('id_mean',0.)
    def use_anchors():st.session_state.update(id_sigma=round(max(anchors['sigma_kms'],5.),1),id_mean=round(anchors['mean_kms'],1))
    a,b,c=st.columns([1,1,1.4])
    id_sigma=a.number_input('Line velocity scatter σ (km/s)',5.,500.,step=5.,key='id_sigma',help='Spread of true lines around the systemic velocity, before adding each centroid error.')
    id_mean=b.number_input('Mean velocity offset (km/s)',-500.,500.,step=5.,key='id_mean',help='Corrects a systematic offset, for example from the redshift or reference frame.')
    if anchors:
        c.caption(f"{anchors['n']} reviewed, significant and reliable anchor lines: mean offset {anchors['mean_kms']:+.1f} km/s, scatter {anchors['sigma_kms']:.1f} km/s.")
        c.button('Use anchor calibration',on_click=use_anchors,disabled=anchors['n']<3,help='Needs at least three anchors.')
    else:c.caption('Review a few secure lines to calibrate σ and the mean offset from your own data.')
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
        st.caption(f'Each candidate’s probability assumes the feature is real: Gaussian in velocity offset (mean {id_mean:+.0f} km/s, σ {id_sigma:.0f} km/s combined with the centroid error), against an “unidentified” alternative with the weight of one line spread over ±{halfwidth:.0f} km/s. Combined = detection reliability × identification probability. Chance match is the probability that at least one unrelated catalogue line lies within ±2σ at this catalogue density; when it is high, a match on its own means little.')
match_settings.update(region_radius_arcsec=region_radius,min_reliability=min_reliability,include_low_reliability=include_low)
st.session_state['identification_matches']=matches
if identifications:
    st.caption('Live query time: '+identifications['queried_utc'])
    if identifications['status']:st.dataframe(pd.DataFrame(identifications['status']),hide_index=True,width='stretch')
    else:st.info('No significant features to match. Review provisional features or adjust the detection settings.')
    with st.expander('Catalogue filter audit'):
        st.json({k:a['reason_counts'] for k,a in identifications['audits'].items()})
        st.download_button('Download catalogue query and filter audit',json.dumps(identifications,indent=2),'live_catalogue_audit.json','application/json')
review_scores={ident:dict(reliability=scores.get(ident,{}).get('reliability'),region=scores.get(ident,{}).get('region'),method_note=method_notes.get(ident),**id_scores.get(ident,{})) for ident in lookup}
provenance=dict(queried_utc=identifications['queried_utc'] if identifications else None,method=method,settings=match_settings) if (identifications or method_suggestions is not None) else None
assignments,solution=identification_review(lookup,matches,solution_context,feature_keys,loaded,extractions,{r.ID:(int(r.Start),int(r.End)) for _,r in edited.iterrows()},provenance,revision,Path(os.environ.get('CONLINES_SOLUTION_DIR',ROOT/'solutions')),scores=review_scores,review_ids=[i for i in lookup if passes_reliability(i)],suggestions=method_suggestions,method=method)

st.subheader('5 · Make integrated maps')
adopted_ids=st.multiselect('Adopt features for mapping',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['sign']}",key='adopt_'+signature)
adopted=edited[edited.ID.isin(adopted_ids)]
if adopted.empty:st.info('Adopt at least one feature above to enable map export.');st.stop()
map_zoom=st.slider('Map preview zoom',1.,8.,3.,step=.5)
baseline=st.selectbox('Map baseline subtraction',['linear','constant','none'],help='Fits local sidebands per pixel. Inspect sidebands for blends.')
choices={ident:assignments[ident]['transition'] for ident in adopted_ids}
st.caption('Maps use the assignments from Step 4. Automatic suggestions remain tentative until reviewed; unidentified features can still be mapped.')
if st.button('Generate maps & download bundle',type='primary'):
    output=io.BytesIO();manifest=dict(detection=analysis,identifications=dict(identifications or {},matches=matches,settings=match_settings,method=method,amase=st.session_state.get('amase') if method=='amase' else None),solution=solution,apertures=configs,map_preview_zoom=map_zoom,lines=[])
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        for _,r in adopted.iterrows():
            k,q=lookup[r.ID];cube=loaded[k]['cube'];line=choices[r.ID]
            try:
                lo,hi=int(r.Start),int(r.End)
                im,h=moment(cube,lo,hi,line['rest_GHz'] if line else matching_frequency(q),z if line else 0.,baseline)
                label=line['species'] if line else 'Unidentified '+r.ID
                if line is None:
                    for keyword in ['RESTFRQ','RESTFREQ','REDSHIFT']:
                        if keyword in h:del h[keyword]
                    h['VREFGHZ']=(matching_frequency(q),'Observed reference for velocity channel widths')
                h['IDSTAT']='REVIEWED' if assignments[r.ID]['reviewed'] else 'TENTATIVE'
                h['LINEID']=label.encode('ascii','backslashreplace').decode();h['SRCFILE']=loaded[k]['name'].encode('ascii','backslashreplace').decode();h.add_history('Aperture: '+json.dumps(configs[k]))
                archive.writestr(r.ID+'.fits',fits_bytes(im,h))
                manifest['lines'].append(dict(id=r.ID,source=loaded[k]['name'],transition=line,lo=lo,hi=hi,baseline=baseline,frame=cube.frame,unit=h['BUNIT']))
                fig=chart_image(im,f"{loaded[k]['name']} · {label} · channels {lo}–{hi}",center=(configs[k]['x'],configs[k]['y']),zoom=map_zoom)
                fig.data[0].colorscale='RdBu_r';fig.data[0].zmid=0;fig.data[0].colorbar.title=h['BUNIT'];st.plotly_chart(fig,width='stretch',key='map_'+r.ID)
                archive.writestr(r.ID+'.html',fig.to_html(include_plotlyjs=True))
            except (ValueError,TypeError) as e:st.error(f'{r.ID}: {e}')
        archive.writestr('identification_solution.json',json.dumps(solution,indent=2));archive.writestr('manifest.json',json.dumps(manifest,indent=2));archive.writestr('features.csv',edited.to_csv(index=False))
    if manifest['lines']:
        st.download_button('Download FITS maps + interactive previews + provenance',output.getvalue(),'adopted_line_maps.zip','application/zip',on_click='ignore')
        st.success(f"Generated {len(manifest['lines'])} signed maps. Absorption retains its negative sign; no significance clipping was applied.")
