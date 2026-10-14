"""Local, multi-SPW aperture spectroscopy workbench. Run with streamlit run app.py."""
from pathlib import Path
import hashlib, io, json, zipfile, os
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import streamlit as st
from conlines.espada_backend import run_espada, aperture_features, match_features
from conlines.plots import chart_image, spectrum_plot
from conlines.solutions import cube_fingerprint, feature_fingerprint, digest
from conlines.review_ui import identification_review
from astropy.wcs import WCS
from conlines.core import read_cube,extract,query_catalogue,select_priors,moment,fits_bytes,C

ROOT=Path(__file__).resolve().parent
st.set_page_config(page_title='CON Line Search',page_icon='🔭',layout='wide')
st.markdown('''<style>.block-container{max-width:1450px;padding-top:2rem}h1{letter-spacing:-.04em} [data-testid="stMetric"]{background:#eaf0f5;padding:16px;border-radius:12px}</style>''',unsafe_allow_html=True)
st.caption('CON LINE SEARCH  /  SPECTRAL WORKBENCH')
st.title('From a cube to an adopted line')
st.write('Extract comparable apertures across spectral windows, detect emission and absorption with ESPADA, identify features with live Splatalogue, and export signed integrated maps.')

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
@st.cache_data(max_entries=8,show_spinner=False)
def preview_image(key,_cube):return np.nanmedian(_cube.data,axis=0)


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

st.subheader('2 · Choose an aperture and extract')
keys=list(loaded);reference=st.selectbox('Reference SPW',keys,format_func=lambda k:loaded[k]['name']);ref=loaded[reference]['cube']
a,b,c,d=st.columns(4)
with a:shape=st.selectbox('Aperture',['Circle','Ellipse','Rectangle','Single pixel'])
with b:radius=st.number_input('Radius / half-width (arcsec)',min_value=.0001,value=float(ref.beam_arcsec),format='%.4f',key='radius_'+reference)
with c:ratio=st.number_input('Minor / major ratio',.1,1.,1.,step=.1,disabled=shape in ['Circle','Single pixel'])
with d:angle=st.number_input('Angle from image x axis (deg)',-180.,180.,0.,step=5.,disabled=shape in ['Circle','Single pixel'])
if shape in ['Circle','Single pixel']:ratio=1.;angle=0.
centre_options=(['Continuum peak (example)'] if 'NUCX' in ref.header else [])+['Image centre','Brightness peak','Custom pixel position']
center_mode=st.radio('Centre',centre_options,horizontal=True)
preview=preview_image(loaded[reference]['fingerprint'],ref)
if center_mode=='Continuum peak (example)':y0,x0=ref.header['NUCY'],ref.header['NUCX']
elif center_mode=='Brightness peak':
    # Useful for continuum cubes; inspect the preview for continuum-subtracted inputs.
    y0,x0=np.unravel_index(np.nanargmax(abs(preview)),preview.shape)
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

st.subheader('3 · Detect and review features with ESPADA')
st.caption('ESPADA/SoFiA searches the 3D cube without a line catalogue. Components intersecting your aperture are screened using its extracted spectrum and spatial noise controls. Green features pass the integrated S/N and baseline checks; amber features require review. Neither threshold is a search-wide false-alarm probability.')
detection_settings=dict(engine='espada-v2-fullfield',threshold=threshold,widths=sorted(widths),mode=mode,clip=clip,continuum_subtracted=continuum_subtracted)
detection_context=dict(datasets=sorted([dict(fingerprint=v['fingerprint'],aperture={name:round(value,8) if isinstance(value,float) else value for name,value in configs[k].items()}) for k,v in loaded.items()],key=lambda v:v['fingerprint']),settings=detection_settings)
signature=digest(dict(context=detection_context,bindings=list(loaded)))
if st.button('Run detection across all loaded SPWs',type='primary',disabled=not extractions or not widths):
    results={};status=[];runs={}
    for i,(k,ex) in enumerate(extractions.items()):
        try:
            with st.spinner(f"Running ESPADA on {loaded[k]['name']}…"):
                run=search_cube(loaded[k]['fingerprint'],tuple(widths),mode,clip,continuum_subtracted,loaded[k]['cube'])
                found=aperture_features(loaded[k]['cube'],ex,run,threshold)
                prefix='SPW'+k.removeprefix('example') if k.startswith('example') else f'Cube{i+1}'
                for n,q in enumerate(found):q['id']=f'{prefix}-F{n+1:03}'
                results[k]=found;runs[k]={name:value for name,value in run.items() if name!='components'}
                status.append(dict(SPW=loaded[k]['name'],Components=len(run['components']),Aperture_features=len(found),Significant=sum(q['significant'] for q in found),Status='Completed'))
        except Exception as e:status.append(dict(SPW=loaded[k]['name'],Status=str(e)))
    st.session_state['analysis']=dict(signature=signature,results=results,status=status,runs=runs)
    st.session_state.pop('identifications',None)
analysis=st.session_state.get('analysis')
if not analysis or analysis['signature']!=signature:
    if analysis:st.info('Inputs or detection settings changed. Rerun detection before adopting features or making maps.')
    st.stop()
st.dataframe(pd.DataFrame(analysis['status']),hide_index=True,width='stretch')
rows=[];lookup={}
for k,features in analysis['results'].items():
    for q in features:
        lookup[q['id']]=(k,q)
        rows.append(dict(ID=q['id'],SPW=loaded[k]['name'],Sign=q['sign'],SNR=round(q['snr'],2) if q['snr'] is not None else None,Significant=q['significant'],Observed_GHz=q['frequency_GHz'],Start=q['lo'],End=q['hi'],Flags=', '.join(q['flags'])))
if not rows:st.info('No ESPADA components intersect this aperture. Try another aperture or inspect detection settings.');st.stop()
if not any(r['Significant'] for r in rows):st.info('No features pass the current aperture significance screen. ESPADA components remain below for provisional review; enable “Also match provisional features” to inspect their possible identifications.')
solution_context=dict(detection=detection_context,redshift=z)
feature_keys={ident:feature_fingerprint(loaded[k]['fingerprint'],q) for ident,(k,q) in lookup.items()}
solution_scope=digest(dict(context=solution_context,feature_keys=feature_keys))
restored_intervals=st.session_state.get('solution_intervals_'+solution_scope,{})
for row in rows:
    if row['ID'] in restored_intervals:row['Start'],row['End']=restored_intervals[row['ID']]
st.write('Feature IDs describe detections, not molecule assignments. Integration channels are zero-based and inclusive. Editing them changes maps only; the reported significance belongs to the original ESPADA interval.')
edited=st.data_editor(pd.DataFrame(rows).style.apply(lambda row: ['background-color: #d9f0e7; color: #075e46' if row['Significant'] else '' for _ in row],axis=1),hide_index=True,width='stretch',key='review_'+signature+'_'+str(st.session_state.get('interval_revision_'+solution_scope,0)),disabled=['ID','SPW','Sign','SNR','Significant','Observed_GHz','Flags'],column_config={'Start':st.column_config.NumberColumn(min_value=0,step=1),'End':st.column_config.NumberColumn(min_value=0,step=1)})
inspect=st.selectbox('Inspect feature',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['sign']}")
k,q=lookup[inspect];st.plotly_chart(spectrum_plot(loaded[k]['cube'],extractions[k],analysis['results'][k],selected_id=inspect),width='stretch',key='inspect_'+inspect)

st.subheader('4 · Identify features with live Splatalogue')
st.caption('Matches use measured feature centroids and your redshift/velocity tolerance. They are frequency-compatible alternatives, sorted by absolute velocity offset—not confirmed identifications or probability rankings. Other transitions and excitation consistency are needed to distinguish species and blends.')
include_provisional=st.checkbox('Also match provisional features',True)
match_settings=dict(detection=signature,z=z,tolerance=halfwidth,filters=selected_filters,energy=energy,exclude=exclude,min_log_a=min_log_a,exclude_multiple=exclude_multiple,exclude_excited=exclude_excited,include_provisional=include_provisional)
match_signature=hashlib.sha256(json.dumps(match_settings,sort_keys=True).encode()).hexdigest()
if st.button('Match features with live Splatalogue',disabled=not selected_filters):
    matches={};audits={};catalogues={};query_status=[]
    for k,features in analysis['results'].items():
        targets=[q for q in features if q['significant'] or include_provisional]
        if not targets:continue
        low=min(q['frequency_GHz'] for q in targets);high=max(q['frequency_GHz'] for q in targets)
        try:
            with st.spinner(f"Fetching live Splatalogue for {loaded[k]['name']}…"):
                raw=query_catalogue(low*(1+z)/(1+halfwidth/C),high*(1+z)/(1-halfwidth/C),selected_filters)
                selected,audit=select_priors(raw,low,high,z,halfwidth,energy,selected_filters,exclude,min_log_a,exclude_multiple,exclude_excited,return_audit=True)
                matches.update(match_features(targets,selected,z,halfwidth));audits[k]=audit;catalogues[k]=raw
                query_status.append(dict(SPW=loaded[k]['name'],Status='Live query completed',Transitions=len(selected)))
        except Exception as e:query_status.append(dict(SPW=loaded[k]['name'],Status=f'Live query failed: {e}'))
    st.session_state['identifications']=dict(signature=match_signature,matches=matches,audits=audits,catalogues=catalogues,status=query_status,queried_utc=datetime.now(timezone.utc).isoformat(),settings=match_settings)
identifications=st.session_state.get('identifications')
if identifications and identifications['signature']!=match_signature:
    st.info('Identification settings changed. Fetch live matches again. Feature detections remain valid.');identifications=None
matches=identifications['matches'] if identifications else {}
if identifications:
    st.caption('Live query time: '+identifications['queried_utc'])
    if identifications['status']:st.dataframe(pd.DataFrame(identifications['status']),hide_index=True,width='stretch')
    else:st.info('No significant features to match. Review provisional features or adjust the detection settings.')
    with st.expander('Catalogue filter audit'):
        st.json({k:a['reason_counts'] for k,a in identifications['audits'].items()})
        st.download_button('Download catalogue query and filter audit',json.dumps(identifications,indent=2),'live_catalogue_audit.json','application/json')
assignments,solution=identification_review(lookup,matches,solution_context,feature_keys,loaded,extractions,{r.ID:(int(r.Start),int(r.End)) for _,r in edited.iterrows()},dict(queried_utc=identifications['queried_utc'],settings=match_settings) if identifications else None,identifications['queried_utc'] if identifications else None,Path(os.environ.get('CONLINES_SOLUTION_DIR',ROOT/'solutions')))

st.subheader('5 · Make integrated maps')
adopted_ids=st.multiselect('Adopt features for mapping',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['sign']}",key='adopt_'+signature)
adopted=edited[edited.ID.isin(adopted_ids)]
if adopted.empty:st.info('Adopt at least one feature above to enable map export.');st.stop()
map_zoom=st.slider('Map preview zoom',1.,8.,3.,step=.5)
baseline=st.selectbox('Map baseline subtraction',['linear','constant','none'],help='Fits local sidebands per pixel. Inspect sidebands for blends.')
choices={ident:assignments[ident]['transition'] for ident in adopted_ids}
st.caption('Maps use the assignments from Step 4. Automatic suggestions remain tentative until reviewed; unidentified features can still be mapped.')
if st.button('Generate maps & download bundle',type='primary'):
    output=io.BytesIO();manifest=dict(detection=analysis,identifications=identifications,solution=solution,apertures=configs,map_preview_zoom=map_zoom,lines=[])
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        for _,r in adopted.iterrows():
            k,q=lookup[r.ID];cube=loaded[k]['cube'];line=choices[r.ID]
            try:
                lo,hi=int(r.Start),int(r.End)
                im,h=moment(cube,lo,hi,line['rest_GHz'] if line else q['frequency_GHz'],z if line else 0.,baseline)
                label=line['species'] if line else 'Unidentified '+r.ID
                if line is None:
                    for keyword in ['RESTFRQ','RESTFREQ','REDSHIFT']:
                        if keyword in h:del h[keyword]
                    h['VREFGHZ']=(q['frequency_GHz'],'Observed reference for velocity channel widths')
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
