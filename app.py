"""Local, multi-SPW aperture spectroscopy workbench. Run with streamlit run app.py."""
from pathlib import Path
import hashlib, io, json, zipfile
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from astropy.wcs import WCS
from conlines.core import read_cube,extract,query_catalogue,select_priors,detect,moment,fits_bytes,C

ROOT=Path(__file__).resolve().parent
st.set_page_config(page_title='CON Line Search',page_icon='🔭',layout='wide')
st.markdown('''<style>.block-container{max-width:1450px;padding-top:2rem}h1{letter-spacing:-.04em} [data-testid="stMetric"]{background:#eaf0f5;padding:16px;border-radius:12px}</style>''',unsafe_allow_html=True)
st.caption('CON LINE SEARCH  /  SPECTRAL WORKBENCH')
st.title('From a cube to an adopted line')
st.write('Extract comparable apertures across spectral windows, inspect catalogue-guided candidates, and export signed integrated maps.')

@st.cache_resource(max_entries=8)
def load_path(path,mtime):return read_cube(path)
@st.cache_resource(max_entries=8)
def load_upload(digest,_file):
    _file.seek(0)
    return read_cube(_file)
@st.cache_data(ttl=86400,show_spinner=False)
def catalogue(low,high,filters):return query_catalogue(low,high,filters)

@st.cache_data(max_entries=32,show_spinner=False)
def measure(key,config,_cube):return extract(_cube,**config)
@st.cache_data(max_entries=8,show_spinner=False)
def preview_image(key,_cube):return np.nanmedian(_cube.data,axis=0)

def chart_image(image,title,mask=None,controls=None):
    fig=go.Figure(go.Heatmap(z=image,colorscale='Viridis',colorbar=dict(title='brightness')))
    if mask is not None:fig.add_trace(go.Contour(z=mask.astype(int),contours=dict(start=.5,end=.5,size=1,coloring='none'),line=dict(color='#ffb649',width=2),showscale=False,showlegend=False,hoverinfo='skip'))
    if controls:fig.add_trace(go.Scatter(x=[c[0] for c in controls],y=[c[1] for c in controls],mode='markers',marker=dict(size=4,color='white'),name='Noise controls'))
    fig.update_layout(title=title,height=430,margin=dict(l=10,r=10,t=45,b=10),xaxis_title='x pixel (zero-based)',yaxis_title='y pixel (zero-based)',showlegend=False,yaxis=dict(scaleanchor='x',constrain='domain',range=[-.5,image.shape[0]-.5]),xaxis=dict(constrain='domain',range=[-.5,image.shape[1]-.5]))
    return fig

def spectrum_plot(cube,extraction,features=()):
    fig=go.Figure(go.Scatter(x=cube.frequency,y=extraction['spectrum'],mode='lines',line=dict(color='#167f83',shape='hv'),name='Aperture mean'))
    for y,name in [(3*extraction['noise'],'+3σ controls'),(-3*extraction['noise'],'−3σ controls')]:fig.add_trace(go.Scatter(x=cube.frequency,y=y,line=dict(color='#b5bec8',dash='dot'),name=name))
    for q in features:fig.add_vrect(x0=min(cube.frequency[q['lo']],cube.frequency[q['hi']]),x1=max(cube.frequency[q['lo']],cube.frequency[q['hi']]),fillcolor='#dfac47',opacity=.2,line_width=0)
    fig.update_layout(height=350,xaxis_title=f'Observed frequency (GHz; {cube.frame})',yaxis_title=f'Aperture mean ({cube.unit})',margin=dict(l=10,r=10,t=45,b=50),legend=dict(orientation='h',y=1.15,x=0,font=dict(size=10)))
    return fig

with st.sidebar:
    st.subheader('Source & catalogue')
    preset=st.selectbox('Source preset',['IC 860','Custom source'])
    z=st.number_input('Redshift z',min_value=-.99,value=3880/C if preset=='IC 860' else 0.,format='%.7f',key='redshift_'+preset)
    st.caption('Use a redshift appropriate to the cube reference frame. The IC 860 preset is a working value; its precise frame remains to be verified.')
    filters=st.multiselect('Astronomical filters',['Extragalactic','Hot Cores'],default=['Extragalactic'])
    selected_filters=tuple('extragalactic' if v=='Extragalactic' else 'ism-hot-core' for v in filters)
    energy=st.number_input('Maximum upper-state energy (K)',1.,10000.,300.,step=50.)
    exclude=st.checkbox('Exclude deuterated / trans species',True)
    catalogue_source=st.selectbox('Catalogue source',['Saved IC 860 SPW0 snapshot','Live Splatalogue'])
    st.divider();st.subheader('Detection settings')
    halfwidth=st.number_input('Search ± velocity (km/s)',20.,2000.,200.,step=20.)
    threshold=st.number_input('Minimum integrated S/N',3.,50.,8.,step=1.)
    widths=st.multiselect('Trial widths (channels)',[3,5,7,9,13,17],default=[3,5,9])
    mode=st.selectbox('Search for',['Both','Emission','Absorption'])
    st.caption('Two local baselines, ≥2 adjacent channels at 3σ, and ≤30% baseline sensitivity. S/N is not a search-wide false-alarm probability.')

st.subheader('1 · Load spectral windows')
input_mode=st.radio('Input source',['Upload FITS cubes','Local file paths','IC 860 prepared example'],horizontal=True)
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
    base=ROOT/'targets/ic860/alma/espada_b5_all/inputs'
    chosen=st.multiselect('Example SPWs',[0,1,2,3],default=[0])
    for spw in chosen:
        p=base/f'spw{spw}_beam060.fits'
        try:loaded[f'example{spw}']=dict(name=f'IC 860 · SPW{spw}',cube=load_path(str(p),p.stat().st_mtime_ns))
        except Exception as e:errors.append(f'SPW{spw}: prepared local cube unavailable. See README data setup. ({e})')
for error in errors:st.error(error)
if not loaded:
    st.info('Add a cube to begin. You can also open the prepared IC 860 example without uploading files.')
    st.stop()
summary=[]
for k,v in loaded.items():
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
center_mode=st.radio('Centre',['Image centre','Brightness peak','Custom pixel position'],horizontal=True)
preview=preview_image(reference,ref)
if center_mode=='Brightness peak':
    # Useful for continuum cubes; inspect the preview for continuum-subtracted inputs.
    y0,x0=np.unravel_index(np.nanargmax(abs(preview)),preview.shape)
else:y0,x0=(np.array(ref.data.shape[1:])-1)/2
if center_mode=='Custom pixel position':
    a,b=st.columns(2);x0=a.number_input('Centre x (zero-based)',0.,float(ref.data.shape[2]-1),float(x0));y0=b.number_input('Centre y (zero-based)',0.,float(ref.data.shape[1]-1),float(y0))
linked=st.checkbox('Use the same sky position and angular aperture for every SPW',True)
sky=WCS(ref.header).celestial.pixel_to_world(x0,y0)
st.caption(f'Reference centre: {sky.to_string("hmsdms")} · Spectra are aperture means, not total fluxes. White dots show spatial noise-control apertures; inspect them for source contamination.')
extractions={};configs={}
for k,v in loaded.items():
    cube=v['cube']
    x,y=WCS(cube.header).celestial.world_to_pixel(sky) if linked else ((cube.data.shape[2]-1)/2,(cube.data.shape[1]-1)/2)
    if not linked:
        with st.expander(f"Centre for {v['name']}"):
            x=st.number_input('x pixel',0.,float(cube.data.shape[2]-1),float(x),key=k+'x');y=st.number_input('y pixel',0.,float(cube.data.shape[1]-1),float(y),key=k+'y')
    cfg=dict(x=float(x),y=float(y),kind=shape,radius=radius,ratio=ratio,angle=angle);configs[k]=cfg
    try:
        ex=measure(k,cfg,cube);extractions[k]=ex
        with st.expander(f"{v['name']} · {ex['pixels']} aperture pixels · {len(ex['controls'])} controls",expanded=k==reference):
            for note in cube.notes:st.info(note)
            col1,col2=st.columns([1,1.4]);col1.plotly_chart(chart_image(preview_image(k,cube),'Median-channel image',ex['mask'],ex['control_centers']),width='stretch',key=k+'image');col2.plotly_chart(spectrum_plot(cube,ex),width='stretch',key=k+'spectrum')
            csv=pd.DataFrame(dict(frequency_GHz=cube.frequency,aperture_mean=ex['spectrum'],control_noise=ex['noise'])).to_csv(index=False)
            st.download_button('Download extracted spectrum',csv,file_name=f'spectrum_{k[:12]}.csv',key=k+'csv')
    except ValueError as e:st.error(f"{v['name']}: {e}")

st.subheader('3 · Detect and review candidates')
st.caption('Catalogue-guided 1D aperture detector. This is a new backend, separate from ESPADA’s 3D source finder. Blank baselines and crowded lines still require inspection; a frequency match is not a molecular identification.')
signature=hashlib.sha256(json.dumps(dict(engine='aperture-v1.1',files=list(loaded),apertures=configs,z=z,filters=selected_filters,energy=energy,exclude=exclude,source=catalogue_source,halfwidth=halfwidth,threshold=threshold,widths=widths,mode=mode),sort_keys=True).encode()).hexdigest()
if st.button('Run detection across all loaded SPWs',type='primary',disabled=not extractions or not selected_filters or not widths):
    results={};status=[];prior_cache={}
    progress=st.progress(0.)
    for i,(k,ex) in enumerate(extractions.items()):
        cube=loaded[k]['cube'];low=cube.frequency.min()*(1+z)*(1-halfwidth/C);high=cube.frequency.max()*(1+z)*(1+halfwidth/C)
        try:
            with st.spinner(f"Querying and searching {loaded[k]['name']}…"):
                if catalogue_source=='Saved IC 860 SPW0 snapshot':
                    request=json.loads((ROOT/'targets/ic860/alma/catalogue_spw0/catalogue/extragalactic_request.json').read_text());limits=(float(request['frequenciesFrom'][0]),float(request['frequenciesTo'][0]))
                    if low<limits[0]-1e-7 or high>limits[1]+1e-7:raise ValueError('Saved catalogue does not cover this SPW/redshift. Select Live Splatalogue.')
                    rows=json.loads((ROOT/'targets/ic860/alma/catalogue_spw0/catalogue/normalized.json').read_text())
                else:rows=catalogue(float(low),float(high),selected_filters)
                priors=select_priors(rows,cube.frequency.min(),cube.frequency.max(),z,halfwidth,energy,selected_filters,exclude)
                if len(priors)>500:raise ValueError(f'{len(priors)} priors selected. Reduce energy/filter scope to at most 500 for this interactive search.')
                prior_cache[k]=priors
                found=detect(cube,ex,priors,halfwidth,threshold,tuple(widths),mode)
                for n,q in enumerate(found):q['id']=k[:12]+f'-{n+1:03}'
                results[k]=found;status.append(dict(SPW=loaded[k]['name'],Priors=len(priors),Candidates=len(found),Significant=sum(q['significant'] for q in found),Status='Completed'))
        except Exception as e:status.append(dict(SPW=loaded[k]['name'],Status=str(e)))
        progress.progress((i+1)/len(extractions))
    st.session_state['analysis']=dict(signature=signature,results=results,status=status,priors=prior_cache)
analysis=st.session_state.get('analysis')
if not analysis or analysis['signature']!=signature:
    if analysis:st.info('Inputs or settings changed. Rerun detection before adopting lines or making maps.')
    st.stop()
st.dataframe(pd.DataFrame(analysis['status']),hide_index=True,width='stretch')
rows=[];lookup={}
for k,features in analysis['results'].items():
    for q in features:
        lookup[q['id']]=(k,q)
        labels=list(dict.fromkeys(p['species'] for p in q['matches']))
        rows.append(dict(ID=q['id'],SPW=loaded[k]['name'],Matches=' / '.join(labels),Sign=q['sign'],SNR=round(q['snr'],2),Significant=q['significant'],Observed_GHz=q['frequency_GHz'],Start=q['lo'],End=q['hi'],Flags=', '.join(q['flags'])))
if not rows:st.info('No candidates reached the provisional S/N ≥3 reporting threshold. This is not an abundance constraint.');st.stop()
st.write('Select features to adopt below the table. Edit integration channels if needed; channel indices are zero-based and inclusive. Provisional candidates remain available for deliberate review.')
edited=st.data_editor(pd.DataFrame(rows),hide_index=True,width='stretch',key='review_'+signature,disabled=['ID','SPW','Matches','Sign','SNR','Significant','Observed_GHz','Flags'],column_config={'Start':st.column_config.NumberColumn(min_value=0,step=1),'End':st.column_config.NumberColumn(min_value=0,step=1)})
inspect=st.selectbox('Inspect candidate',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['species']} · S/N {lookup[ident][1]['snr']:.1f}")
k,q=lookup[inspect];st.plotly_chart(spectrum_plot(loaded[k]['cube'],extractions[k],[q]),width='stretch',key='inspect_'+inspect)
st.dataframe(pd.DataFrame([{pkey:p.get(pkey) for pkey in ['species','rest_GHz','qns','Eu_K','log10_A_s']} for p in q['matches']]),hide_index=True,width='stretch')
st.caption('Listed alternatives are from the selected catalogue filters and energy range. Absence from this table does not exclude other species or blends.')

st.subheader('4 · Make integrated maps')
adopted_ids=st.multiselect('Adopt features for mapping',list(lookup),format_func=lambda ident:f"{ident} · {lookup[ident][1]['species']} · S/N {lookup[ident][1]['snr']:.1f}",key='adopt_'+signature)
adopted=edited[edited.ID.isin(adopted_ids)]
if adopted.empty:st.info('Adopt at least one candidate above to enable map export.');st.stop()
baseline=st.selectbox('Map baseline subtraction',['linear','constant','none'],help='Fits local sidebands per spatial pixel, outside the adopted integration channels. Inspect sidebands for blends. None is appropriate only if residual continuum is negligible.')
choices={}
for _,r in adopted.iterrows():
    k,q=lookup[r.ID];members=q['matches']
    choice=st.selectbox(f'Adopted transition for {r.ID}',range(len(members)),format_func=lambda j,mm=members:f"{mm[j]['species']} · {mm[j]['rest_GHz']:.7f} GHz · {mm[j].get('qns','')}",key='line_'+r.ID+'_'+signature)
    choices[r.ID]=members[choice]
if st.button('Generate maps & download bundle',type='primary'):
    output=io.BytesIO();manifest=dict(redshift=z,settings=dict(halfwidth=halfwidth,threshold=threshold,widths=widths,filters=selected_filters,energy=energy,catalogue_source=catalogue_source),apertures=configs,catalogue_priors=analysis['priors'],lines=[])
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        for _,r in adopted.iterrows():
            k,q=lookup[r.ID];cube=loaded[k]['cube'];line=choices[r.ID]
            try:
                lo,hi=int(r.Start),int(r.End);im,h=moment(cube,lo,hi,line['rest_GHz'],z,baseline);h['LINEID']=line['species'].encode('ascii','backslashreplace').decode();h['SRCFILE']=loaded[k]['name'].encode('ascii','backslashreplace').decode();h.add_history('Aperture: '+json.dumps(configs[k]))
                archive.writestr(r.ID+'.fits',fits_bytes(im,h))
                meta=dict(id=r.ID,source=loaded[k]['name'],transition=line,lo=lo,hi=hi,baseline=baseline,detector_candidate=q,frame=cube.frame,unit=h['BUNIT']);manifest['lines'].append(meta)
                fig=chart_image(im,f"{loaded[k]['name']} · {line['species']} · channels {lo}–{hi}")
                fig.data[0].colorscale='RdBu_r';fig.data[0].zmid=0;fig.data[0].colorbar.title=h['BUNIT'];st.plotly_chart(fig,width='stretch',key='map_'+r.ID)
                archive.writestr(r.ID+'.html',fig.to_html(include_plotlyjs=True))
            except (ValueError,TypeError) as e:st.error(f'{r.ID}: {e}')
        archive.writestr('manifest.json',json.dumps(manifest,indent=2));archive.writestr('candidates.csv',edited.to_csv(index=False))
    if manifest['lines']:
        st.download_button('Download FITS maps + interactive previews + provenance',output.getvalue(),'adopted_line_maps.zip','application/zip',on_click='ignore')
        st.success(f"Generated {len(manifest['lines'])} signed maps. Absorption retains its negative sign; no significance clipping was applied.")
