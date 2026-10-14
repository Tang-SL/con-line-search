"""Step 4: batch suggestions, sequential manual review and reusable solutions."""
from pathlib import Path
import json
import re
import uuid
import pandas as pd
import streamlit as st
from .plots import spectrum_plot
from .solutions import digest,transition_key,suggest_assignments,make_solution,restore_solution


def identification_review(lookup,matches,context,feature_keys,loaded,extractions,intervals,provenance,revision,solution_dir):
    scope=digest(dict(context=context,feature_keys=feature_keys));state_key='solution_'+scope;ids=list(lookup)
    if state_key not in st.session_state:
        st.session_state[state_key]=dict(assignments=suggest_assignments(ids,{}),cursor=None,revision=None,widget_revision=0,provenance=None)
    state=st.session_state[state_key]
    if revision and revision!=state['revision']:
        state['assignments']=suggest_assignments(ids,matches,state['assignments'])
        state.update(revision=revision,provenance=provenance,widget_revision=state['widget_revision']+1)
    assignments=state['assignments']
    if state['provenance']:
        st.caption('Assignment catalogue provenance: '+str(state['provenance'].get('queried_utc','saved solution'))+'. Reviewed choices are retained when live matches are refreshed.')
    st.caption('Nearest-frequency matches are tentative suggestions. Review them here without returning to Step 3. A question mark on the plot means the assignment has not been manually reviewed.')
    summary=[]
    for ident,(k,q) in lookup.items():
        a=assignments[ident];line=a['transition']
        summary.append(dict(Feature=ident,SPW=loaded[k]['name'],Sign=q['sign'],Significant=q['significant'],Assignment=line['species'] if line else 'Unidentified',Rest_GHz=line['rest_GHz'] if line else None,Offset_kms=line.get('velocity_offset_kms') if line else None,Review='Reviewed' if a['reviewed'] else 'Tentative'))
    st.dataframe(pd.DataFrame(summary),hide_index=True,width='stretch')
    reviewed=sum(a['reviewed'] for a in assignments.values())
    st.progress(reviewed/len(ids),text=f'{reviewed} of {len(ids)} features manually reviewed')
    show_spw=st.selectbox('SPW identification overview',list(dict.fromkeys(k for k,q in lookup.values())),format_func=lambda k:loaded[k]['name'],key='overview_'+scope)
    st.plotly_chart(spectrum_plot(loaded[show_spw]['cube'],extractions[show_spw],[q for k,q in lookup.values() if k==show_spw],assignments=assignments),width='stretch',key='assignplot_'+scope)

    def start_review():state['cursor']=0
    st.button('Review identifications one by one',on_click=start_review,key='start_'+scope)
    cursor=state['cursor']
    if cursor is not None and cursor<len(ids):
        ident=ids[cursor];k,q=lookup[ident];current=assignments[ident]
        with st.container(border=True):
            st.subheader(f'Feature {cursor+1} of {len(ids)} · {ident}')
            st.write(f"{loaded[k]['name']} · {q['sign']} · {q['frequency_GHz']:.7f} GHz · "+('Significant' if q['significant'] else 'Provisional'))
            members=list(matches.get(ident,[]))
            if current['transition'] is not None and not any(transition_key(r)==transition_key(current['transition']) for r in members):
                members.append(current['transition'])
                st.caption('The saved/current assignment is retained below. Fetch live matches to compare against the current catalogue filters.')
            selected=next((j for j,r in enumerate(members) if transition_key(r)==transition_key(current['transition'])),-1)
            choice_key=f'choice_{scope}_{state["widget_revision"]}_{ident}'
            choice=st.selectbox('Identification choice',[-1,*range(len(members))],index=selected+1,format_func=lambda j:'Unidentified / no reliable assignment' if j==-1 else f"{members[j]['species']} · {members[j]['rest_GHz']:.7f} GHz · Δv {members[j].get('velocity_offset_kms',0):+.1f} km/s · {members[j].get('qns','')}",key=choice_key)
            chosen=members[choice] if choice>=0 else None
            fig=spectrum_plot(loaded[k]['cube'],extractions[k],[q],selected_id=ident)
            f=loaded[k]['cube'].frequency;left,right=sorted([f[q['lo']],f[q['hi']]])
            pad=max((right-left),abs(f[1]-f[0])*8)
            fig.update_xaxes(range=[left-pad,right+pad])
            if chosen:
                fig.add_vline(x=chosen['rest_GHz']/(1+context['redshift']),line_color='#7356ad',line_dash='dash',annotation_text='Catalogue position',annotation_position='bottom right')
            st.plotly_chart(fig,width='stretch',key='reviewplot_'+ident)
            if members:
                with st.expander(f'Compare {len(members)} possible transitions'):
                    st.dataframe(pd.DataFrame(members)[[c for c in ['species','qns','rest_GHz','velocity_offset_kms','Eu_K','log10_A_s'] if c in members[0]]],hide_index=True,width='stretch')
            def confirm():
                value=st.session_state[choice_key]
                assignments[ident]=dict(transition=members[value] if value>=0 else None,reviewed=True,origin='manual')
                state['cursor']=cursor+1
            def previous():state['cursor']=max(0,cursor-1)
            a,b=st.columns(2)
            a.button('Previous feature',disabled=cursor==0,on_click=previous,key='prev_'+scope)
            b.button('Confirm & next' if cursor<len(ids)-1 else 'Confirm & finish',on_click=confirm,type='primary',key='confirm_'+scope)
    elif cursor is not None:
        st.success('Review complete. Your assignments are ready to save and use for maps.')

    st.markdown('**Save or reuse an identification solution**')
    name=st.text_input('Solution name',value='IC 860 identification' if any('IC 860' in v['name'] for v in loaded.values()) else 'Line identification',max_chars=100,key='solution_name_'+scope)
    payload=make_solution(name,context,feature_keys,assignments,intervals,state['provenance'])
    data=json.dumps(payload,indent=2,allow_nan=False)
    a,b=st.columns(2)
    a.download_button('Download solution (JSON)',data,'identification_solution.json','application/json',key='download_solution_'+scope)
    if b.button('Save solution locally',key='save_solution_'+scope):
        folder=Path(solution_dir);folder.mkdir(parents=True,exist_ok=True)
        slug=re.sub('[^a-zA-Z0-9_-]+','-',name).strip('-')[:60] or 'solution'
        path=folder/f'{slug}-{uuid.uuid4().hex[:8]}.json';path.write_text(data)
        st.success(f'Saved {path.name}')
    st.caption('Solutions store assignments, review status, integration channels and analysis settings. Reload them with the same cubes, aperture, redshift and detection settings. Reusing assignments does not query a cached catalogue; new candidate searches always use live Splatalogue.')
    with st.expander('Load a saved solution'):
        saved=sorted(Path(solution_dir).glob('*.json'),key=lambda p:p.stat().st_mtime,reverse=True)
        chosen_file=st.selectbox('Locally saved solutions',saved,format_func=lambda p:p.name,key='saved_solution_'+scope) if saved else None
        upload=st.file_uploader('Or upload a solution JSON',type=['json'],key='upload_solution_'+scope)
        if st.button('Load solution',disabled=upload is None and chosen_file is None,key='load_solution_'+scope):
            try:
                raw=upload.getvalue().decode('utf-8') if upload is not None else chosen_file.read_text()
                if len(raw)>5_000_000:raise ValueError('Solution exceeds 5 MB.')
                restored=json.loads(raw)
                aa,ii=restore_solution(restored,context,feature_keys,{ident:len(loaded[k]['cube'].frequency) for ident,(k,q) in lookup.items()})
                state.update(assignments=aa,cursor=None,widget_revision=state['widget_revision']+1,provenance=restored.get('catalogue_provenance'))
                st.session_state['solution_intervals_'+scope]=ii
                st.session_state['interval_revision_'+scope]=st.session_state.get('interval_revision_'+scope,0)+1
                st.session_state['solution_notice_'+scope]='Loaded solution: '+str(restored.get('name','Untitled'))
                st.rerun()
            except (ValueError,TypeError,KeyError,OSError) as e:st.error(f'Could not load solution: {e}')
    notice=st.session_state.get('solution_notice_'+scope)
    if notice:st.info(notice)
    return assignments,payload
