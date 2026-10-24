"""Step 4: batch suggestions, sequential manual review and reusable solutions."""
from pathlib import Path
import copy
import json
import re
import uuid
import pandas as pd
import streamlit as st
from .plots import spectrum_plot
from .solutions import digest,transition_key,suggest_assignments,make_solution,restore_solution
from .measurements import matching_frequency
from .core import C


def evidence_panel(feature, chosen, members, loaded, extractions, ident, score=None):
    st.markdown('**Evidence for this identification**')
    score = score or {}
    if score.get('reliability') is not None:
        st.write(f"Detection reliability: **{score['reliability']:.2f}** (catalogue-independent; from the off-source artefact rate).")
    elif score.get('region') == 'Control':
        st.warning('This component lies outside the source region, where the app treats components as the artefact control sample.')
    if score.get('method_note'):
        st.write(score['method_note'])
    if score.get('chance') is not None:
        st.caption(f"Chance that at least one unrelated catalogue line falls this close: {score['chance']:.0%}. Probability that the line is none of the listed transitions: {score.get('p_unidentified', 0):.0%}.")
    measurement=feature.get('aperture_measurement',{})
    if measurement.get('reliable'):
        st.write(f"Aperture centroid: {measurement['frequency_GHz']:.7f} GHz · estimated error ±{measurement['centroid_error_kms']:.1f} km/s.")
        if measurement.get('sigma_kms') is not None:
            st.write(f"Aperture moment width σ: {measurement['sigma_kms']:.1f} ± {measurement['sigma_error_kms']:.1f} km/s.")
    else:
        st.info('Matching uses the ESPADA centroid because an aperture centroid could not be measured reliably. '+ '; '.join(measurement.get('flags',[])))
    st.caption('Moments use the original detection interval and both local baselines, preserving emission and absorption signs. Widths can be truncated by that interval or affected by blends. Errors are conditional estimates from noise controls and baseline changes; they exclude reference-frame/redshift uncertainty.')
    if chosen is None:
        st.caption('Choose a candidate to inspect its catalogue and companion-transition evidence.')
        return
    st.write(f"Velocity offset: {chosen['velocity_offset_kms']:+.1f} km/s relative to the entered systemic redshift.")
    fresh_row = next((r for r in members if transition_key(r) == transition_key(chosen)), {})
    if fresh_row.get('id_probability') is not None:
        combined = fresh_row.get('combined_probability')
        st.write(f"Identification probability (if real): **{fresh_row['id_probability']:.2f}**"+(f" · combined with detection reliability: **{combined:.2f}**" if combined is not None else ''))
    if fresh_row.get('amase_global') is not None:
        st.write(f"Astro AMASE scores: global **{fresh_row['amase_global']:.1f}** (threshold 93.5), local **{(fresh_row.get('amase_local') or 0):.2f}** (threshold 0.7); chemical relevance {(fresh_row.get('amase_structural') or 0):.1f}, frequency match {(fresh_row.get('amase_frequency') or 0):.3f}.")
        if fresh_row.get('amase_penalties'):
            st.caption('AMASE penalties: '+'; '.join(fresh_row['amase_penalties']))
    elif fresh_row.get('method_label'):
        st.caption(fresh_row['method_label'])
    st.caption(f"Catalogue: {chosen.get('catalogue') or 'not supplied'} · Species entry: {chosen.get('species_id') or 'not supplied'} · Upper energy: {chosen.get('Eu_K')} K · log₁₀ A: {chosen.get('log10_A_s')}. Energy and Einstein A alone do not predict observed brightness.")
    # Evidence is only taken from the current successful live request, never
    # from an assignment restored from a saved solution.
    fresh=next((r for r in members if transition_key(r)==transition_key(chosen)),None)
    evidence=fresh.get('evidence') if fresh else None
    if evidence is None:
        st.info('No current live companion evidence for this choice. Fetch live matches; saved assignments remain available.')
        return
    ns=len(evidence['other_significant_features']);np_=len(evidence['other_provisional_features'])
    st.write(f"Other compatible features: **{ns} significant**, **{np_} provisional**.")
    st.caption(f"Same catalogue species entry and line sign, within ±{evidence['velocity_tolerance_kms']:g} km/s of this candidate's velocity offset. Counts use distinct features, not the number of catalogue rows in a blend. Chance coincidences remain possible; this is not a probability or excitation fit.")
    if evidence['unavailable_spws']:
        st.warning('Live catalogue evidence unavailable for: '+', '.join(evidence['unavailable_spws']))
    entries=evidence['transitions']
    if entries:
        table=pd.DataFrame(entries).drop(columns=['dataset_key'])
        st.dataframe(table,hide_index=True,width='stretch')
        st.caption('“No ESPADA feature” is not an exclusion: a transition may be weak, blended or below the detection threshold. This inventory includes only transitions passing the current catalogue filters. Cross-SPW checks require linked apertures and matching known reference frames; line fluxes are not compared.')
        with st.expander('Inspect a predicted transition in its spectrum'):
            index=st.selectbox('Predicted transition',range(len(entries)),format_func=lambda i:f"{entries[i]['SPW']} · {entries[i]['Rest_GHz']:.7f} GHz · {entries[i]['Transition']} · {entries[i]['Status']}",key='companion_'+ident)
            entry=entries[index];key=entry['dataset_key']
            if entry['Status']=='Frame comparison unavailable':
                st.info('Align the spectral reference frames before comparing this predicted position.')
            elif key in extractions:
                cube=loaded[key]['cube'];expected=entry['Expected_GHz']
                fig=spectrum_plot(cube,extractions[key])
                pad=max(expected*300/C,abs(cube.frequency[1]-cube.frequency[0])*5)
                fig.update_xaxes(range=[expected-pad,expected+pad])
                fig.add_vline(x=expected,line_color='#7356ad',line_dash='dash',annotation_text='At candidate velocity')
                st.plotly_chart(fig,width='stretch',key='companionplot_'+ident)


def identification_review(lookup,matches,context,feature_keys,loaded,extractions,intervals,provenance,revision,solution_dir,scores=None,review_ids=None,suggestions=None,method='builtin'):
    scores=scores or {}
    scope=digest(dict(context=context,feature_keys=feature_keys));state_key='solution_'+scope;ids=list(lookup)
    # Reviewed choices belong to a detected component (its feature key), not to an
    # aperture, so they carry over when the aperture or redshift changes.
    store=st.session_state.setdefault('assignments_by_feature',{})
    if state_key not in st.session_state:
        carried={ident:store[feature_keys[ident]] for ident in ids if feature_keys[ident] in store}
        st.session_state[state_key]=dict(assignments=suggest_assignments(ids,{},carried),cursor=None,revision=None,widget_revision=0,provenance=None)
    state=st.session_state[state_key]
    if revision and revision!=state['revision']:
        state['assignments']=suggest_assignments(ids,matches,state['assignments'],suggestions,method)
        state.update(revision=revision,provenance=provenance,widget_revision=state['widget_revision']+1)
    assignments=state['assignments']
    for ident,(k,q) in lookup.items():
        line=assignments[ident]['transition']
        if line is not None:
            line['velocity_offset_kms']=float(C*(1-matching_frequency(q)/(line['rest_GHz']/(1+context['redshift']))))
    if state['provenance']:
        st.caption('Assignment catalogue provenance: '+str(state['provenance'].get('queried_utc','saved solution'))+'. Reviewed choices are retained when live matches are refreshed.')
    st.caption({'builtin':'Nearest-frequency matches','admit':'ADMIT choices','amase':'Astro AMASE assignments'}.get(method,'Suggestions')+' are tentative suggestions. Review them here without returning to Step 3. A question mark on the plot means the assignment has not been manually reviewed.')
    queue=[i for i in ids if review_ids is None or i in set(review_ids)]
    summary=[]
    for ident,(k,q) in lookup.items():
        a=assignments[ident];line=a['transition']
        score=scores.get(ident,{})
        fresh=next((r for r in matches.get(ident,[]) if line is not None and transition_key(r)==transition_key(line)),{})
        probability=fresh.get('id_probability') if line is not None else (score.get('p_unidentified') if matches.get(ident) else None)
        summary.append(dict(Feature=ident,SPW=loaded[k]['name'],Sign=q['sign'],Significant=q['significant'],Reliability=score.get('reliability'),Assignment=line['species'] if line else 'Unidentified',Rest_GHz=line['rest_GHz'] if line else None,Offset_kms=line.get('velocity_offset_kms') if line else None,P_assignment=probability,Combined=probability*score['reliability'] if probability is not None and score.get('reliability') is not None else None,Chance_match=score.get('chance'),Method_note=score.get('method_note'),Review='Reviewed' if a['reviewed'] else 'Tentative',In_review_queue=ident in queue))
    st.dataframe(pd.DataFrame(summary),hide_index=True,width='stretch')
    reviewed=sum(assignments[i]['reviewed'] for i in queue)
    if queue:st.progress(reviewed/len(queue),text=f'{reviewed} of {len(queue)} features in the review queue manually reviewed')
    if len(queue)<len(ids):st.caption(f'{len(ids)-len(queue)} features below the detection-reliability cut are left out of the review queue; they stay in the table and can still be mapped.')
    show_spw=st.selectbox('SPW identification overview',list(dict.fromkeys(k for k,q in lookup.values())),format_func=lambda k:loaded[k]['name'],key='overview_'+scope)
    st.plotly_chart(spectrum_plot(loaded[show_spw]['cube'],extractions[show_spw],[q for k,q in lookup.values() if k==show_spw],assignments=assignments),width='stretch',key='assignplot_'+scope)

    def start_review():state['cursor']=0
    st.button('Review identifications one by one',on_click=start_review,key='start_'+scope)
    cursor=state['cursor']
    if cursor is not None and cursor<len(queue):
        ident=queue[cursor];k,q=lookup[ident];current=assignments[ident]
        with st.container(border=True):
            st.subheader(f'Feature {cursor+1} of {len(queue)} · {ident}')
            st.write(f"{loaded[k]['name']} · {q['sign']} · {matching_frequency(q):.7f} GHz · "+('Significant' if q['significant'] else 'Provisional'))
            members=list(matches.get(ident,[]))
            if current['transition'] is not None and not any(transition_key(r)==transition_key(current['transition']) for r in members):
                members.append(current['transition'])
                st.caption('The saved/current assignment is retained below. Fetch live matches to compare against the current catalogue filters.')
            selected=next((j for j,r in enumerate(members) if transition_key(r)==transition_key(current['transition'])),-1)
            choice_key=f'choice_{scope}_{state["widget_revision"]}_{ident}'
            choice=st.selectbox('Identification choice',[-1,*range(len(members))],index=selected+1,format_func=lambda j:('Unidentified / no reliable assignment'+(f" · P {scores[ident]['p_unidentified']:.2f}" if matches.get(ident) and scores.get(ident,{}).get('p_unidentified') is not None else '')) if j==-1 else f"{members[j]['species']} · {members[j]['rest_GHz']:.7f} GHz · Δv {members[j].get('velocity_offset_kms',0):+.1f} km/s · {members[j].get('qns','')}"+(f" · P {members[j]['id_probability']:.2f}" if members[j].get('id_probability') is not None and members[j] in matches.get(ident,[]) else '')+(f" · {members[j]['method_label']}" if members[j].get('method_label') and members[j] in matches.get(ident,[]) else ''),key=choice_key)
            chosen=members[choice] if choice>=0 else None
            fig=spectrum_plot(loaded[k]['cube'],extractions[k],[q],selected_id=ident)
            f=loaded[k]['cube'].frequency;left,right=sorted([f[q['lo']],f[q['hi']]])
            pad=max((right-left),abs(f[1]-f[0])*8)
            fig.update_xaxes(range=[left-pad,right+pad])
            if chosen:
                fig.add_vline(x=chosen['rest_GHz']/(1+context['redshift']),line_color='#7356ad',line_dash='dash',annotation_text='Catalogue position',annotation_position='bottom right')
            st.plotly_chart(fig,width='stretch',key='reviewplot_'+ident)
            evidence_panel(q,chosen,list(matches.get(ident,[])),loaded,extractions,ident,scores.get(ident))
            if members:
                with st.expander(f'Compare {len(members)} possible transitions'):
                    comparison=[]
                    for member in members:
                        row={c:member.get(c) for c in ['species','qns','rest_GHz','velocity_offset_kms','Eu_K','log10_A_s','catalogue']}
                        fresh=next((r for r in matches.get(ident,[]) if transition_key(r)==transition_key(member)),{})
                        ev=fresh.get('evidence')
                        row['ID_probability']=fresh.get('id_probability');row['Combined']=fresh.get('combined_probability')
                        row['Other_significant_features']=len(ev['other_significant_features']) if ev else None
                        row['Other_provisional_features']=len(ev['other_provisional_features']) if ev else None
                        comparison.append(row)
                    st.dataframe(pd.DataFrame(comparison),hide_index=True,width='stretch')
            def confirm():
                value=st.session_state[choice_key]
                assignments[ident]=dict(transition=members[value] if value>=0 else None,reviewed=True,origin='manual')
                state['cursor']=cursor+1
            def previous():state['cursor']=max(0,cursor-1)
            a,b=st.columns(2)
            a.button('Previous feature',disabled=cursor==0,on_click=previous,key='prev_'+scope)
            b.button('Confirm & next' if cursor<len(queue)-1 else 'Confirm & finish',on_click=confirm,type='primary',key='confirm_'+scope)
    elif cursor is not None:
        st.success('Review complete. Your assignments are ready to save and use for maps.')

    st.markdown('**Save or reuse an identification solution**')
    name=st.text_input('Solution name',value='IC 860 identification' if any('IC 860' in v['name'] for v in loaded.values()) else 'Line identification',max_chars=100,key='solution_name_'+scope)
    payload=make_solution(name,context,feature_keys,assignments,intervals,state['provenance'])
    payload['measurement_method']='signed aperture moments v1'
    for record,record_id in zip(payload['features'],feature_keys):
        record['aperture_measurement']=lookup[record_id][1].get('aperture_measurement')
    payload['evidence_snapshot']=dict(queried_utc=provenance.get('queried_utc') if provenance else None,
                                    selected={feature_keys[ident]:next((r.get('evidence') for r in matches.get(ident,[]) if transition_key(r)==transition_key(a['transition'])),None) for ident,a in assignments.items()})
    data=json.dumps(payload,indent=2,allow_nan=False)
    a,b=st.columns(2)
    a.download_button('Download solution (JSON)',data,'identification_solution.json','application/json',key='download_solution_'+scope)
    if b.button('Save solution locally',key='save_solution_'+scope):
        folder=Path(solution_dir);folder.mkdir(parents=True,exist_ok=True)
        slug=re.sub('[^a-zA-Z0-9_-]+','-',name).strip('-')[:60] or 'solution'
        path=folder/f'{slug}-{uuid.uuid4().hex[:8]}.json';path.write_text(data)
        st.success(f'Saved {path.name}')
    st.caption('Solutions store assignments, review status, integration channels and analysis settings. Reload them with the same cubes, aperture, redshift and detection settings. Within a session, reviewed choices follow each detected component when you move the aperture. Reusing assignments does not query a cached catalogue; new candidate searches always use live Splatalogue.')
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
    # Distinct names: the confirm callback above reads `ident` when it runs.
    for feature_id,choice_record in assignments.items():store[feature_keys[feature_id]]=copy.deepcopy(choice_record)
    return assignments,payload
