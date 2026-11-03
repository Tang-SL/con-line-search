"""Step 4 › Blended features: fit a feature with several transitions sharing one line shape."""
import copy
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from . import blends
from .solutions import digest, transition_key, describe_assignment

COLOURS = ['#c2410c', '#7356ad', '#0f766e', '#b45309', '#1d4ed8', '#be185d', '#4d7c0f', '#6b7280']


def component_key(ident, row):
    return f"{ident}|{row['species']}|{row['rest_GHz']:.7f}|{row.get('qns', '')}"


def fit_figure(fit, region_index, freq, spectrum, unit, names):
    model = fit['models'][region_index]
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[.74, .26], vertical_spacing=.04)
    fig.add_trace(go.Scatter(x=freq, y=spectrum, mode='lines', line=dict(color='#167f83', shape='hv'), name='Aperture spectrum'), 1, 1)
    fig.add_trace(go.Scatter(x=freq, y=model['total'], mode='lines', line=dict(color='#111827', width=2), name='Total fit'), 1, 1)
    for n, (key, values) in enumerate(model['components'].items()):
        fig.add_trace(go.Scatter(x=freq, y=np.array(values)+model['offset'], mode='lines',
                                 line=dict(color=COLOURS[n % len(COLOURS)], dash='dash'), name=names[key]), 1, 1)
    fig.add_trace(go.Scatter(x=freq, y=model['residual'], mode='lines', line=dict(color='#6b7280', shape='hv'), name='Residual', showlegend=False), 2, 1)
    fig.add_hline(y=0, line_color='#cbd5e1', row=2, col=1)
    fig.update_yaxes(title_text=unit, row=1, col=1); fig.update_yaxes(title_text='Residual', row=2, col=1)
    fig.update_xaxes(title_text='Observed frequency (GHz)', row=2, col=1)
    fig.update_layout(height=470, margin=dict(l=10, r=10, t=30, b=40), legend=dict(orientation='h', y=-.22, x=0, font=dict(size=10)))
    return fig


def blend_section(lookup, loaded, extractions, z, inventories, feature_keys, state, assignments):
    """lookup: ident -> (spw key, feature); state: the review state (its assignments are updated on request)."""
    st.markdown('**Blended features**')
    st.caption('A broad or odd-shaped feature may contain several lines. Here it is fitted with one copy of a shared line shape '
               'at each candidate transition (lines from one nucleus share its kinematics); only their strengths differ. '
               'The fit tells you which lines are needed and how much each contributes.')
    st.caption('Vibrationally excited lines (e.g. HCN v₂=1, E_up ≈ 1000 K) are only listed if the maximum upper-state energy in the catalogue filters '
               '(above, advanced) allows them; the default is 300 K.')
    with_candidates = [i for i in lookup if inventories.get(i)]
    if not with_candidates:
        st.info('No catalogue transitions to fit yet. Search Splatalogue (or run Astro AMASE) first.')
        return
    width = lambda i: lookup[i][1]['hi']-lookup[i][1]['lo']
    default = max(with_candidates, key=lambda i: (sum(r.get('inside', True) for r in inventories[i]) > 1, width(i)))
    if st.session_state.get('blend_feature') not in with_candidates:
        st.session_state['blend_feature'] = default
    feature = st.selectbox('Feature to decompose', with_candidates, key='blend_feature',
                           format_func=lambda i: f"{i} · {lookup[i][1]['sign']} · {width(i)+1} channels · {len(inventories[i])} candidate lines · now: {describe_assignment(assignments.get(i))}")
    others = [i for i in with_candidates if i != feature]
    with_key = 'blend_with_'+feature
    if any(v not in others for v in st.session_state.get(with_key, [])):
        st.session_state[with_key] = [v for v in st.session_state[with_key] if v in others]
    together = st.multiselect('Fit together with other features', others, key=with_key,
                              format_func=lambda i: f"{i} · {loaded[lookup[i][0]]['name']} · {describe_assignment(assignments.get(i))}",
                              help='Features fitted together share one line shape. Needed to tie a line to its partner in another feature or spectral window.')
    group = [feature, *together]
    gkey = digest(group)

    # Candidate lines.
    store = st.session_state.setdefault('blend_include', {})
    rows = []
    for g in group:
        current = assignments.get(g) or {}
        chosen = [current.get('transition'), *current.get('blend', [])]
        ranked = sorted(inventories[g], key=lambda r: (not r.get('inside', True), r.get('Eu_K') if r.get('Eu_K') is not None else 1e9))
        auto = {transition_key(r) for r in ranked[:4]}
        for r in inventories[g]:
            key = component_key(g, r)
            include = store[gkey].get(key) if gkey in store and key in store[gkey] else (
                any(c is not None and transition_key(c) == transition_key(r) for c in chosen) or transition_key(r) in auto)
            rows.append(dict(Fit=bool(include), Line=blends.label(r), Feature=g, Offset_kms=round(r.get('velocity_offset_kms', 0.), 1),
                             Eu_K=r.get('Eu_K'), log10_A=r.get('log10_A_s'), Inside='yes' if r.get('inside', True) else 'near edge', key=key))
    revision = st.session_state.get('blend_revision', 0)
    table = st.data_editor(pd.DataFrame(rows).drop(columns=['key']), hide_index=True, width='stretch',
                           key=f'blend_table_{gkey}_{revision}', disabled=['Line', 'Feature', 'Offset_kms', 'Eu_K', 'log10_A', 'Inside'],
                           column_config={'Fit': st.column_config.CheckboxColumn('Fit', help='Include this line in the fit'),
                                          'Offset_kms': st.column_config.NumberColumn('Δv from centroid (km/s)'),
                                          'Eu_K': st.column_config.NumberColumn('E_up (K)'), 'log10_A': st.column_config.NumberColumn('log₁₀ A')})
    store[gkey] = {r['key']: bool(f) for r, f in zip(rows, table['Fit'])}
    included = [(r, next(x for x in inventories[r['Feature']] if component_key(r['Feature'], x) == r['key'])) for r, f in zip(rows, table['Fit']) if f]
    names = {r['key']: f"{row['species']} {row.get('qns', '')}".strip() for r, row in included}

    a, b = st.columns([1.3, 1])
    profile = a.radio('Shared line shape', list(blends.PROFILES), format_func=blends.PROFILES.get, key='blend_profile',
                      help='Start with one Gaussian. Use two Gaussians for double-peaked or self-absorbed lines. A template copies the measured shape of another, unblended line.')
    template = None
    if profile == 'template':
        candidates = [i for i in lookup if (assignments.get(i) or {}).get('transition') and i not in group]
        if not candidates:
            b.warning('Identify a clean, unblended line first (in the review above); its shape is then available as a template.')
        else:
            pick = b.selectbox('Template line', candidates, key='blend_template',
                               format_func=lambda i: f"{i} · {assignments[i]['transition']['species']}")
            k, q = lookup[pick]
            template = blends.make_template(loaded[k]['cube'].frequency, extractions[k]['spectrum'],
                                            assignments[pick]['transition']['rest_GHz'], z, q['lo'], q['hi'])
    labels = [r['key'] for r, _ in included]
    b.caption('**Tie lines** whose strength ratio is known, e.g. HCN v₂=1 l=1e equal to l=1f (ratio 1). '
              'A tie to a line in another feature needs both features in the fit.')
    ties = b.data_editor(pd.DataFrame(st.session_state.get('blend_ties_'+gkey, pd.DataFrame(columns=['Line', 'Same as', 'Ratio']))),
                         num_rows='dynamic', hide_index=True, width='stretch', key=f'blend_ties_editor_{gkey}',
                         column_config={'Line': st.column_config.SelectboxColumn('Line', options=labels),
                                        'Same as': st.column_config.SelectboxColumn('Same strength as', options=labels),
                                        'Ratio': st.column_config.NumberColumn('× ratio', min_value=0.01, default=1.0)})
    with st.expander('Fit settings (advanced)'):
        c1, c2, c3 = st.columns(3)
        margin = c1.number_input('Extra channels on each side', 0, 200, 10, key='blend_margin', help='Line-free channels beside the feature anchor the baseline.')
        vmax = c2.number_input('Velocity search range ± (km/s)', 20., 2000., 300., step=20., key='blend_vmax')
        baseline = c3.checkbox('Fit a constant offset per feature', True, key='blend_baseline')

    tie_map, problems = {}, []
    for _, t in ties.dropna(subset=['Line', 'Same as']).iterrows():
        if t['Line'] == t['Same as'] or t['Line'] not in labels or t['Same as'] not in labels:
            problems.append('Each tie needs two different lines that are both in the fit.'); continue
        tie_map[t['Line']] = (t['Same as'], float(t['Ratio']) if pd.notna(t['Ratio']) else 1.)
    for problem in set(problems):
        st.warning(problem)

    regions = []
    for g in group:
        k, q = lookup[g]; cube = loaded[k]['cube']; n = len(cube.frequency)
        lo, hi = max(0, q['lo']-int(margin)), min(n-1, q['hi']+int(margin))
        regions.append(dict(id=g, spw=k, feature_key=feature_keys[g], freq=np.asarray(cube.frequency[lo:hi+1], float),
                            spectrum=np.asarray(extractions[k]['spectrum'][lo:hi+1], float),
                            noise=np.asarray(extractions[k]['noise'][lo:hi+1], float), lo=lo, hi=hi))
    components = []
    for r, row in included:
        region = group.index(r['Feature'])
        components.append(dict(key=r['key'], row=row, expected_GHz=row['rest_GHz']/(1+z), region=region,
                               sign=lookup[r['Feature']][1]['sign'], tie=tie_map.get(r['key'])))
    inputs = digest(dict(group=group, comps=[(c['key'], c['tie']) for c in components], profile=profile, z=z, margin=int(margin),
                         vmax=float(vmax), baseline=baseline, template=template,
                         data=[digest([float(x) for x in r['spectrum']]) for r in regions]))
    sigma0 = None
    m = lookup[feature][1].get('aperture_measurement', {})
    if m.get('reliable') and m.get('sigma_kms'):
        sigma0 = float(m['sigma_kms'])
    if st.button('Fit the blend', type='primary', disabled=not components or (profile == 'template' and template is None), key='blend_fit'):
        try:
            with st.spinner('Fitting…'):
                fit = blends.decompose(regions, components, profile, z, template, float(vmax), baseline, sigma0)
            for r, reg in zip(fit['regions'], regions):
                r['feature_key'] = reg['feature_key']
            fit['names'] = names
            st.session_state['blend_last'] = dict(inputs=inputs, fit=fit)
        except ValueError as e:
            st.error(f'The fit failed: {e}')
    last = st.session_state.get('blend_last')
    if not last:
        st.caption('Choose the lines to fit (the table starts with the lines inside the feature with the lowest upper-state energies) and press **Fit the blend**.')
        return
    if last['inputs'] != inputs:
        st.info('The lines or settings changed. Press **Fit the blend** again.')
        return
    fit = last['fit']
    st.caption(f"Shared shape: {blends.describe_shape(fit)} · reduced χ² {fit['reduced_chi2']:.2f} over {fit['n']} channels.")
    for i, reg in enumerate(regions):
        st.plotly_chart(fit_figure(fit, i, reg['freq'], reg['spectrum'], loaded[reg['spw']]['cube'].unit, names), width='stretch', key=f'blend_plot_{i}')
    unit = loaded[regions[0]['spw']]['cube'].unit
    result = pd.DataFrame([dict(Line=names[q['key']], Feature=fit['regions'][q['region']]['id'],
                                Peak=q['amplitude'], Peak_error=q['amplitude_error'], Flux=q['flux'], Flux_error=q['flux_error'],
                                Share=q['fraction'], Delta_BIC=q['delta_bic'], Verdict=q['verdict'],
                                Trades_flux_with=', '.join(names[d] for d in q['degenerate_with']),
                                Tied_to=(f"{names[q['tie'][0]]} × {q['tie'][1]:g}" if q['tie'] else '')) for q in fit['components']])
    st.dataframe(result, hide_index=True, width='stretch',
                 column_config={'Peak': st.column_config.NumberColumn(f'Peak ({unit})', format='%.4g'),
                                'Peak_error': st.column_config.NumberColumn('±', format='%.2g'),
                                'Flux': st.column_config.NumberColumn(f'Integrated ({unit} km/s)', format='%.4g'),
                                'Flux_error': st.column_config.NumberColumn('±', format='%.2g'),
                                'Share': st.column_config.NumberColumn('Share of feature', format='percent'),
                                'Delta_BIC': st.column_config.NumberColumn('ΔBIC', format='%.1f', help='How much worse the fit gets without this line (refitting the shape). >10: needed; 2–10: weak; <2: not needed.')})
    if fit.get('at_bounds'):
        st.warning('The line shape reached the edge of its allowed range ('+', '.join(fit['at_bounds'])+'), so it does not describe this line well. '
                   'Try two Gaussians or a template, or widen the velocity range (advanced).')
    if fit['reduced_chi2'] > 2:
        st.warning('The fit leaves significant structure (reduced χ² > 2): try two Gaussians, a template, or more candidate lines.')
    if any(q['verdict'] == 'Needed, but not separable' for q in fit['components']):
        st.warning('Some lines are needed together but cannot be told apart with this spectrum alone: the fit can move flux between them '
                   '(see “Trades flux with”). Tie one of them to a partner line in another feature or window, use a template, '
                   'or report their combined flux.')
    if any(q['tie'] for q in fit['components']):
        st.caption('“Set by tie”: the line’s strength follows its partner, so the data in its own window need not require it on their own (ΔBIC shown for reference).')
    st.caption('Uncertainties come from the full fit (line strengths and shape together) and the per-channel noise of the aperture controls; '
               'channel-to-channel correlations are ignored. Self-absorbed lines are fitted as shapes: their strengths are not optically thin '
               'line ratios.')

    def drop_unneeded():
        for q in fit['components']:
            if q['verdict'] == 'Not needed':
                st.session_state['blend_include'][gkey][q['key']] = False
        st.session_state['blend_revision'] = st.session_state.get('blend_revision', 0)+1

    def adopt():
        chosen = blends.assignments_from_fit(fit)
        stored = st.session_state.setdefault('blend_fits', {})
        record = copy.deepcopy(fit)
        for reg in record['regions']:
            stored[reg['feature_key']] = record
        for ident, a in chosen.items():
            state['assignments'][ident] = a
        state['widget_revision'] += 1
        st.session_state['blend_notice'] = 'Applied: ' + '; '.join(f"{i} = {describe_assignment(a)}" for i, a in chosen.items())

    c1, c2 = st.columns(2)
    c1.button('Drop the lines that are not needed', on_click=drop_unneeded, disabled=not any(q['verdict'] == 'Not needed' for q in fit['components']), key='blend_drop')
    c2.button('Use this fit as the identification', on_click=adopt, type='primary', key='blend_adopt',
              help='Main line = the strongest needed line in each feature; the other needed lines are recorded as blended. Marked as reviewed; maps of each line become available in Step 5.')
    if st.session_state.get('blend_notice'):
        st.success(st.session_state['blend_notice'])
