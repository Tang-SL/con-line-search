"""Interactive scientific previews; viewport changes never crop exported data."""
import html
import numpy as np
import plotly.graph_objects as go
from .measurements import matching_frequency

def chart_image(image,title,mask=None,controls=None,center=None,zoom=1.):
    fig=go.Figure(go.Heatmap(z=image,colorscale='Viridis',colorbar=dict(title='brightness')))
    if mask is not None:fig.add_trace(go.Contour(z=mask.astype(int),contours=dict(start=.5,end=.5,size=1,coloring='none'),line=dict(color='#ffb649',width=2),showscale=False,showlegend=False,hoverinfo='skip'))
    if controls:fig.add_trace(go.Scatter(x=[c[0] for c in controls],y=[c[1] for c in controls],mode='markers',marker=dict(size=4,color='white'),name='Noise controls'))
    cx,cy=center if center is not None else ((image.shape[1]-1)/2,(image.shape[0]-1)/2)
    halfx=image.shape[1]/(2*zoom);halfy=image.shape[0]/(2*zoom)
    xmin=max(-.5,min(cx-halfx,image.shape[1]-.5-2*halfx));ymin=max(-.5,min(cy-halfy,image.shape[0]-.5-2*halfy))
    fig.update_layout(title=title,height=430,margin=dict(l=10,r=10,t=45,b=10),xaxis_title='x pixel (zero-based)',yaxis_title='y pixel (zero-based)',showlegend=False,yaxis=dict(scaleanchor='x',constrain='domain',range=[ymin,ymin+2*halfy]),xaxis=dict(constrain='domain',range=[xmin,xmin+2*halfx]))
    return fig

def spectrum_plot(cube,extraction,features=(),selected_id=None,assignments=None):
    fig=go.Figure(go.Scatter(x=cube.frequency,y=extraction['spectrum'],mode='lines',line=dict(color='#167f83',shape='hv'),name='Aperture mean'))
    for y,name in [(3*extraction['noise'],'+3σ controls'),(-3*extraction['noise'],'−3σ controls')]:fig.add_trace(go.Scatter(x=cube.frequency,y=y,line=dict(color='#b5bec8',dash='dot'),name=name))
    shown=set()
    for i,q in enumerate(sorted(features,key=lambda q:q['frequency_GHz'])):
        frequency=matching_frequency(q)
        significant=q['significant'];color='#087f5b' if significant else '#b38337'
        selected=q.get('id')==selected_id
        lo,hi=sorted([cube.frequency[q['lo']],cube.frequency[q['hi']]])
        fig.add_vrect(x0=lo,x1=hi,fillcolor=color,opacity=.22 if significant else .07,line_width=0)
        fig.add_vline(x=frequency,line_color=color,line_width=3 if selected or significant else 1,line_dash='solid' if significant else 'dot')
        label=q.get('id','feature')+' '+q.get('sign','')
        assignment=(assignments or {}).get(q.get('id'),{})
        line=assignment.get('transition')
        if assignments is not None:
            species=line['species'] if line else 'Unidentified'
            label=q.get('id','feature')+' '+(species[:23]+'…' if len(species)>24 else species)+(' ?' if not assignment.get('reviewed') else '')
        snr='unavailable' if q['snr'] is None else f"{q['snr']:.2f}"
        status='Significant' if significant else 'Provisional'
        idx=int(np.argmin(abs(cube.frequency-frequency)))
        method='Aperture centroid' if q.get('aperture_measurement',{}).get('reliable') else 'ESPADA centroid (fallback)'
        fig.add_trace(go.Scatter(x=[frequency],y=[extraction['spectrum'][idx]],mode='markers',marker=dict(size=11 if significant else 8,color=color,symbol='diamond' if significant else 'circle-open'),name=status,legendgroup=status,showlegend=status not in shown,hovertemplate=html.escape(label)+f'<br>{status} · S/N {snr}<br>{method}: %{{x:.6f}} GHz<extra></extra>'))
        shown.add(status)
        fig.add_annotation(x=frequency,y=1.02+(i%3)*.12,xref='x',yref='paper',text=('<b>'+html.escape(label)+'</b>') if significant else html.escape(label),showarrow=False,font=dict(color=color,size=10),textangle=-20,xanchor='center')
    fig.update_layout(height=440 if features else 350,xaxis_title=f'Observed frequency (GHz; {cube.frame})',yaxis_title=f'Aperture mean ({cube.unit})',margin=dict(l=10,r=10,t=120 if features else 45,b=50),legend=dict(orientation='h',y=-.25,x=0,font=dict(size=10)))
    return fig
