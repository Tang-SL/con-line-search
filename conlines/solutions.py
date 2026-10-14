"""Portable, explicit identification choices; saved solutions are not catalogues."""
import copy
import hashlib
import json
import math
from datetime import datetime, timezone
import numpy as np
from astropy.wcs import WCS


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def cube_fingerprint(cube):
    """Same scientific data loaded from another path receive the same identity."""
    h=hashlib.sha256()
    h.update(np.ascontiguousarray(cube.data,dtype='<f4').data)
    h.update(np.asarray(cube.frequency,dtype='<f8').tobytes())
    metadata=dict(shape=cube.data.shape,unit=cube.unit,frame=cube.frame,wcs=dict(WCS(cube.header).celestial.to_header()),beam={k:cube.header.get(k) for k in ['BMAJ','BMIN','BPA']})
    h.update(json.dumps(metadata,sort_keys=True).encode())
    return h.hexdigest()


def feature_fingerprint(cube_id,feature):
    return digest(dict(cube=cube_id,sign=feature['sign'],lo=feature['lo'],hi=feature['hi'],frequency_GHz=round(feature['frequency_GHz'],9),x=round(feature.get('x',0.),6),y=round(feature.get('y',0.),6)))


def transition_key(row):
    return None if row is None else (row['species'],round(float(row['rest_GHz']),9),row.get('qns',''))


def suggest_assignments(feature_ids,matches,previous=None):
    out={}
    for ident in feature_ids:
        options=matches.get(ident,[])
        old=(previous or {}).get(ident)
        if old and old.get('reviewed'):
            out[ident]=copy.deepcopy(old)
            updated=next((r for r in options if transition_key(r)==transition_key(old['transition'])),None)
            if updated is not None:out[ident]['transition']=copy.deepcopy(updated)
        else:
            nearest=min(options,key=lambda r:abs(r['velocity_offset_kms'])) if options else None
            out[ident]=dict(transition=copy.deepcopy(nearest),reviewed=False,origin='nearest' if nearest else 'unidentified')
    return out


def make_solution(name,context,feature_keys,assignments,intervals,catalogue_provenance=None):
    return dict(schema='conlines-identification-solution',version=1,name=name,
                saved_utc=datetime.now(timezone.utc).isoformat(),context=context,
                catalogue_provenance=catalogue_provenance,features=[
                    dict(feature_key=feature_keys[i],display_id=i,assignment=copy.deepcopy(assignments[i]),
                         lo=int(intervals[i][0]),hi=int(intervals[i][1])) for i in feature_keys])


def restore_solution(payload,context,feature_keys,channel_counts):
    """Reject mismatched analyses; never silently attach a choice to another line."""
    if not isinstance(payload,dict) or payload.get('schema')!='conlines-identification-solution' or payload.get('version')!=1:
        raise ValueError('Unsupported identification solution format.')
    if digest(payload.get('context'))!=digest(context):
        raise ValueError('Solution does not match the current cubes, aperture, redshift or detection settings. Restore those settings before loading it.')
    records=payload.get('features')
    if not isinstance(records,list):raise ValueError('Solution features must be a list.')
    by_key={}
    for record in records:
        if not isinstance(record,dict) or not isinstance(record.get('feature_key'),str):raise ValueError('Invalid solution feature.')
        if record['feature_key'] in by_key:raise ValueError('Duplicate solution feature.')
        by_key[record['feature_key']]=record
    if set(by_key)!=set(feature_keys.values()):raise ValueError('Saved features differ from the current detection. Rerun with the saved settings.')
    assignments={};intervals={}
    for ident,key in feature_keys.items():
        record=by_key[key];a=record.get('assignment')
        if not isinstance(a,dict) or type(a.get('reviewed')) is not bool or 'transition' not in a:raise ValueError('Invalid saved assignment.')
        row=a['transition']
        if row is not None:
            if not isinstance(row,dict) or not isinstance(row.get('species'),str) or not isinstance(row.get('qns',''),str):raise ValueError('Invalid saved transition.')
            freq=row.get('rest_GHz')
            if not isinstance(freq,(int,float)) or not math.isfinite(freq) or freq<=0:raise ValueError('Invalid saved rest frequency.')
        lo,hi=record.get('lo'),record.get('hi')
        if type(lo) is not int or type(hi) is not int or not 0<=lo<=hi<channel_counts[ident]:raise ValueError('Invalid saved integration interval.')
        assignments[ident]=copy.deepcopy(a);intervals[ident]=(lo,hi)
    return assignments,intervals
