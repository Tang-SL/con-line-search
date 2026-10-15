"""Catalogue compatibility evidence, without excitation or probability claims."""
import copy
import re
from .core import C
from .measurements import matching_frequency
from .solutions import transition_key


def species_key(row):
    # Keep isotopologues, conformers and excitation entries separate. Do not
    # infer a chemical family by deleting numbers or vibrational-state labels.
    return ('id', str(row['species_id'])) if row.get('species_id') is not None else ('label', re.sub(r'\s+', '', row['species']))


def add_evidence(matches, lookup, selected_catalogues, loaded, z,
                 velocity_tolerance=50., linked_apertures=True, analysed_spws=None):
    """Count distinct other features, never catalogue rows in the same blend.

    Each companion is evaluated at the candidate's velocity offset. Missing
    features are displayed without a penalty. Only successful live queries are
    passed in selected_catalogues; failures remain explicit in coverage metadata.
    """
    result = copy.deepcopy(matches)
    inventory = {}
    for key, rows in selected_catalogues.items():
        for row in rows:
            inventory.setdefault(species_key(row), {}).setdefault(key, {})[transition_key(row)] = row
    for ident, options in result.items():
        source_key, feature = lookup[ident]
        source_frame = loaded[source_key]['cube'].frame
        source_measurement = feature.get('aperture_measurement', {})
        for candidate in options:
            offset = C*(1-matching_frequency(feature)/(candidate['rest_GHz']/(1+z)))
            candidate['velocity_offset_kms'] = float(offset)
            entries = []; significant = set(); provisional = set()
            for key, transitions in inventory.get(species_key(candidate), {}).items():
                cube = loaded[key]['cube']; frequencies = cube.frequency
                for row in transitions.values():
                    expected = row['rest_GHz']/(1+z)*(1-offset/C)
                    if not min(frequencies) <= expected <= max(frequencies):
                        continue
                    selected = transition_key(row) == transition_key(candidate)
                    entry = dict(dataset_key=key, SPW=loaded[key]['name'], Species=row['species'],
                                 Transition=row.get('qns', ''), Rest_GHz=row['rest_GHz'],
                                 Expected_GHz=float(expected), Features='', Status='No ESPADA feature',
                                 Velocity_difference_kms=None)
                    if analysed_spws is not None and key not in analysed_spws:
                        entry['Status'] = 'Detection unavailable'
                    elif key != source_key and (cube.frame != source_frame or source_frame == 'UNKNOWN'):
                        entry['Status'] = 'Frame comparison unavailable'
                    elif key != source_key and not linked_apertures:
                        entry['Status'] = 'Apertures are not linked'
                    else:
                        compatible = []
                        for other_id, (other_key, other) in lookup.items():
                            if other_key != key:
                                continue
                            centre = matching_frequency(other)
                            delta = C*(1-centre/(row['rest_GHz']/(1+z)))-offset
                            lo, hi = sorted([frequencies[other['lo']], frequencies[other['hi']]])
                            half_channel = abs(frequencies[1]-frequencies[0])/2
                            if abs(delta) > velocity_tolerance or not lo-half_channel <= expected <= hi+half_channel:
                                continue
                            measured = other.get('aperture_measurement', {}).get('reliable', False)
                            if other_id == ident:
                                status = 'Selected feature' if selected else 'Same feature / possible blend'
                            elif selected:
                                status = 'Same transition in another feature'
                            elif other['sign'] != feature['sign']:
                                status = 'Opposite-sign feature'
                            elif not measured or not source_measurement.get('reliable'):
                                status = 'Uncertain aperture centroid'
                            elif other['significant']:
                                status = 'Significant feature, same sign'; significant.add(other_id)
                            else:
                                status = 'Provisional feature, same sign'; provisional.add(other_id)
                            compatible.append((other_id, status, float(delta)))
                        if compatible:
                            entry.update(Features=', '.join(c[0] for c in compatible),
                                         Status='; '.join(dict.fromkeys(c[1] for c in compatible)),
                                         Velocity_difference_kms=min(compatible, key=lambda c:abs(c[2]))[2])
                    entries.append(entry)
            candidate['evidence'] = dict(other_significant_features=sorted(significant),
                                         other_provisional_features=sorted(provisional-significant),
                                         transitions=entries, velocity_tolerance_kms=velocity_tolerance,
                                         queried_spws=[loaded[k]['name'] for k in selected_catalogues],
                                         unavailable_spws=[v['name'] for k,v in loaded.items() if k not in selected_catalogues])
    return result
