"""Signed aperture moments, conditional on the original detector interval."""
import numpy as np
from .core import C, residual


def matching_frequency(feature):
    measurement = feature.get('aperture_measurement', {})
    return measurement.get('frequency_GHz') if measurement.get('reliable') else feature['frequency_GHz']


def _scatter(values):
    return float(max(np.std(values, ddof=1), 1.4826*np.median(abs(values-np.median(values)))))


def measure_feature(cube, extraction, feature):
    """No positive clipping or Gaussian assumption; retain the ESPADA identity.

    Propagate whole control spectra through the moment derivatives, preserving
    their channel correlations. Add the constant/linear baseline difference in
    quadrature. Errors are conditional estimates, not calibrated confidence bounds.
    """
    out = dict(method='signed aperture moments v1', reliable=False,
               frequency_GHz=None, centroid_error_kms=None, sigma_kms=None,
               sigma_error_kms=None, flags=[], models=[])
    if feature.get('snr') is None or feature['snr'] < 3 or 'aperture sign mismatch' in feature.get('flags', []):
        out['flags'].append('Aperture S/N or sign does not support a centroid')
        return out
    lo, hi = feature['lo'], feature['hi']
    if hi-lo < 2:
        out['flags'].append('Fewer than three channels in the feature')
        return out
    f = np.asarray(cube.frequency[lo:hi+1], float)
    all_f = np.asarray(cube.frequency, float)
    edges = np.r_[all_f[0]-(all_f[1]-all_f[0])/2,
                  (all_f[1:]+all_f[:-1])/2, all_f[-1]+(all_f[-1]-all_f[-2])/2]
    df = abs(np.diff(edges))[lo:hi+1]
    sign = 1 if feature['sign'] == 'Emission' else -1
    values = np.vstack([extraction['spectrum'], extraction['controls']])
    try:
        for model in ['constant', 'linear']:
            r, _, complete = residual(values, lo, hi, model)
            if not complete:
                raise ValueError('Incomplete baseline coverage')
            controls = r[1:]-np.median(r[1:], axis=0)
            signal = sign*(r[0]-np.median(r[1:], axis=0))
            integral = np.sum(signal*df)
            if not np.isfinite(integral) or integral <= 0:
                raise ValueError('Non-positive signed aperture integral')
            # Centre the coordinates before taking moments to avoid cancellation.
            origin = f.mean()
            centre = origin+np.sum((f-origin)*signal*df)/integral
            if not f.min() <= centre <= f.max():
                raise ValueError('Aperture centroid lies outside the feature')
            variance = np.sum((f-centre)**2*signal*df)/integral
            centroid_error = _scatter(np.sum(controls*(f-centre)*df, axis=1)/integral)
            width = float(np.sqrt(variance)) if variance > 0 else None
            width_error = (_scatter(np.sum(controls*((f-centre)**2-variance)*df, axis=1)/integral)/(2*width)
                           if width is not None else None)
            out['models'].append(dict(baseline=model, frequency_GHz=float(centre),
                                      centroid_error_GHz=centroid_error, sigma_GHz=width,
                                      sigma_error_GHz=width_error))
    except ValueError as error:
        out['flags'].append(str(error))
        return out
    constant, linear = out['models']
    centre = linear['frequency_GHz']
    error = np.hypot(max(m['centroid_error_GHz'] for m in out['models']),
                     constant['frequency_GHz']-centre)
    out.update(frequency_GHz=centre, centroid_error_kms=float(C*error/centre), reliable=True)
    if error > abs(f[-1]-f[0])/2 or feature.get('variation', 0) is None or feature.get('variation', 0) > .3:
        out['reliable'] = False
        out['flags'].append('Centroid is too uncertain or baseline sensitive')
    if all(m['sigma_GHz'] is not None for m in out['models']):
        width_error = np.hypot(max(m['sigma_error_GHz'] for m in out['models']),
                              constant['sigma_GHz']-linear['sigma_GHz'])
        out.update(sigma_kms=float(C*linear['sigma_GHz']/centre),
                   sigma_error_kms=float(C*width_error/centre))
    else:
        out['flags'].append('Signed second moment does not give a positive width')
    return out
