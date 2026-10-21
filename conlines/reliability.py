"""Two separate scores: is a feature real, and which identification is likely.

Detection reliability never uses the line catalogue, and identification
probabilities never change whether a feature counts as detected. Neither score
is a calibrated false-alarm probability; both are explicit, reproducible
estimates whose assumptions are listed next to each function.
"""
import math
import numpy as np
from .core import C


def channel_noise(cube):
    """Robust per-channel noise (1.4826 x MAD) over the full field.

    Source pixels are a small fraction of a full-field cube, so the median
    statistics are dominated by noise. One channel at a time keeps memory low.
    """
    noise = np.empty(len(cube.frequency))
    for i, plane in enumerate(cube.data):
        values = plane[np.isfinite(plane)]
        centre = np.median(values)
        noise[i] = 1.4826*np.median(np.abs(values-centre))
    return noise


def component_properties(cube, components, noise):
    """Catalogue-free strength of each ESPADA component.

    peak_snr: the largest signed, channel-integrated S/N over the pixels of the
    component's own SoFiA mask (sign follows emission/absorption). Channels are
    treated as independent, so this ranks components; it is not an aperture S/N.
    """
    ny, nx = cube.data.shape[1:]
    out = []
    for c in components:
        lo, hi = int(c['lo']), int(c['hi'])
        sign = 1. if c['sign'] == 'Emission' else -1.
        pixels = np.asarray(c.get('footprint') or [], int)
        if pixels.size == 0:
            ix = int(np.clip(round(c['x']), 0, nx-1)); iy = int(np.clip(round(c['y']), 0, ny-1))
            pixels = np.array([iy*nx+ix])
        block = cube.data[lo:hi+1].reshape(hi-lo+1, -1)[:, pixels]
        sigma = math.sqrt(float(np.sum(noise[lo:hi+1]**2)))
        snr = sign*np.nansum(block, axis=0)/sigma if sigma > 0 else np.zeros(pixels.size)
        out.append(dict(id=c['id'], sign=c['sign'], channels=hi-lo+1, pixels=int(pixels.size),
                        peak_snr=float(np.max(snr))))
    return out


def poisson_upper(k):
    """Approximate 84% (1 sigma) Poisson upper limit (Gehrels 1986)."""
    return k+1+math.sqrt(k+.75)


def detection_reliability(props, radius, region_radius, areas, near_factor=3.):
    """Empirical reliability from the artefact rate outside the source region.

    props: component_properties() rows. radius: each component's distance from
    the source centre (arcsec). areas: dict(inside, near, far) in pixels, where
    'near' is region_radius < r <= near_factor*region_radius and 'far' is beyond.

    For a component inside the region, count components of the same sign that
    are at least as strong in BOTH peak S/N and number of channels: n_obs inside,
    k_near and k_far in the two control zones. Each zone predicts
    k*A_inside/A_zone false components inside; the larger prediction is used, so
    artefacts concentrated around a bright nucleus are not underestimated.
    reliability = 1 - (84% upper limit of that prediction)/n_obs, clipped to [0, 1];
    reliability_best uses the counts themselves. Components outside the region
    are the control sample and receive no reliability.
    Assumes artefacts are spread evenly within each control zone and that
    real features outside the region are rare; both make the estimate cautious.
    """
    inside = [r <= region_radius for r in radius]
    near = [region_radius < r <= near_factor*region_radius for r in radius]
    far = [r > near_factor*region_radius for r in radius]
    out = []
    for p, ins in zip(props, inside):
        if not ins:
            out.append(dict(region='Control', reliability=None, reliability_best=None,
                            observed=None, control_near=None, control_far=None, expected_false=None))
            continue
        def at_least(q):
            return q['sign'] == p['sign'] and q['peak_snr'] >= p['peak_snr'] and q['channels'] >= p['channels']
        n_obs = sum(1 for q, i in zip(props, inside) if i and at_least(q))
        counts = {}
        predictions = []
        for zone, member in [('near', near), ('far', far)]:
            k = sum(1 for q, m in zip(props, member) if m and at_least(q))
            counts[zone] = k
            if areas.get(zone, 0) > 0:
                ratio = areas['inside']/areas[zone]
                predictions.append((k*ratio, poisson_upper(k)*ratio))
        if not predictions:
            out.append(dict(region='Source', reliability=None, reliability_best=None, observed=n_obs,
                            control_near=counts.get('near'), control_far=counts.get('far'), expected_false=None))
            continue
        best = max(p_[0] for p_ in predictions); upper = max(p_[1] for p_ in predictions)
        out.append(dict(region='Source', reliability=float(np.clip(1-upper/n_obs, 0, 1)),
                        reliability_best=float(np.clip(1-best/n_obs, 0, 1)), observed=n_obs,
                        control_near=counts['near'], control_far=counts['far'], expected_false=float(best)))
    return out


def neighbour_counts(offsets_arcsec, distance_arcsec):
    """Other components (any sign or window) within distance on the sky."""
    xy = np.asarray(offsets_arcsec, float).reshape(-1, 2)
    if not len(xy):
        return []
    d = np.hypot(*(xy[:, None, :]-xy[None, :, :]).transpose(2, 0, 1))
    return [int(v) for v in (d <= distance_arcsec).sum(axis=1)-1]


def identification_probabilities(offsets_kms, sigma_kms, mean_kms=0., window_kms=200.,
                                 centroid_error_kms=None, priors=None):
    """Probabilities over {each catalogue candidate, unidentified}.

    A candidate's weight is prior x Gaussian(offset; mean, s) with
    s^2 = sigma^2 + centroid_error^2. 'Unidentified' (a real line that is
    not one of the listed transitions) has the weight of one candidate spread
    uniformly over the +/-window. Probabilities assume the feature is real.
    """
    s = math.hypot(sigma_kms, centroid_error_kms or 0.)
    priors = priors or [1.]*len(offsets_kms)
    weights = [w*math.exp(-.5*((v-mean_kms)/s)**2)/(math.sqrt(2*math.pi)*s)
               for v, w in zip(offsets_kms, priors)]
    none = 1/(2*window_kms)
    total = sum(weights)+none
    return [w/total for w in weights], none/total


def chance_match_probability(line_count, coverage_GHz, frequency_GHz, sigma_kms, mean_kms=0., k=2.):
    """P(at least one unrelated catalogue line within mean +/- k*sigma).

    Treats the filtered catalogue lines as uniformly spread over the queried
    coverage. A high value means a match is expected by chance alone.
    """
    if coverage_GHz <= 0:
        return None
    density = line_count/coverage_GHz
    width = 2*k*sigma_kms/C*frequency_GHz
    return float(1-math.exp(-density*width))


def anchor_statistics(offsets_kms):
    """Mean and robust scatter of velocity offsets of reviewed anchor lines."""
    v = np.asarray([o for o in offsets_kms if o is not None and np.isfinite(o)], float)
    if len(v) < 2:
        return None
    centre = float(np.median(v))
    scatter = float(max(np.std(v, ddof=1), 1.4826*np.median(np.abs(v-centre))))
    return dict(n=int(len(v)), mean_kms=centre, sigma_kms=scatter)
