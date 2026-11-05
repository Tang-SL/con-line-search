"""Describe every cube: coverage, resolution, noise, blank pixels, continuum, brightest emission."""
import numpy as np
from astropy.wcs import WCS
from ..core import read_cube, C
from . import figs


def robust_sigma(values):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(1.4826*np.median(np.abs(values-np.median(values))))


def central_box(cube, fraction=.5):
    ny, nx = cube.data.shape[1:]
    hy, hx = int(ny*fraction/2), int(nx*fraction/2)
    return slice(ny//2-hy, ny//2+hy), slice(nx//2-hx, nx//2+hx)


def describe(record):
    cube = read_cube(record['path'])
    f = cube.frequency; nchan = len(f)
    df = abs(f[1]-f[0])
    sy, sx = central_box(cube)
    centre = cube.data[:, sy, sx]
    noise = np.array([robust_sigma(plane) for plane in centre])          # per channel, central half of the field
    finite = float(np.isfinite(cube.data[nchan//2]).mean())
    # Continuum: a real continuum raises every channel, including both band edges, while line
    # wings rarely reach both edges at once. The estimate is the lower of the two band-edge means
    # (first and last ~8% of channels). Percentiles over channels are not used: wherever lines
    # fill part of the band they shift every quantile upwards and mimic continuum. Significance
    # is measured against how that map scatters across the (mostly empty) field.
    median_map = np.full(cube.data.shape[1:], np.nan)
    ne = max(5, nchan//12)
    low_map = np.minimum(np.nanmean(centre[:ne], axis=0), np.nanmean(centre[-ne:], axis=0))
    median_map[sy, sx] = low_map
    spread = robust_sigma(low_map)
    level = float(np.nanmedian(low_map))
    cy, cx = np.unravel_index(np.nanargmax(np.nan_to_num(low_map-level, nan=-np.inf)), low_map.shape)
    continuum_snr = float((low_map[cy, cx]-level)/spread) if spread > 0 else None
    cy += sy.start; cx += sx.start
    # Brightest line emission: peak |S/N| over channels after removing the continuum estimate.
    snr_cube = (centre-low_map)/noise[:, None, None]
    peak_snr_map = np.nanmax(np.abs(snr_cube), axis=0)
    ly, lx = np.unravel_index(np.nanargmax(np.nan_to_num(peak_snr_map)), peak_snr_map.shape)
    lchan = int(np.nanargmax(np.abs(snr_cube[:, ly, lx])))
    ly += sy.start; lx += sx.start
    w = WCS(cube.header).celestial
    ny, nx = cube.data.shape[1:]
    centre_world = w.pixel_to_world((nx-1)/2, (ny-1)/2)
    line_world = w.pixel_to_world(lx, ly)
    info = dict(key=record['key'], path=record['path'], object=str(cube.header.get('OBJECT', '')) or record.get('target'),
                frame=cube.frame, unit=cube.unit, shape=list(cube.data.shape),
                freq_min_GHz=float(f.min()), freq_max_GHz=float(f.max()), nchan=nchan,
                channel_MHz=float(df*1e3), channel_kms=float(C*df/f.mean()),
                beam_major_arcsec=float(cube.header['BMAJ']*3600), beam_minor_arcsec=float(cube.header['BMIN']*3600),
                beam_pa_deg=float(cube.header.get('BPA', 0.)), pixel_arcsec=float(max(cube.pixel_arcsec)),
                finite_fraction=finite, noise_median=float(np.nanmedian(noise)), noise_per_channel=noise,
                continuum_peak=float(median_map[cy, cx]-level), continuum_snr=continuum_snr,
                continuum_pixel=[int(cx), int(cy)], line_peak_pixel=[int(lx), int(ly)], line_peak_channel=lchan,
                line_peak_GHz=float(f[lchan]), line_peak_snr=float(np.nanmax(peak_snr_map)),
                centre_ra=float(centre_world.ra.deg), centre_dec=float(centre_world.dec.deg),
                line_peak_ra=float(line_world.ra.deg), line_peak_dec=float(line_world.dec.deg),
                notes=list(cube.notes))
    flags = []
    if continuum_snr is not None and continuum_snr > 10:
        flags.append('continuum present: the cube does not look continuum-subtracted (continuum peak %.0f sigma above the field); '
                     'use --no-contsub when detecting, or subtract it first' % continuum_snr)
    if finite < .99:
        flags.append(f'{100*(1-finite):.0f}% of each channel is blank (primary-beam cutoff); the search uses the area with data')
    edge = noise[[0, -1]]
    if np.nanmax(edge) > 2*np.nanmedian(noise):
        flags.append('noisier channels at the band edges')
    info['flags'] = flags
    return cube, info, median_map, _full(peak_snr_map, cube, sy, sx)


def _full(sub, cube, sy, sx):
    out = np.full(cube.data.shape[1:], np.nan)
    out[sy, sx] = sub
    return out


def run(run, z=None):
    records = run.cubes()
    rows = []
    for key, record in records.items():
        cube, info, median_map, peak = describe(record)
        lx, ly = info['line_peak_pixel']
        figs.inspection(run.plot(f'inspect_{key}'), cube, info, median_map, peak, z)
        rows.append(info)
        del cube
    figs.coverage(run.plot('coverage'), rows, z)
    run.save('inventory', dict(cubes=rows))
    run.log('inspect', f'{len(rows)} cubes described', z=z)
    return rows


def summary(rows):
    lines = []
    for r in rows:
        lines.append(f"{r['key']}: {r['freq_min_GHz']:.3f}-{r['freq_max_GHz']:.3f} GHz, {r['nchan']} ch x {r['channel_kms']:.1f} km/s, "
                     f"beam {r['beam_major_arcsec']:.2f}x{r['beam_minor_arcsec']:.2f}\", noise {r['noise_median']:.2e} {r['unit']}/ch; "
                     f"brightest line S/N {r['line_peak_snr']:.0f} at {r['line_peak_GHz']:.4f} GHz, pixel {r['line_peak_pixel']}; "
                     f"continuum S/N {r['continuum_snr']:.0f}" + ''.join(f"\n    ! {f}" for f in r['flags']))
    return '\n'.join(lines)
