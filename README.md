# con-line-search

Catalogue-guided spectral-line searches in ALMA cubes, starting with IC 860
Band 5 SPW0. This research prototype uses Splatalogue transition positions to
prepare targeted windows for ESPADA/SoFiA, then checks feature robustness and
retains competing molecular identifications.

**Current scope:** IC 860-specific scripts and a worked example. This is not yet
a general source-analysis package or an excitation/radiative-transfer model.

![IC 860 catalogue positions and nuclear spectrum](targets/ic860/alma/catalogue_spw0/qa/catalogue_overview.png)

## Example results

The pilot searched 18 windows in emission and absorption. Eight distinct
components were recovered; two absorption features pass the operational
robustness screen. One is consistent with H2S. The other has ambiguous
catalogue associations, including CH3OH and 30SiO. Feature detection does not
establish a molecular identification.

- [All prior positions and results](targets/ic860/alma/catalogue_spw0/qa/results.md)
- [Window spectra](targets/ic860/alma/catalogue_spw0/qa/targeted_windows.png)
- [Detailed method and limitations](targets/ic860/alma/catalogue_spw0/README.md)
- [Earlier four-SPW baseline analysis](targets/ic860/alma/espada_b5_all/README.md)

The +/-200 km/s windows cover 96% of SPW0; there is little search-bandwidth
reduction. S/N thresholds are robustness criteria, not calibrated false-alarm
probabilities. The published systemic redshift's exact reference frame remains
to be verified before interpreting small velocity offsets.

## Setup

Use Python 3.12 and a separate virtual environment:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install [SoFiA 2](https://gitlab.com/SoFiA-Admin/SoFiA-2) separately and put its
`sofia` executable on PATH. The pilot used SoFiA 2.7.1, commit
`a7434d9f91605cd738f2b13cd31b1450304fbc11`, with WCSLIB 8.2.2 and ESPADA 1.0.1.
ESPADA's SIP stage is disabled; no SIP installation is required for this workflow.

## Local data

Large FITS data, detector outputs, logs and environment files are excluded from
Git. The repository includes catalogue query snapshots, parameters, plots and
compact measurement results. The catalogue data remain third-party data from
[Splatalogue](https://splatalogue.online/#/basic); consult its upstream sources
and citation requirements when publishing results.

To use an existing CONs workspace without copying its cubes:

```sh
python scripts/link_ic860_data.py /path/to/CONs
```

This creates ignored local symlinks for the prepared common-beam cube directory
and the spectrum/noise archive. Alternatively, place the original Sabine Band 5
files under `targets/ic860/alma/Sabine/b5/`, then run:

```sh
python scripts/ic860_espada_b5_all_prepare.py
```

That preparation script needs the continuum map and paired continuum-subtracted
and +continuum cubes for all four SPWs. It writes the common-beam inputs and
noise spectra consumed by the SPW0 example. The archived all-SPW review also
requires its own ESPADA run; it is not needed for the catalogue pilot.

## Reproduce the catalogue-guided pilot

From the repository root, after local data setup:

```sh
# Optional: refresh the public catalogue queries (otherwise use saved snapshots).
python scripts/ic860_spw0_catalogue_query.py
python scripts/ic860_spw0_catalogue_prepare.py
MPLBACKEND=Agg espada -c targets/ic860/alma/catalogue_spw0/config/espada.yaml
python scripts/ic860_spw0_catalogue_review.py
```

Preparation regenerates absolute paths for the current checkout; no originating
machine paths are required. The review reuses baseline measurement functions
from `scripts/ic860_espada_b5_rank.py` without running its full four-SPW analysis.
The copied snapshots retain their original scientific provenance; regenerated
results can differ if the live catalogue or software versions change.

## Development direction

Separate source redshift/frame metadata, catalogue selection, targeted detection,
blend hypotheses, joint multi-transition physical models, and injection/recovery
validation. Retain non-detections and quantify upper limits under stated linewidth
and excitation assumptions. Do not lower detection thresholds solely because a
catalogue predicts a transition.

No redistribution license for the project's original code has been selected yet.

## Local web interface

```sh
source .venv/bin/activate
python -m pip install -r requirements-ui.txt
python -m streamlit run app.py
```

Open http://127.0.0.1:8501. The server binds to localhost; it is not deployed to
a public service. Cubes stay local. Live Splatalogue sends frequency ranges and
catalogue filter choices, not cube data. The saved SPW0 snapshot works offline.

1. Upload one or more FITS cubes, supply local paths (better for large files),
   or select the prepared IC 860 SPWs. Upload limit is 4 GB per file; RAM limits
   the actual batch size. Cubes are loaded in memory and up to eight are cached.
2. Choose a circle (default radius: one beam major-axis FWHM), ellipse,
   rectangle or single pixel. Choose the image centre, median-image brightness
   peak or custom zero-based pixel position. A shared sky position is projected
   onto each SPW's WCS. Angles are relative to each image's x axis; rotated grids
   need care for noncircular apertures. Spectra are aperture means in the cube's
   brightness units, not total integrated fluxes.
3. Inspect the aperture and noise-control positions. The detector needs at
   least eight non-overlapping translated apertures outside the central source;
   it uses at most 64. Source emission in these controls affects the noise
   estimate. Large apertures or small fields can fail this check explicitly.
4. Run the catalogue-guided search across the batch. Use **Live Splatalogue**
   for other SPWs; the saved snapshot is restricted to its actual frequency
   coverage. The default is Extragalactic, Eu <=300 K, +/-200 km/s, and an
   integrated S/N screen of 8. Hot Cores and other settings are adjustable.
5. Review significant and provisional candidates, inspect competing transitions,
   edit integration channels, and select which features to adopt. Adoption is
   manual; none is selected automatically. Choose the adopted transition for
   each feature before map generation.
6. Export a ZIP with celestial-WCS FITS moment-0 maps, interactive HTML previews,
   the candidate table and a JSON manifest containing source/aperture settings,
   selected catalogue priors, adopted transitions, intervals and baseline choices.
   Maps retain negative absorption and use positive channel widths in km/s.

### Interface detector and supported data

The interface uses a **new 1D aperture detector**, not the archived 3D ESPADA
backend. It searches 3/5/9-channel boxcars by default only around catalogue
positions. Two baseline estimates (constant, five sideband channels per side;
linear, ten per side) start after a two-channel gap. The same operation is
applied to the control spectra, retaining their spectral correlations. The
reported S/N is the smaller absolute value from the two models, requiring a
consistent sign. Significant candidates additionally require full sidebands,
two adjacent >=3-sigma channels, <=30% baseline variation, and no channel with
noise >3 times the median. Candidates >=3 S/N are available for manual review.

The search scans positions and widths and merges overlapping same-sign windows;
these S/N values are **not calibrated false-alarm probabilities**. No FDR,
injection-completeness or LTE abundance inference is claimed. Sidebands can be
contaminated by blends; inspect them and the map-baseline setting. Map generation
refits the chosen baseline per spatial pixel outside the edited line interval.
Manual interval changes do not silently recompute the original detection S/N.
Changing source/aperture/search inputs invalidates earlier detections and adoption.

Supported FITS ordering: RA, DEC, FREQ/VRAD/VOPT, optionally singleton Stokes.
Multiple polarization planes and other axis orders fail explicitly. Velocity
axes require RESTFRQ. Units: Jy/beam, mJy/beam (converted to Jy/beam), K, Jy/pixel.
Per-channel BEAMS tables are homogenized to a circular beam for each cube
(square pixels and finite data required); cubes are **not matched to a shared
beam across SPWs**. No frame conversion, primary-beam correction, spectral
stacking or cross-SPW flux combination is performed. Unknown/different reference
frames are shown explicitly. FITS maps preserve celestial WCS and beam metadata.

### Verification

```sh
python -m pytest -q
```

Synthetic tests exercise injected absorption, frequency/velocity conversion,
positive integration widths on reversed axes, linear-continuum subtraction,
invalid apertures, and FITS export. Streamlit tests exercise the empty app and,
when local IC 860 data are linked, detection, manual adoption, map generation,
and stale-result invalidation. The original ESPADA pilot remains unchanged.
