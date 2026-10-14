# con-line-search

An ALMA spectral-line workbench with an IC 860 Band 5 example. The interactive
workflow detects emission/absorption features with ESPADA/SoFiA, then queries
live Splatalogue for possible identifications and exports adopted feature maps.
Historical catalogue-guided scripts and results are retained for reproducibility.
This research prototype is not an excitation/radiative-transfer model.

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
# Put the compiled SoFiA executable on PATH, or set SOFIA_EXECUTABLE.
python -m streamlit run app.py
```

Open http://127.0.0.1:8501. Cubes stay local. The interface now uses the real
ESPADA/SoFiA 3D finder first, followed by **live Splatalogue only** for identification.
Archived catalogue snapshots remain solely for reproducing the historical pilot.

1. Upload multiple FITS cubes, supply local paths, or select the IC 860 examples.
2. Choose a circle, ellipse, rectangle or single pixel and inspect the extracted
   aperture-mean spectra and spatial noise controls. Default radius is one beam
   major-axis FWHM. Shared sky positions are projected onto each SPW separately.
3. Run ESPADA on each cube, independently of redshift and molecular catalogues.
   Both emission and absorption are enabled by default. Components whose masks
   intersect the aperture are listed with short feature IDs. All are labelled
   on the spectrum; green features pass the additional aperture significance screen.
4. Fetch live Splatalogue matches for significant features (optionally provisional
   ones too). Every explicit matching request fetches fresh data. Network failure
   leaves detections available without substituting a saved catalogue. Changing
   catalogue filters or redshift does not invalidate detection or its significance.
5. Review possible identifications separately, sorted by absolute velocity offset
   within the adjustable default ±200 km/s tolerance. These are alternatives, not
   probability rankings or confirmed molecular identifications. Detection centroids
   come from the full ESPADA 3D component, not a new fit to the aperture spectrum.
6. Adopt features, optionally assign transitions, edit integration intervals and
   generate signed integrated maps. The default assignment is **Unidentified**.
   Unidentified maps use the observed centroid as the velocity-width reference;
   their FITS header has VREFGHZ instead of an assumed rest frequency/redshift.
   Previews default to 3× zoom; exported FITS maps retain the full field.

The ZIP export includes maps, interactive previews, feature table and provenance:
ESPADA parameter files/logs, aperture settings, live query time/results/filter audit,
and adopted transitions. Edited map intervals do not recompute detection significance.

### Detection and significance

ESPADA runs SoFiA in separate emission and absorption modes. Defaults: 5σ voxel
threshold; unsmoothed and 3/5/9-channel spectral kernels; unsmoothed and one-beam
spatial kernels. Source linking retains the pilot's minimum sizes (5 pixels per
spatial axis, 2 channels, 50 voxels). Thus spatially tiny sources can be missed.
The continuum-subtracted checkbox defaults on; turn it off for SoFiA's linear
continuum subtraction. Cached runs can be reused when only the aperture changes.

Each overlapping component is independently assessed in the chosen aperture.
Two local baseline estimates (constant with five sideband channels per side,
linear with ten, after a two-channel gap) are applied to both source and control
spectra. Reported S/N is the smaller absolute result and must agree with the
ESPADA sign. Significant features require S/N≥8, complete sidebands, two adjacent
channels at ≥3σ, ≤30% baseline variation, and no spectral-edge/high-noise flag.
All other ESPADA components remain available as provisional features. Large
apertures can dilute compact lines below the significance threshold.

These thresholds are **not search-wide false-alarm probabilities**. No FDR,
completeness calibration, LTE modelling or abundance inference is claimed.
Blends can contaminate sidebands and alter centroids. Identification requires
checking other transitions and physically consistent excitation.

### Catalogue filtering

Default Extragalactic selection, Eu≤300 K, deuterated/trans species excluded,
multiple explicit rare-isotope substitutions excluded, log10(A/s⁻¹)≥−5 required.
Hot Cores and all screening choices are adjustable. Single rare-isotope
substitutions and excited vibrational/torsional lines are retained by default.
These label-based screens are incomplete and can exclude real weak transitions;
Einstein A is not a flux prediction. Download the filter audit to inspect exclusions.
Catalogue settings affect identification only.

### Supported data and verification

FITS axis order: RA, DEC, FREQ/VRAD/VOPT, optionally singleton Stokes. Velocity
axes require RESTFRQ. Brightness units: Jy/beam, mJy/beam, K, Jy/pixel. Variable
BEAMS are homogenized per cube, not across SPWs. This ESPADA integration requires
finite pixels; crop blank borders first. At least eight translated noise-control
apertures must fit in the field. No frame conversion, primary-beam correction,
spectral stacking or cross-SPW flux combination is performed.

```sh
python -m pytest -q
```

Local-data UI integration tests require ESPADA and SoFiA. Tests cover catalogue-free
detection, independent identification settings, failed live queries, unidentified
map export, stale aperture invalidation, labels, spectral conversion and baselines.
