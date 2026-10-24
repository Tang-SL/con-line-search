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

That links the **historical cropped pilot** inputs. The current interface uses
separate full-field examples. Prepare them from the original Sabine Band 5 data:

```sh
python scripts/ic860_fullfield_prepare.py /path/to/CONs/targets/ic860/alma/Sabine/b5
```

This writes four ignored `targets/ic860/alma/fullfield/inputs/` FITS cubes,
retaining all 1024 × 1024 spatial pixels and convolving each channel to the same
0.060″ circular beam used in the pilot. Originals and historical 301 × 301
cutouts are untouched. The continuum map supplies the default nucleus position.
The preparation command refuses to overwrite existing outputs; `--spws 0` can
prepare one SPW at a time. Allow about 2.6 GB disk space for all four outputs.

To regenerate the archived pilot instead, place originals under
`targets/ic860/alma/Sabine/b5/` and run
`python scripts/ic860_espada_b5_all_prepare.py`.

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

### Starting and stopping (macOS)

Double-click **CON Line Search.app** in the project folder. It starts the app in
the background and opens it in your browser once it is ready; no terminal window
needs to stay open. You can drag it to the Dock. A copy in /Applications also
works once the app has been started from the project folder.

- If the app is already running, the icon asks whether to **Open** or **Stop** it.
- **Start app.command** does the same in a Terminal window and shows progress.
  **Stop app.command** stops the app.
- On the first start, a missing `.venv` is created with Python 3.12 or newer, and
  `requirements-ui.txt` is installed (a few minutes, needs internet). Starting
  from the app icon opens Terminal for this step so you can follow it.
- Install SoFiA 2 separately (see Setup) and place the compiled `sofia` executable
  in `.local-tools/`, or set `SOFIA_EXECUTABLE`. Without it the app opens, but
  detection is unavailable.

The same launcher works from a shell on macOS or Linux:

```sh
scripts/launch.sh start    # set up if needed, start in the background, open the browser
scripts/launch.sh open     # open the running app
scripts/launch.sh stop
scripts/launch.sh status   # exit code 0 if running
```

The server log is `.local-tools/app.log`. `CONLINES_PORT` changes the port
(default 8501). The address `http://127.0.0.1:8501` points to your own computer,
so it works only while the app is running; after a restart, start it again.

To run in the foreground instead:

```sh
source .venv/bin/activate
python -m pip install -r requirements-ui.txt
python -m streamlit run app.py
```

Open http://127.0.0.1:8501. Cubes stay local. The interface now uses the real
ESPADA/SoFiA 3D finder first, followed by **live Splatalogue only** for identification.
Archived catalogue snapshots remain solely for reproducing the historical pilot.

1. Upload multiple FITS cubes, supply local paths, or select the IC 860 full-field examples.
2. Run ESPADA on each full cube, independently of any aperture, redshift or
   molecular catalogue. Both emission and absorption are enabled by default.
   Every component in the field is shown on a map and in a table, with a feature
   ID that belongs to the cube's detection (for example `SPW0-F003`) and stays
   fixed when the aperture changes. Detection reruns only when the cubes or the
   finder settings (voxel threshold, trial widths, mode, continuum flag) change.
3. Choose a circle, ellipse, rectangle or single pixel and inspect the extracted
   aperture-mean spectra and spatial noise controls. Default radius is one beam
   major-axis FWHM. The centre can also be set on any detected component.
   Shared sky positions are projected onto each SPW separately.
   The example defaults to the continuum nucleus position. Step 3 has a 3× preview
   zoom; extraction and noise estimation always use the full loaded field. Controls
   are integer-pixel translations of the exact source stencil, avoiding boundary
   round-off changing their pixel counts. Set zoom to 1 to inspect the full field.
   Components whose masks intersect the aperture are screened automatically: all
   are labelled on the spectrum, and green features pass the additional aperture
   significance screen. Moving the aperture re-screens without rerunning ESPADA.
4. Fetch live Splatalogue matches for all detected features by default (uncheck the provisional
   option to restrict to significant ones). Every explicit matching request fetches fresh data. Network failure
   leaves detections available without substituting a saved catalogue. Changing
   catalogue filters or redshift does not invalidate detection or its significance.
   Moving the aperture does not need a new query: matches and companion evidence
   are recomputed from the same live response for the new aperture features.
   Reviewed choices follow each detected component within a session.
   Queries now cover each loaded SPW's full frequency range, including SPWs with
   no detections. Species IDs, catalogue names, lower energies and the original
   frequency text are retained alongside the existing transition metadata.
5. Step 4 proposes the closest frequency match for every queried feature, within
   the adjustable default ±200 km/s tolerance. The overview table and spectrum
   label these **tentative** suggestions. Use **Review identifications one by one**
   to choose a transition or leave a feature unidentified, then **Confirm & next**.
   This panel has its own spectrum and does not depend on the Step 3 inspector.
   Reviewed choices survive catalogue refreshes; a choice absent from the current
   filtered catalogue is kept and identified as such in the review panel.
   The **Evidence for this identification** panel shows the aperture centroid,
   its estimated error, the signed second-moment width, catalogue metadata and
   other transitions of the same catalogue species entry. Companion matches must
   lie inside another feature's interval and agree with the candidate's velocity
   offset within the adjustable default 50 km/s. Significant and provisional
   features are counted separately; multiple transitions in one feature count
   once. Opposite-sign features and uncertain centroids are shown separately.
   A companion spectrum viewer marks the predicted position at that velocity.
   These are frequency-compatibility checks, not identification probabilities.
6. Download a named solution as JSON or save it locally (ignored `solutions/`).
   Load a local solution or upload its JSON later. Solutions include assignments,
   review status, integration intervals, data fingerprints, aperture/detection
   settings, redshift and original live-query provenance. Only matching analyses
   are accepted; files may move to another path because identity is based on data
   and WCS, not filenames. Reusing chosen assignments is distinct from querying
   a catalogue: fresh alternatives always come from live Splatalogue.
7. Select features for maps in Step 5. It uses Step 4's assignments automatically,
   with no repeated transition selection. Unidentified maps use the observed
   centroid as the velocity-width reference, recording VREFGHZ instead of an
   assumed rest frequency/redshift. IDSTAT records TENTATIVE or REVIEWED.
   Previews default to 3× zoom; FITS maps retain the full spatial field.

Nearest frequency is not a physical identification or a probability ranking.
Matching now uses signed aperture moments inside the original ESPADA interval
when measurable. The original 3D centroid is retained for stable feature identity
and is an explicitly labelled fallback when aperture S/N is below 3, the line
sign disagrees, baseline coverage is incomplete, or the centroid is unstable.
Moment errors propagate whole control spectra (including channel correlations)
and constant/linear baseline differences. They are conditional estimates, not
calibrated confidence intervals, and omit systemic-redshift/frame uncertainty.
The width is a second-moment sigma within the detector interval, not a fitted
Gaussian FWHM or an instrumental-resolution-corrected width. Blending and interval
truncation can affect both measurements.

Companion evidence uses only successful current live queries and active catalogue
filters. Missing catalogue queries and failed detection runs are labelled
unavailable. Cross-SPW comparisons require linked apertures and identical known
spectral frames. No intensity ratios are compared, and a missing detection incurs
no penalty: an excitation/opacity model would be needed to predict detectability.
No fixed number of companion lines is required for a suggestion.

Provisional features remain provisional even after a molecule is assigned.
Solutions now also include aperture measurements and an evidence snapshot, while
remaining compatible with earlier version-1 solutions. Loading a solution does
not make its historical evidence current; fetch live matches again. Solutions do
not restore cubes or rerun detection for you; load the same data and settings,
run detection, then load the solution.

The ZIP export includes maps, interactive previews, feature table and provenance:
ESPADA parameter files/logs, aperture settings, live query time/results/filter audit,
adopted transitions, and the reusable identification solution. Edited map intervals do not recompute detection significance.

### Detection and significance

ESPADA runs SoFiA in separate emission and absorption modes. Defaults: 5σ voxel
threshold; unsmoothed and 3/5/9-channel spectral kernels; unsmoothed and one-beam
spatial kernels. Source linking retains the pilot's minimum sizes (5 pixels per
spatial axis, 2 channels, 50 voxels). Thus spatially tiny sources can be missed.
The continuum-subtracted checkbox defaults on; turn it off for SoFiA's linear
continuum subtraction. Detection never depends on the aperture, so it is not rerun when only the aperture changes.

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

### Detection reliability and identification probability

Many ESPADA components are artefacts, and with a ±200 km/s window most of them
still find a catalogue line by chance (about 92% at the IC 860 SPW0
Extragalactic density of ~11 lines/GHz). The app therefore keeps two separate
scores, defined in `conlines/reliability.py`.

**Detection reliability (Step 3)** asks whether a component is real and never uses
the catalogue. Components outside a source region (default: five beams around
the aperture centre) are the artefact control sample. For each component inside,
the app counts same-sign components that are at least as strong in both peak
channel-integrated S/N (over the component's own SoFiA mask) and number of
channels: inside the region, in a near control zone (up to 3× the radius) and in
the far field. Each zone predicts how many such artefacts should fall inside the
region by area; the larger prediction is used, so artefacts concentrated around a
bright nucleus are not diluted by the far field. Reliability is
1 − (84% Poisson upper limit of the prediction) / (observed number). Narrow,
off-centre, isolated artefacts thus lower the reliability of similar components at
the nucleus. Features below the chosen minimum (default 0.5) are not matched or put
in the review queue unless you include them; they stay in the tables and can still
be mapped. The number of neighbouring components within one beam is shown for
context. Assumptions: artefacts are spread evenly within each control zone and
real features outside the region are rare (both make the estimate cautious); it is
not a calibrated false-alarm probability.

**Identification probability (Step 4)** asks which transition a real feature is.
For each catalogue candidate the weight is a Gaussian in velocity offset (mean and
scatter σ set in Step 4, combined with the measured centroid error); an
“unidentified” alternative has the weight of one line spread uniformly over the
match window. Weights are normalised over candidates plus unidentified. Calibrate
the mean and σ from reviewed, significant, reliable anchor lines (at least three)
with **Use anchor calibration**. Combined = detection reliability × identification
probability. Each feature also shows the chance-match probability: the probability
that at least one unrelated catalogue line lies within ±2σ at the current catalogue
density. Species priors, line-strength (LTE) companion predictions and width
consistency are planned next; until then all candidates have equal prior weight.

### Identification methods (Step 4)

Choose one of three methods. Each one only suggests assignments, which stay
tentative until you review them, and detection reliability applies to all three.

- **Velocity likelihood (built-in):** the candidates within the match window,
  ranked by a velocity-offset probability against an “unidentified”
  alternative (see above).
- **ADMIT LineID rules:** a port of the identification logic of ADMIT's
  `LineID_AT` (ALMA Data Mining Toolkit, University of Maryland, MIT licence;
  `conlines/admit_lineid.py`). The original needs CASA, so its rules are applied
  to the ESPADA features:
  1. ADMIT's Tier 1 list of bright, common lines (copied unchanged to
     `conlines/data/admit_tier1_lines.list`) is matched first, within ±40 km/s,
     or ±300 km/s when the systemic velocity exceeds 150 km/s.
  2. Catalogue lines inside each remaining feature are filtered as ADMIT does:
     NRAO-recommended entries with a line strength, no exotic atoms unless
     allowed, and recombination lines limited.
  3. The candidates are narrowed to one: drop rare isotopologues of heavier
     molecules, take the lowest upper-state energy (unless a much lighter
     species has a strong line), and pick the highest line strength. Other
     entries of the same species at the same energy are listed as blends.
  4. A feature with no candidate is unidentified ("U").

  ADMIT's own line finding is not ported, nor its search for equally spaced
  peaks from rotating sources. It uses the live Splatalogue rows of the selected
  astronomical filters, without this app's energy and Einstein A screens; Hot
  Cores plus Extragalactic is closest to ADMIT's full catalogue.
- **Astro AMASE:** runs the real package (Fried et al., arXiv:2608.18221;
  https://github.com/zfried/astro_amase) on the aperture spectrum, converted to K.
  It does LTE simulation of each candidate, checks that its predicted strong
  lines are present and none are predicted too strong, scores chemical relevance
  to the molecules already assigned, and fits column densities at the end.
  - Its line assignments are mapped to the ESPADA features they fall in.
  - A molecule removed by AMASE's best-fit model is not suggested.
  - Absorption features are analysed in a separate run on the sign-flipped
    spectrum. For optically thin lines from one excitation temperature this
    keeps the relative line strengths, but column densities and temperatures
    from that run are not physical.
  - The systemic velocity is passed as v = cz/(1+z), which is the radio
    convention AMASE uses for frequency shifts. You can let AMASE refine it,
    with T_ex, within ±250 km/s.
  - AMASE assumes one T_ex, linewidth and velocity, and optically thin Gaussian
    lines. Expect more false positives in line-confused CON spectra.
  - Its default elements exclude Si; add it for SiO and its isotopologues.

#### Setting up Astro AMASE

AMASE pins different numpy/scipy/astropy versions, so give it its own environment:

```sh
conda create -n amase python=3.11 -y
conda activate amase
git clone https://github.com/zfried/astro_amase
cd astro_amase && pip install -e .
python -c "import sys; print(sys.executable)"   # copy this path
```

Then download the AMASE database folder from the Dropbox link in its README.
The first run also downloads the chemical-embedding model weights.

In Step 4, enter the Python path and the database folder, or set the
`AMASE_PYTHON` and `AMASE_DATA` environment variables. The app remembers the
paths in `.local-tools/amase_paths.json`, which is git-ignored.

Each run takes minutes, so it starts only when you press **Run Astro AMASE**.
The interactive AMASE fit can be downloaded from Step 4, and the run settings
and results go into the export ZIP and the solution provenance.

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

UI tests use synthetic cubes and controlled detector/catalogue results. Tests cover
sequential review, closest-match defaults, unidentified assignments, solution save/load,
context mismatch rejection, map export, zoom independence, stale aperture invalidation,
labels, spectral conversion and baselines. Real ESPADA/SoFiA is checked separately
against the local full-field IC 860 data.
