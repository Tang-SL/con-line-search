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

Double-click **CON Line Search.app** in the project folder; you can also drag it
to the Dock. A Terminal window opens and starts the app in the background, then
your browser opens once the app is ready. You can close that window; the app keeps
running. A copy in /Applications also works once the app has been started from
the project folder.

Starting goes through Terminal because macOS privacy protection stops a
home-made, unsigned app from reading files in Documents, while Terminal already
has that access. The icon itself only checks, opens and stops the app, so it
needs no folder access.

- If the app is already running, the icon asks whether to **Open**, **Restart** or **Stop** it.
  Restart after pulling or switching branches: Streamlit keeps the helper modules it
  imported at start, so an old server would mix old and new code. `Start app.command`
  and `scripts/launch.sh start` detect changed code and restart automatically.
- **Start app.command** does the same in a Terminal window and shows progress.
  **Stop app.command** stops the app.
- On the first start, a missing `.venv` is created with Python 3.12 or newer, and
  `requirements-ui.txt` is installed (a few minutes, needs internet).
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

Open http://127.0.0.1:8501. Cubes stay local.

### Using the app: five steps

The app is a guided sequence. Only the current step is open, with a short
**What to do** line and one main button. A finished step folds into a single line
with a ✓ and a summary (for example *2 · Find lines — 8 features found across
1 cube*); click it to look at or change it. Later steps follow the change: new
cubes or detection settings send you back to Step 2, while moving the aperture or
changing the redshift just updates the later steps. Each setting appears in the
step that uses it, and expert settings sit in folded **(advanced)** sections with
sensible defaults. **Start over** clears everything.

1. **Load data.** Upload FITS cubes, give local file paths (better for large
   cubes), or pick the IC 860 full-field example; one cube per spectral window.
   A table shows each cube's channels, frequency range, beam and frame. Press
   **Continue**.
2. **Find lines.** Press **Find lines in all cubes**. ESPADA/SoFiA searches each
   full cube for emission and absorption, independent of any aperture, redshift or
   catalogue. Every feature is shown on a map and in a table, with an ID that
   belongs to the cube's detection (for example `SPW0-F003`) and stays fixed when
   the aperture moves. *Advanced:* voxel threshold, trial widths, emission/absorption,
   continuum subtraction. Detection reruns only when the cubes or these settings change.
3. **Choose region.** Choose where to extract a spectrum: image centre, brightness
   peak, the continuum nucleus (example), a detected feature, or a pixel position;
   then the aperture shape and radius (default one beam). The spectrum and the
   features at that position appear at once, with counts of features, significant
   ones, reliable ones and likely artefacts set aside (see *Detection reliability*
   below). *Advanced:* linking the sky position across SPWs, preview zoom, the
   significance threshold, the source-region radius and the reliability cut.
   Moving the aperture re-screens without rerunning ESPADA. Press **Continue to
   line identification**.
4. **Identify lines.** Set the redshift (IC 860 example: cz = 3880 km/s), then
   pick a method (see *Identification methods* below):
   - *Velocity likelihood* and *ADMIT* need a Splatalogue search: choose the
     astronomical filters (and, for the built-in method, the match window), then
     press **Search Splatalogue**. Catalogue filters (energy, Einstein A,
     isotopologues, excited states) are under *advanced*; ADMIT applies its own
     rules, so it uses only the filters.
   - *Astro AMASE* uses its own CDMS/JPL database, not Splatalogue. If it is not
     installed, the step shows **Set up Astro AMASE…**. Then press **Run Astro AMASE**.
     A Splatalogue search is optional and only feeds the companion-line evidence.

   The suggestions are tentative. Press **Review identifications one by one**,
   choose a transition (or *Unidentified*) for each feature with the help of the
   evidence panel, and **Confirm & next**. Reviewed choices are kept when you
   search again, switch method or move the aperture. *Save or load an
   identification solution* (folded) stores assignments, review status, channel
   ranges, data fingerprints, settings and query provenance as JSON; a solution
   loads only onto the same data and analysis. *Advanced:* which features to
   identify (provisional, low-reliability) and the companion velocity tolerance.
   Press **Continue to maps**.
5. **Make maps.** Choose the features to map (default: those you identified and
   reviewed), check their channel ranges (editable; this changes the maps only,
   not the detection), choose the baseline subtraction and press **Generate maps &
   download bundle**. Unidentified features can be mapped too: they use the
   observed centroid as the velocity-width reference (VREFGHZ) instead of an
   assumed rest frequency. IDSTAT records TENTATIVE or REVIEWED. Absorption keeps
   its negative sign.

Details on matching and evidence:
the closest velocity match is proposed within the default ±200 km/s window.
Queries cover each loaded SPW's full frequency range, including SPWs with no
detections, and every search fetches fresh live data; a network failure leaves
detections intact without substituting a saved catalogue. Archived catalogue
snapshots are used only to reproduce the historical pilot. The **Evidence for this
identification** panel shows the aperture centroid and its error, the signed
second-moment width, catalogue metadata and other transitions of the same species
entry. Companion matches must lie inside another feature's interval and agree in
velocity within the default 50 km/s; significant and provisional features are
counted separately, and these are frequency-compatibility checks, not probabilities.

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

In Step 4, choose Astro AMASE and press **Set up Astro AMASE…**, or double-click
**Set up Astro AMASE.command** in the project folder. It installs everything into
`.local-tools/`, which is git-ignored:

- the Astro AMASE source;
- a private Python 3.11 environment (Astro AMASE cannot use 3.12: its
  `astrochem_embedding` dependency needs `rdkit-pypi`, which has no 3.12 builds);
- its database, from the Dropbox link in the AMASE README.

It uses uv, which supplies Python 3.11 with Tk support. Tk is needed because the
`group_selfies` dependency imports `tkinter`.

The downloads (PyTorch, RDKit, the database) are large and take several minutes.
When Terminal says “Done”, press **Use installed Astro AMASE** and the paths fill in.
Re-running the setup skips finished steps and updates the source.

To use your own installation instead, enter its Python and database folder, or set
`AMASE_PYTHON` and `AMASE_DATA`. The folder must contain
`all_cdms_final_official.csv` and `cdms_pkl/`.

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
