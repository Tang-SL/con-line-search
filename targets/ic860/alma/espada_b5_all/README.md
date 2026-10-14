# IC 860 — all four Band 5 SPWs: detection and ranking

All four Sabine spectral windows were searched in emission and absorption.
Spatial scope is the central **1.204 arcsec square**, not the entire original
4.096 arcsec field. The input cubes were convolved channel-by-channel to a
common **0.060 arcsec circular beam**, retaining Jy/beam through beam-area
scaling. Original FITS files are unchanged. No new primary beam correction
was possible; its status remains unknown.

## Results

| SPW | Absorption components | Emission components |
|---|---:|---:|
| 0 | 2 | 5 |
| 1 | 2 | 9 |
| 2 | 2 | 5 |
| 3 | 0 | 5 |

There are **30 candidate components**, not 30 identified transitions.
Emission and absorption at different positions/velocities can belong to the
same molecular transition. Molecular identifications were not used in ranking.

The strict shortlist contains three nuclear absorption components:

| Candidate | SPW | Observed LSRK centroid (GHz) | Lower integrated S/N | Baseline change |
|---|---:|---:|---:|---:|
| S0-A2 | 0 | 167.326679 | 14.4 | 8% |
| S1-A1 | 1 | 168.143224 | 12.7 | 6% |
| S0-A1 | 0 | 166.599973 | 9.5 | 1% |

S0-A1 is associated with the H2S-labelled transition. The other molecular
identifications remain unassigned in this analysis.

**S1-A2 at 169.324203 GHz deserves priority follow-up:** its lower provisional
S/N is 14.7, but the wider baseline model has only 16 usable channels rather
than the required 20 (10 on each side). It fails the full baseline-coverage
rule. This is a quality flag, not evidence that the absorption is spurious.

No component falls in the clean 5–8 S/N follow-up class. Five components cannot
be scored by both models because there are fewer than three uncontaminated
channels on at least one side: S2-E1, S2-E5, S3-E1, S3-E4, S3-E5. They are
listed as N/A, not as zero-significance detections. In particular, **the water
candidate S3-E1 requires a dedicated baseline/blend assessment**. Its earlier
pilot estimate was based on chosen sidebands and does not meet the uniform
all-SPW sideband requirements.

## Selection criteria

Exact rules and the documented QA amendment are in
`config/selection_rules.md`. The strict shortlist requires:

1. Integrated signal/noise >=8 under **both** local baseline models.
2. At least two adjacent channels with same-sign contrast >=3 times the
   empirical channel noise under both models.
3. Integral changes by no more than 30% between models.
4. Full uncontaminated baseline coverage and no specified spatial/spectral
   boundary or high-noise-channel flags.

The lower S/N sets the rank; molecule labels and distance from the nucleus do
not affect it. Strong but flagged components remain visible in the full table.
This ranks **measurement robustness**, not astrophysical importance or novelty.

The candidate catalogue comes from the same ESPADA/SoFiA smooth-and-clip settings
as the pilot (threshold 5, spatial kernels 0/15 pixels, spectral kernels 0/3/5
channels). The ranking then measures each feature at its rounded catalogue
centroid, without moving that position to maximise the signal. Integrals sum
the candidate's entire contiguous channel extent and are per-beam brightness,
not total source fluxes.

Baseline A is a constant median of the nearest five usable channels on each
side; baseline B is a straight-line fit to ten on each side. Both start after
a two-channel gap and search at most 30 channels outward. Known components
within 0.09 arcsec of the measurement position, with one-channel spectral guards,
are excluded. Undetected lines can still contaminate these sidebands.

Noise is measured by applying the identical baseline and channel-sum operation
to surrounding sky. We take the largest of standard deviation and 1.4826*MAD
in either of two annuli (0.25–0.40 and 0.40–0.55 arcsec), excluding a 0.09-arcsec
circle around the candidate. No outlier clipping is used. This retains actual
channel correlations and baseline noise; correlated pixels are not counted
as independent samples to improve significance. Real emission in the control
annuli can inflate this estimate and make the selection conservative.

**These S/N values are not calibrated search-wide probabilities.** The 8 threshold
is an operational conservative screen following a multi-window blind search.
We have not measured a false-discovery rate or completeness with injection/
recovery or independent observations. Broad weak lines, one-channel lines,
extended low-surface-brightness emission and crowded spectra can be excluded
by these rules even if real. A failure is not a physical non-detection.

## Why the baseline check changed the first-pass ranking

The first pass compared fixed median and linear sidebands. SPW1 emission S1-E7
appeared to pass at S/N 9.6, but visual inspection showed the sidebands included
nearby absorption S1-A2. Both models shared the bias. Excluding neighbouring
detected features reduces its lower S/N to about 4.3. Its baseline coverage and
channel support also fail the rules. Agreement between two contaminated models
was therefore insufficient.

The exclusion policy was applied uniformly to all candidates, with thresholds
unchanged. `qa/firstpass_*` and `qa/fixed_window_*` preserve earlier results.
This is a transparently documented QA revision, not a preregistered statistical
test. The current authoritative products have no such prefix.

## Files

* `qa/ranked_candidates.md`: complete ranked table, including flags and N/A.
* `qa/ranked_candidates.json`: numerical measurements, sidebands, model results
  and matched checks in the separately imaged +cont products. The +cont data
  are not independent observations and are not used to claim confirmation.
* `qa/ranking.png`: all-component ranking, with class colours.
* `qa/shortlist_1.png`: spectra and signed integrated maps for the strict shortlist.
* `qa/preflight.png`, `qa/preflight.json`, `qa/spectra.npz`: preparation QA.
* `inputs/`: single-HDU common-beam cubes, with and without +cont.
* `run/espada_spw*/`: source catalogues, masks, moments and cubelets.
* `config/pilot.yaml`, `config/abs.par`, `config/emi.par`: actual search settings.

SPW0/SPW3 pilot outputs remain separately in `../espada_b5_pilot/`.

## Reproduction and checks

Use `scripts/ic860_espada_b5_all_prepare.py`, then ESPADA with
`config/pilot.yaml`, then `scripts/ic860_espada_b5_rank.py` from the CONs project.
See the repository root README for environment setup. Put SoFiA on PATH and
set MPLBACKEND=Agg for headless plotting.

ESPADA 1.0.1, SoFiA 2.7.1 (commit
`a7434d9f91605cd738f2b13cd31b1450304fbc11`), WCSLIB 8.2.2; the same Python
environment recorded in the repository requirements and root README.
All eight prepared cubes passed finite-pixel, single-HDU, common-beam and FITS
checksum checks (`qa/input_checks.json`). The convolution method was previously
checked for point-source peak and integrated-flux preservation in the pilot.

The search succeeded for all pairs with detections. SPW3 absorption returns
SoFiA code 8 (no candidates), not a processing failure. Automatic ESPADA HTML
thumbnails still inherit the SoFiA EPS `nan` issue described in the pilot; use
the independently generated and visually inspected PNGs above. No candidate
selection relies on those thumbnails.

SoFiA's positive-versus-negative reliability model was not used because real
absorption is present; see its [official tutorial](https://gitlab.com/SoFiA-Admin/SoFiA-2/-/wikis/SoFiA-Tutorial).
