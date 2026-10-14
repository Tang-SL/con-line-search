# Selection rules fixed before the four-SPW search

Scope: all four spectral windows, both signs, central 1.204 arcsec square,
common circular 0.06 arcsec beam. Same SoFiA settings as pilot: 5-sigma
smooth-and-clip with XY kernels 0/15 pixels and Z kernels 0/3/5 channels.
This is not a full-field survey or a completeness-limited catalogue.

Measure each candidate at its rounded SoFiA spatial centroid, summing its
detected contiguous channel extent. Do not move the position to maximise S/N.
This measures brightness per beam, not total source flux.

Two deterministic local baselines:
A: median of the 5 channels on each side, after a 2-channel gap.
B: straight-line least-squares fit to 10 channels on each side, same gap.
Neither sideband set is assumed proven line-free; agreement is a sensitivity
check, not proof of a clean baseline. Require both sides in full for the shortlist.

For each model, perform exactly the same channel sum and baseline operation at
background positions. Estimate integrated-map noise in both annuli 0.25–0.40
and 0.40–0.55 arcsec. Exclude positions within 0.09 arcsec of the candidate
centroid. Use the largest of ordinary standard deviation and 1.4826*MAD in
either annulus (no clipping of potential real off-source features). Subtract
the pooled background median. This propagates observed spectral covariance
and baseline noise without assuming independent channels or pixels.

Rank by the lower sign-corrected integrated S/N of models A and B.

**Strict shortlist** requires all:
- worst-case integrated S/N >= 8;
- at least two adjacent channels with sign-corrected brightness >= 3 times
  empirical channel noise after baseline removal, for both models;
- integral changes by <=30% between models, relative to the larger absolute value;
- full baseline coverage, centroid >=0.06 arcsec from spatial boundaries;
- no channels in its detection interval with noise >3 times SPW median or
  within two channels of the spectral boundaries.

**Follow-up candidate**: worst-case S/N >=5 and the same quality/continuity/
stability checks, but below 8. All others remain in the table with explicit
flags. Borderline candidates are not rejected as physically nonexistent.

No molecular-ID, nuclear-location or prior-science preference enters the rank.
Emission/absorption components of one transition may be counted separately.
Neither input pair nor different smoothing scales constitute independent data.

The threshold of 8 is conservative but operational. No catalogue FDR or formal
search-wide p-values are claimed: spatial correlations, uncertain noise tails,
baseline contamination and selection on the same data prevent that calibration.
The two-channel rule also deliberately excludes unresolved one-channel features.
Broad lines/line-rich sidebands can fail these local-baseline tests even if real.

## Documented QA amendment after the first pass

The first-pass SPW1 emission S1-E7 had both baseline windows contaminated by
adjacent absorption S1-A2, artificially enhancing its contrast. The agreement
test alone cannot identify a bias shared by both models. The first-pass results
are retained as `qa/firstpass_*`.

For ALL candidates, baseline channels now exclude the channel ranges (plus one
channel on each side) of ALL detected components within 0.09 arcsec of the
measurement position. Select the nearest 5 (A) or 10 (B) unexcluded channels on
each side, beginning after the same 2-channel gap and looking no farther than
30 channels beyond the candidate interval. The strict shortlist requires the
full requested counts on both sides. A minimum of three available channels
on each side is required even for a provisional measurement. Otherwise the component is unassessed and flagged
`insufficient-clean-sidebands`, rather than promoted or declared absent.
The intermediate fixed-window exclusion results are retained as `qa/fixed_window_*`.
Thresholds and all other rules remain unchanged. This is a documented quality
control revision, not an independent confirmation or a preregistered hypothesis
test. Undetected lines can still contaminate the remaining sidebands.
