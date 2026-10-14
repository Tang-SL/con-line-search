# IC 860 Band 5 SPW0 — catalogue-conditioned ESPADA pilot

## Outcome

18 catalogue-defined windows were independently processed by ESPADA in both
absorption and emission. There are 18 accepted window/component associations,
which deduplicate to eight spatial/spectral features. Two pass the conservative
feature robustness screen: absorption at 166.599931 GHz (S/N 9.45) and
167.326913 GHz (S/N 12.89). These are feature detections, not two unique molecule
identifications. See `qa/results.md`, `qa/features.json`, and the PNG plots.

The first feature agrees with H2S 1(1,0)-1(0,1), rest 168.7627624 GHz,
expected 166.606493 GHz; its offset is +11.8 km/s. Neighbouring SO2, 34SO and
Hot Cores transitions remain blend alternatives. The second is near CH3OH
rest 169.488192 GHz (expected 167.322654 GHz; offset -7.6 km/s), but this
transition has Eu=299.59 K and A=4.09e-8 s^-1. Hot Cores also includes 30SiO
4-3 at 169.4868766 GHz (Eu=20.34 K, A=2.42e-4 s^-1; offset about -9.9 km/s).
These two catalogue frequencies differ by only ~2.3 km/s, far less than a
channel. Neither match establishes the identity; abundances, excitation,
optical depth and companion transitions must be tested. A values alone do
not predict observed absorption strength. All alternative matches, including
high-excitation entries and other isotopologues, are retained in features.json.

The previous blind-search S/N for the second absorber was 14.4. Here the
independent cutout run selects a slightly different measurement position/mask;
the same empirical baseline method gives 12.9. Neither value is a global
false-alarm significance. The six remaining components fail the strict screen.

## Frequency and redshift

Aalto et al. 2019, Sect. 3.2, explicitly adopts cz=3880 +/-20 km/s and z=0.01294:
https://arxiv.org/abs/1905.07275
https://pure.manchester.ac.uk/ws/portalfiles/portal/99010941/ic860_ALMA_revised_final.pdf

We use z=3880/299792.458=0.01294229 and nu_obs=nu_rest/(1+z).
This is optical cz, not radio velocity. The original cube's VRAD axis is
converted to frequency using nu=RESTFRQ*(1-v_radio/c). Its frequency frame is
LSRK. The cited systemic-velocity paragraph does not explicitly specify its
reference frame, so the published z is a working value without an unverified
frame correction. This uncertainty must be resolved before interpreting small
velocity offsets as inflow/outflow. The +/-200 km/s exploratory window is much
larger than the quoted +/-20 km/s systemic uncertainty and permits disk motions.

Observed channel centres: 166.023828537–167.868586537 GHz (167 channels).
Systemic rest-frequency interval approximately 168.1726–170.0412 GHz.
The actual catalogue query expands the interval by +/-200 km/s:
168.060363904–170.154629216 GHz. No double redshift correction is applied:
the catalogue receives rest frequencies and frequencyRedshift=0.

## Catalogue selection

Splatalogue: https://splatalogue.online/#/basic
The site's basic filters are radio buttons, so Hot Cores (`ism-hot-core`) and
Extragalactic (`extragalactic`) were queried separately through its own basic
search endpoint. Requests and complete responses, including server metadata,
are cached in `catalogue/`. Responses contain 2412 and 125 actual line entries
respectively (one additional metadata row each). Neither exceeds the server
row cap. Union/deduplication by species id, frequency and quantum numbers gives
2365 entries; 164 are excluded as deuterated/trans. Other rare isotopologues
are retained as alternatives, not assumed abundant. D filtering acts on species
formula/name, not arbitrary text in the chemical description or quantum numbers.

Primary search: 23 non-deuterated, non-trans Extragalactic entries with Eu<=300 K,
grouped into 18 windows by merging same-species entries within one channel.
Members and quantum numbers remain in `catalogue/priors.json`; the group centre
is their unweighted mean, a resolution-level convenience rather than a physical
hyperfine-intensity model. The 300 K cut is an adjustable first-pass screening
choice, not a temperature measurement. Higher excitation transitions and Hot
Cores-only transitions were not exhaustively searched as independent priors.
The full Hot Cores set is used for alternative/blend assessment.

No intensity threshold is imposed: catalogue A and Eu are descriptive priors,
not an abundance/excitation/radiative-transfer model. Some lines are intrinsically
weak. No absence claim or abundance upper limit follows from this experiment.
L01 and L18 systemic centres lie outside the observed band; L01/L02/L17/L18
search intervals are edge-truncated. The results list them for wing coverage,
not as fully observed transitions.

## Injecting priors into ESPADA

ESPADA 1.0.1 has no molecular-catalogue inference stage in this workflow. The
wrapper creates one FITS spectral cutout per prior: +/-200 km/s plus five
channels of padding on each side. It updates CRPIX3 to preserve the original
frequency WCS, retains the common 0.060 arcsec beam, and searches the same
central 1.204 arcsec square. It then runs ordinary ESPADA/SoFiA detection with
use_mask=false and only accepts component centroids inside the prior interval.
This genuinely reruns detection in preselected windows; no source mask is
forced and no line is declared detected merely because a catalogue predicts it.
Original Sabine files are unchanged.

The broad overlapping windows cover 160/167 channels (95.8%). Thus this pilot
adds explicit transition hypotheses and reproducible associations, but provides
little reduction in searched bandwidth. It does not justify lowering detection
thresholds or claiming a large statistical trials advantage. A tighter velocity
prior must be supported independently, not adjusted to make features significant.

SoFiA settings match the prior all-SPW run: S+C threshold 5; spatial kernels
0/15 pixels; spectral kernels 0/3/5 channels; minSizeZ=2, minPixels=50;
no reliability model because real absorption is present. Noise is re-estimated
within each cutout, so the candidate masks need not equal a slice of the blind
run. Overlapping detections are deduplicated within 0.024 arcsec and 3 channels,
keeping the realization furthest from a cutout edge; all triggering priors are
retained. Nearby source substructure may need manual deblending later.

Measurement uses the full SPW cube and the earlier reviewed baseline functions:
constant median with five clean channels per side versus linear fit with ten,
starting beyond a two-channel gap, searching at most 30 channels. Neighbouring
features within 0.09 arcsec are excluded from sidebands. Integrated noise is
the largest std or MAD scale from two off-nuclear annuli, with the identical
baseline operation applied to sky pixels, preserving spectral correlations.
Robust requires both integrated S/N>=8, two adjacent channels>=3 sigma, <=30%
baseline variation, full sidebands and no edge/noisy-channel flags. This is
an operational robustness screen, not a calibrated FDR or molecular identity
probability. Sidebands can still contain undetected lines.

## Files and reproduction

- `catalogue/*_request.json`, `*_raw.json`: public API requests and full responses.
- `catalogue/normalized.json`: both filters and exclusions, Eu, A, quantum numbers.
- `catalogue/priors.json`: primary line groups, predicted positions and channel ranges.
- `config/settings.json`: working source and selection settings.
- `config/espada.yaml`, `abs.par`, `emi.par`: actual detector configuration.
- `qa/results.md`: all primary hypotheses, including no accepted components.
- `qa/features.json`: measurements, ambiguity lists and flags.
- `qa/catalogue_overview.png`, `targeted_windows.png`: inspected scientific plots.
- `qa/validation.json`: FITS/WCS and run completion checks.
- `inputs/`, `run/`: cutouts and independent ESPADA outputs.

See the repository root README for portable setup and reproduction commands.
The pipeline remains IC 860-specific. Changing source/redshift requires matching
catalogue queries, cube preparation and validated frame metadata.

A general successor should separate source/frame metadata, catalogue priors,
detection, joint multi-transition excitation/absorption modelling, and validation
with injection/recovery. It should retain competing blends and report upper
limits when sensitivity and line-width assumptions support them.
