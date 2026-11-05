# Playbook: identifying molecular lines in galaxies, and in compact obscured nuclei

Distilled from the literature. Cite these papers in the weblog when a rule is applied. Before
relying on any number here for a specific source, check the source's own literature.

## 1. How line identification is done in extragalactic spectral surveys

- **Fit species, not single lines.** In line-rich sources with broad lines, fit a synthetic
  (LTE) spectrum of each species, add species one by one from bright to faint, and subtract
  blends before fitting the next species. Gaussian fits to individual features are
  unreliable when lines overlap.
  - Arp 220 (SMA 1.3 mm survey): ~1.8 lines per GHz; the band is partly confusion-limited;
    line widths of 250–420 km/s. Velocity and width were fixed to CO, and only N and T_ex
    were fitted per species (Martín et al. 2011, A&A 527, A36; arXiv:1012.3753).
  - NGC 4418 (ALMA spectral scan): 317 lines above 3σ from 45 species. A synthetic LTE and
    non-LTE fit used one 120 km/s line width from blend-free lines, fitting molecules from
    brightest to faintest; blending adds ~20% uncertainty (Costagliola et al. 2015, A&A 582,
    A91; arXiv:1506.09027).
  - NGC 253 (ALCHEMI): MADCUBA/SLIM LTE fits, with 146 species and 42,121 transitions in the
    model and ~6–8 lines per GHz at 125–375 GHz. When unconstrained, the velocity and width
    were fixed, T_ex was set to a typical value, and isotopologues inherited their parent's
    parameters (Martín et al. 2021, A&A 656, A46; arXiv:2109.08638).
- **Detection criteria. A species counts as detected only if:**
  - Its transitions are consistent with one excitation model.
  - Every transition the model predicts to be detectable and unblended (ALCHEMI: above 5σ)
    is actually detected. No predicted line may be missing (Snyder et al. 2005, ApJ 619,
    914; Belloche et al. 2008/2013).
  - The rest frequencies are accurate and agree with the source velocity.
  - The relative intensities are consistent.
- **Tentative.** A species seen in only one transition (Costagliola et al. 2015), unless it is
  an isotopologue of a secure parent. Also any case where the fit is unsatisfactory.
- **Unidentified features are normal.** 40–50% of lines remain unidentified in Sgr B2
  surveys. Most likely they are excited states of known molecules, or species missing from
  the catalogues. Report them; don't force an identification.
- **Known misidentification:**
  - H31α in Arp 220 was really vibrationally excited HC₃N (Martín et al. 2011).
  - In broad-line sources, recombination-line claims need unbiased coverage.

## 2. Is a feature real? (significance and fidelity)

- **Look-elsewhere effect.** A search over N_beams × N_channels independent elements finds
  noise peaks at any S/N.
  - ASPECS measures fidelity from negative detections: fidelity = 1 − N_neg/N_pos at the same
    S/N and kernel width. It keeps candidates with fidelity ≥ 0.9 (González-López et al.
    2019, ApJ 882, 139; arXiv:1903.09161).
  - Broader kernels mean fewer independent elements, so a broad line at a given S/N is more
    significant than a narrow one.
- **In this toolkit:**
  - *significance* = integrated S/N in the aperture under two baselines, ≥2 adjacent channels
    at 3σ, and ≤30% baseline sensitivity.
  - *reliability* = 1 − (artefacts expected from equally strong and wide components away from
    the source)/(observed). This is the same idea as fidelity, using spatial controls; it is
    not a calibrated false-alarm probability.
- **Checks to apply by eye (plots):**
  - The feature is spatially coincident with the source.
  - It spans more than one or two channels, or is consistent with the line width.
  - It is not at a band edge.
  - It doesn't coincide with an atmospheric line or a noisy channel.
  - It is present under either baseline.
  - Its sign makes sense: absorption needs continuum behind it.
- **Interpreting the S/N:**
  - A 3–5σ feature near a bright source can be a sidelobe or a cleaning artefact. Check
    whether it shares the bright line's peak channels.
  - A single-channel spike is never a line in data with Hanning-smoothed channels.

## 3. Is it physically possible? (plausibility)

- **Companion lines (the strongest test).**
  - Under LTE and optically thin emission, integrated flux ∝ ν³ S µ² exp(−E_u/T). The other
    transitions of the species in the observed windows must then appear with the predicted
    ratios at some reasonable T.
  - The `identify` step turns this into f_max(T): the largest fraction of the feature the
    species can explain without a missing line.
  - A temperature at which no companion is testable proves nothing.
  - A companion that falls on another feature is "blended": no support, but it can still be
    missing.
- **Optical depth and self-absorption.** In CONs, HCN, HCO⁺ and CO are often optically thick
  and self-absorbed. Their companions (isotopologues, other J) are then much weaker than
  thin-LTE scaling predicts, so treat companion tests for these species loosely and say so.
- **Abundance prior.** Strong features in galaxies come from the abundant species:
  - Bright: CO and isotopologues, HCN, HCO⁺, HNC, CS, CN, CCH, H₂O, CH₃OH, H₂CO, SiO, N₂H⁺.
  - Common but weaker: SO, SO₂, OCS, HC₃N, CH₃CN, HNCO, H₂S, c-C₃H₂, CH₃CCH, NS, HOC⁺, CH₂NH,
    and isotopologues.
  - An exotic molecule claimed for the strongest feature in a galaxy needs extraordinary
    evidence.
  - Check whether the species has ever been detected outside the Galaxy (CDMS "molecules in
    space" list; McGuire's census of interstellar and extragalactic molecules,
    arXiv:2109.13848).
- **Excitation.** E_u must be reachable:
  - Cold gas (T ~ 10–50 K) doesn't excite E_u ~ 500 K lines.
  - Vibrationally excited lines (E ~ 1000 K) need warm dust and IR pumping (see §4).
- **Kinematics.** One source has one systemic velocity and similar line widths. A candidate
  needing an offset of hundreds of km/s, or a width very different from secure lines, is
  suspect. Exceptions are outflows, absorption against the continuum, and self-absorption.
- **Redshift and frame.**
  - Velocities in the literature may be optical or radio, heliocentric or LSR. At z ≈ 0.04,
    cz(optical) − cz(radio) ≈ 520 km/s.
  - The toolkit uses f_obs = f_rest/(1+z), with offsets in the radio convention relative to
    the predicted frequency.
  - ALMA cubes are LSRK.
  - Always confirm z with a secure bright line (`zfrom`) before trusting offsets.

## 4. Compact obscured nuclei (CONs)

- **Definition.**
  - Σ(HCN-vib J=3–2) > 1 L☉ pc⁻² over a region of radius ≥ 5 pc (CON-quest: Falstad et al.
    2021, A&A 649, A105; arXiv:2102.13563).
  - The earlier criterion was L(HCN-vib)/L_IR > 10⁻⁸.
  - Found in 38% of ULIRGs and 21% of LIRGs, and none of the sub-LIRGs (CON-quest sample).
  - CONs: IRAS 17208-0014, Arp 220 W, IRAS F22491-1808, IRAS 17578-0400, Zw 049.057,
    NGC 4418, ESO 320-G030. Also with HCN-vib: IC 860, NGC 1377, IRAS 20551-4250, Mrk 231
    (Aalto et al. 2015, A&A 584, A42; arXiv:1504.06824).
  - CONs tend to have low IRAS f25/f60.
- **HCN-vib (ν₂=1) needs:** T_dust ≳ 100 K, τ(14 µm) ≳ 5 and N(H₂) ≳ 10²⁴ cm⁻². It is IR
  pumped: E ≈ 1024 K, critical density > 10¹⁰ cm⁻³.
- **Key transitions at J=3–2, with blend traps** (rest frequencies from CDMS; always re-check
  in the catalogue):
  - HCN v=0 at 265.8864 GHz.
  - HCN v₂=1 l=1e at 265.8527 GHz: 34 MHz (≈38 km/s) below HCN 3–2, so inseparable from it
    in broad lines. Constrain it through l=1f, since l=1e ≈ l=1f in thermal conditions.
    Their Einstein A values differ by only 1.5%.
  - HCN v₂=1 l=1f at 267.1993 GHz: ≈ +400 km/s from HCO⁺ 3–2 at 267.5576 GHz, separable in
    LIRGs but partly blended in ULIRGs.
  - HOC⁺ near HCO⁺: a check on HCO⁺ absorption.
  - CH₂NH on the blue side of HCN 3–2 (IC 860, Zw 049.057).
  - CH₃OH near HCN-vib.
  - H¹³CCN between HCN and HCN-vib.
  - At J=4–3: HC₃N J=39–38 overlaps HCN 4–3 (Arp 220 W, IRAS 17208-0014).
  - At 369 GHz: ortho- and para-H₂S are 26 MHz apart.
- **Self-absorption.**
  - HCN and HCO⁺ 3–2 are often double-peaked, with a minimum near systemic.
  - The absorption is blueshifted in Arp 220 W, IRAS 17208-0014 and Zw 049.057 (outflow), and
    redshifted in IC 860 and ESO 320-G030 (inflow).
  - Self-absorbed lines make centroids and widths unreliable.
  - HCN-vib is not self-absorbed: a shared double-peaked shape is a compromise. Test single
    and double shapes and report the difference as a systematic.
- **Other CON lines:** vibrationally excited HC₃N (v₇=1, v₇=2, v₆=1) and HNC v₂=1, all
    detected in NGC 4418, Arp 220 and others. Lines in absorption against the hot continuum
    (CO 6–5, CS 14–13 in IRAS 17208-0014: Baba et al. 2022, arXiv:2202.11105).
- **Continuum.** The 1.1 mm continuum surface brightness tracks Σ(HCN-vib) closely. Pipeline
  continuum subtraction can remove part of the broad line emission in line-rich CONs, so check
  the baseline beside broad lines.

## 5. Ambiguous and faint cases

- **Write down every plausible solution** (hypothesis), including "unidentified" and "no line"
  for weak features.
- **Rank the solutions by:**
  1. *Robustness to noise*: the `robust` win fraction over noise realisations from the
     spatial controls.
  2. *Parsimony*: lower BIC; fewer free lines.
  3. *Physical expectation*: abundance prior, excitation, literature detections in this or
     similar sources.
  4. *Consistency across windows*: companions, ties.
- **When the fit is model-limited** (reduced χ² ≳ 2), noise robustness overstates certainty.
  Compare line shapes (one or two Gaussians, a template) and report the spread as a
  systematic error.
- **Baseline and fit window.** For lines hundreds of km/s wide, a baseline fitted only on a
  few channels beside the feature is degenerate with the broad components. A faint blended
  line then absorbs the difference. Example from the toolkit's IRAS 17208-0014 test: HCN-vib
  l=1f changed by 40% between 10 and 80 line-free channels per side. Fit with enough line-free
  channels (default: a third of the feature width on each side), check the fitted baseline
  against line-free channels measured directly, and report the change as a systematic.
- **Prefer a weaker but robust claim**, e.g. "HCN-vib l=1f is clearly needed; the l=1e/l=1f
  ratio depends on the line-shape model", over a precise but fragile one.
