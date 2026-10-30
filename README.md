# EAGW2026
Tutorials and Hands-On Sessio

# AGN × Gravitational Waves — Working Group

Exploratory project on the association between binary black hole mergers detected in
gravitational waves (GW) and flares from active galactic nuclei (AGN). The working
hypothesis is the AGN-disk formation channel: a fraction of stellar-mass BBH mergers
occur embedded in AGN accretion disks, where the merger remnant's recoil through dense
gas can produce a detectable electromagnetic flare.

Organised around the [EAGW 2026 workshop](https://sites.google.com/view/eagw2026/home).
**All participants should register.**

---

## Objectives

1. Pick one concrete, scoped project rather than several loose threads.
2. Work out the full toolchain required for that project end to end.
3. Flesh out the primary project before adding secondary components.

Working format: 5–7 teams of pairs.

---

## Primary Project — ZTF Flare / GW Association

Cross-match GW events against the ZTF AGN flare catalogue to test whether any flares are
plausibly associated with mergers, and to quantify the significance of such associations.

### Proposed workflow

1. **Ingest GW events** from GWOSC, including sky localisation maps and luminosity
   distance posteriors.
2. **Ingest AGN flares** from the ZTF AGN-flare catalogue; retrieve original light curves
   from the ZTF database as needed.
3. **Spatial and temporal cross-match.** For each GW event, search for flares appearing
   near the same sky position within a short time window after the event.
4. **Score candidate matches.** An association is stronger when the flare lies closer to
   the GW sky localisation, occurs sooner after the event, and has a host-galaxy distance
   consistent with the GW distance estimate.
5. **Null test.** Shuffle flare times and repeat the analysis to establish the rate of
   chance coincidences.
6. **Case study.** Apply the pipeline to GW190521 and the candidate flare ZTF19abanrhr as
   a validation benchmark.

### Known challenges

- **Limited multi-band information.** Per Albert: association is difficult because we have
  little information across different filters, and optical data alone lacks the spatial
  resolution to discriminate. An optical light curve on its own may not be sufficient.
- Possible mitigation: classify by AGN type (Type I, Type II, blazar) to see whether
  viewing geometry affects detectability. *Open question: would an edge-on AGN be
  observable at all?*

---

## Secondary and Exploratory Ideas

- **Quick GW parameter estimation** — lightweight tool for rapid parameter estimates.
- **AGN flare simulation** — forward-model expected flare signatures.
- **Observational scenario construction** — build synthetic flare + GW scenarios and
  compare against real data.
- **Binary orientation vs. formation channel** — test whether spin/orbital orientation
  correlates with an AGN origin.
- **EM signature of (micro)lensing in an AGN** (Shigeo) — electromagnetic analogue of
  lensed GW signals from AGN disks.
- **Feasibility check needed:** are full AGN radiative-transfer simulations too
  computationally expensive for this scope?

---

## Data Sources

| Resource | Use | Link |
|---|---|---|
| GWOSC | GW event catalogue, skymaps, posteriors | https://gwosc.org |
| ZTF AGN flare catalogue | Flare candidates (largely public) | https://github.com/Lyle0831/AGN-Flares |
| ZTF Lightcurve API (IRSA) | Original light curves | https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html |

---

## Software and Tools

| Tool | Purpose | Link |
|---|---|---|
| `gwfast` | Fast GW parameter estimation / Fisher forecasting | https://github.com/CosmoStatGW/gwfast |
| `bbhagn` | BBH–AGN association analysis (companion code to Cabrera et al.) | https://github.com/tomas-cabrera/bbhagn |
| SKIRT | Radiative transfer for AGN/dusty media | https://skirt.ugent.be |
| ASCL 2103.01 | AGN-related simulation code | https://ascl.net/2103.01 |

---

## References

### AGN formation channel — foundational

1. Bartos, Kocsis, Haiman & Márka (2017), *ApJ* **835**, 165 —
   [DOI](https://doi.org/10.3847/1538-4357/835/2/165)
2. Bartos, Haiman, Márka et al. (2017), *Nature Communications* **8**, 831 —
   [DOI](https://doi.org/10.1038/s41467-017-00851-7)
3. *ApJ* (2018) — [DOI](https://doi.org/10.3847/1538-4357/aadae5)
   <!-- TODO: confirm citation; likely McKernan et al. (2018), ApJ 866, 66 -->
4. *ApJL* (2021) — [DOI](https://doi.org/10.3847/2041-8213/ac2cc1)
   <!-- TODO: confirm full citation -->

### Electromagnetic counterparts

5. Graham et al. (2020), *PRL* **124**, 251102 — *Candidate Electromagnetic Counterpart to
   the Binary Black Hole Merger Gravitational-Wave Event S190521g* —
   [DOI](https://doi.org/10.1103/PhysRevLett.124.251102) ·
   [arXiv:2006.14122](https://arxiv.org/abs/2006.14122)
6. *A Standard Siren Cosmological Measurement from the Potential GW190521 Electromagnetic
   Counterpart ZTF19abanrhr* — [arXiv:2009.14057](https://arxiv.org/abs/2009.14057)

### Lensing constraints

7. *Constraining Binary Mergers in AGN Disks Using the Non-observation of Lensed
   Gravitational Waves* — [arXiv:2408.13144](https://arxiv.org/abs/2408.13144) ·
   published as *ApJL*, [DOI](https://doi.org/10.3847/2041-8213/ad9ead)

### Closely related to this project

8. *ApJL* (2025) — [DOI](https://doi.org/10.3847/2041-8213/adf31f) — **the most directly
   relevant prior work.** Citing papers are worth tracking:
   [SciX citation list](https://scixplorer.org/abs/2025ApJ...989L..15Z/citations?p=1)
9. *ApJ* — [DOI](https://doi.org/10.3847/1538-4357/ae4d37) — companion analysis;
   code at [`tomas-cabrera/bbhagn`](https://github.com/tomas-cabrera/bbhagn)

---

## Meeting Log

### 2026-08-24
- Group structure set: 5–7 teams of pairs.
- Initial AGN reference list circulated, covering lensing optical depth, the S190521g
  candidate counterpart, and the GW190521 standard-siren measurement.

### 2026-10-09
- Workshop registration opened — all members to register.
- Agreed to converge on a single concrete project and build out its complete toolchain.
- ZTF flare association identified as the leading candidate; catalogue and light-curve
  access confirmed as largely public.
- Surojit proposed the GWOSC × ZTF cross-match workflow with time-shuffling null test
  (see Primary Project above).
- Noted strong overlap with recent literature; need to position our work relative to it.

### 2026-10-23
- Project discussion and introductions.
- **Action:** Albert to give a tutorial covering ZTF flare + AGN association (with and
  without timing information) and GW + AGN association.
- Confirmed an existing ZTF AGN catalogue is available.
- Shigeo proposed an additional direction: EM signatures of (micro)lensing within an AGN,
  paralleling the lensed-GW analysis.

---

## Open Action Items

- [ ] All members register for EAGW 2026
- [ ] Albert — deliver ZTF flare + GW + AGN association tutorial
- [ ] Confirm citations for references 3, 4, 8 and 9
- [ ] Review the citation list of the closely related ApJL paper to identify gaps our
      project can fill
- [ ] Assess computational cost of AGN simulation tools before committing to that branch
- [ ] Scope and freeze Project 1 before adding secondary componentsn
