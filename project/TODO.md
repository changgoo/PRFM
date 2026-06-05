# TIGRESS-PHANGS-ILI Project TODO

**Branch:** `project/tigress-phangs-ili`
**Last updated:** 2026-06-05

---

## In Progress / Immediate

- [ ] **Refine `design_fields_selection.png`**
    - [ ] Figure refinement with additional design rules (in progress)
    - [ ] Decide whether to show n=64 sample overlay on top of selection figure
    - [ ] Link finalized figure into `results.tex` (currently commented out)

- [ ] **Resolve `\CGK{}` flags in Paper I**
    - [ ] Verify H_star convention: exponential vs sech² in Vijayakumar et al. 2025
          (matters for z_star mapping in `tigress_mapping.tex` eq. 2)
    - [ ] Confirm final aperture count after quality cuts (~46k joined, ~29k with Sigma_gas>0)
    - [ ] Replace `\CGK{Section~3}` and `\CGK{Section~6}` cross-refs with `\autoref{}`
          once section labels are confirmed

- [ ] **Bibliography cleanup**
    - [ ] Add `Scott1992` entry (Multivariate Density Estimation) for KDE bandwidth citation
    - [ ] Verify or remove `2023ApJS..264...10K` (second TIGRESS-NCR ref flagged in sampling.tex)
    - [ ] Fill `CITATION_NEEDED_TBD` placeholders in discussion.tex and paper2 intro
          (wavelet scattering: Allys+2019, Cheng+2020; SBI: Cranmer+2020 already in bib)

---

## Sampling & Code

- [ ] **Wire directional tail extension into `synthesize_expanded_kde_sobol`**
    - Currently only broadens bandwidth (factor 1.5); the `lower_padding_dex` /
      `upper_padding_dex` config fields (Ω: −0.3 dex, H_star: −0.25 dex,
      q: −0.15 dex) are not yet applied to the quantile function bounds
    - Document in `TIGRESS-PHANGS-sobol-sampling.md` §6

- [ ] **Update `project/prfm_phangs_sampling.ipynb`**
    - Replace LHS-based sampling cells with `synthesize_kde_sobol` calls
    - Add fairness comparison: LHS vs Sobol at same n

- [ ] **Run sampling for all final suite sizes**
    - Pilot: n=32 at Σ_gas=10 (already done → `project/output/`)
    - First suite: n=64, 128 at Σ_gas=5, 10, 20
    - Produce final design tables for TIGRESS handoff

---

## Simulations (TIGRESS-NCR)

- [ ] **Run pilot suite: n=32 at Σ_gas=10**
    - Design CSV: `project/output/design_Sgas10.0_n0032.csv`
    - Confirmed mapping:
        - z_star = H_star = Sigma_star / (4 * rho_star_mid)
        - rho_dm = 0.0064 M_sun/pc^3 (R8 fiducial, fixed)
        - Zg = Zd = 1 (fiducial metallicity)

- [ ] **Validate pilot suite against PHANGS**
    - Compare emergent Sigma_SFR, f_mol, P_DE, sigma_eff distributions
    - Fill in `results.tex` §5.2 (Simulated ISM Properties) and §5.3 (Validation)

---

## Paper I (`project/doc/paper1-suite/`)

- [ ] **Finalize results.tex figures**
    - [ ] `fig:field_distributions`: link `design_fields_distributions.png` once finalized
    - [ ] `fig:design_fields`: link `design_fields_corner_*.png` (caption already written)
    - [ ] Add expanded-prior comparison figure

- [ ] **Complete abstract** (currently placeholder)

- [ ] **Add authors** (currently only Chang-Goo Kim)

- [ ] **§4 TIGRESS mapping**: fill in z_star ↔ H_star conversion note once H_star
      convention (Vijayakumar+2025 sech² vs exponential) is confirmed

- [ ] **§5 Results**: fill in once pilot simulations are done

- [ ] **§6 Discussion**: review and expand `CITATION_NEEDED_TBD` references

---

## Paper II (`project/doc/paper2-ili/`)

- [ ] Design ILI/NPE architecture (neural posterior estimator)
- [ ] Build training pipeline (simulation → summary statistics → NPE)
- [ ] Run inference on mock data (posterior recovery check)
- [ ] Fill in `inference.tex` and `results.tex` placeholders
- [ ] Complete abstract once framework is implemented

---

## Completed (this session, 2026-06-05)

- [x] **KDE-Sobol sampling implementation** (`prfm/phangs_sampling.py`)
    - `_compute_marginal_quantiles()`, `synthesize_kde_sobol()`,
      `synthesize_expanded_kde_sobol()`, updated `sample()` dispatch
    - Nesting property S(2^k) ⊂ S(2^{k+1}) verified
    - Unit tests (67 pass)

- [x] **Sobol sampling design document** (`project/TIGRESS-PHANGS-sobol-sampling.md`)

- [x] **Two AASTeX7 paper drafts** (`project/doc/`)
    - Paper I: sampling design + PHANGS data + TIGRESS mapping + results stub
    - Paper II: ILI skeleton
    - Shared: `preamble.sty`, `references.bib` (adstex), `aastex702.cls`,
      `aasjournalv7.bst`
    - Both compile cleanly

- [x] **Figure generation scripts** (`project/scripts/`)
    - `run_sampling.py`: design CSVs + fairness summary
    - `plot_sampling_diagnostics.py`: distribution overlays, pair overlays
    - `plot_design_fields.py`: reduced correlation matrix (selection figure),
      corner plots, standard vs expanded distributions
    - `prfm.mplstyle`: CM font, inward ticks, consistent styling

- [x] **Design output** (`project/output/`)
    - n=32, 64, 128, 256 CSVs for Σ_gas=5, 10, 20
    - Fairness ε < 0.07 dex across all fields and targets

- [x] **Paper I results.tex §5.1**: fairness table + reference pixel counts

- [x] **Paper I tigress_mapping.tex**: confirmed z_star=H_star, ρ_dm=0.0064 M☉/pc³
