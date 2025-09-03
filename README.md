## RedundancyNWT
This repository contains scripts and outputs for the functional resilience and turnover analysis of boreal landbird communities under projected climate and habitat change in the Taiga Plains in the Northwest Territories in Canada. Proposed in the manuscript entitled:Northern landbird communities show geographically constrained resilience to projected climate and forest changes.
The workflow integrates functional group (FG) binary matrices, temporal beta diversity (TBI), body-mass discontinuity analyses, redundancy metrics and Non-Metric Multidimensional Scaling (NMDS) analyses.

📂 Repository Structure
inputs/    # Raw or input data
tables/    # Processed tabular outputs (e.g., FG binary matrices, log mass matrices)
outputs/   # Intermediate outputs (TBI results, discontinuities, redundancy metrics, NMDS)
scripts/   # All R scripts used in the workflow
scripts/
 ├── 01_fg_tbi/
 │    ├── 01_build_fg_binary_matrices.R
 │    ├── 02_calculate_fg_tbi.R
 │    ├── 03_generate_fg_tbi_rasters.R
 │    ├── 04_plot_fg_tbi_changes.R
 │    └── 05_summarize_fg_tbi_changes.R
 ├── 02_bodymass/
 │    ├── 06_build_logmass_matrices.R
 │    ├── 07_calculate_bodymass_discontinuities.R
 │    └── 08_plot_bodymass_aggregations.R
 ├── 03_resilience/
 │    ├── 09_calculate_resilience_metrics.R
 │    ├── 09a_within_scale_redundancy.R
 │    ├── 10_plot_resilience_metrics.R
 │    └── 11_plot_baseline_resilience_metrics.R
 ├── 04_nmds/
 │    └── 12_run_nmds_5regions.R
 └── 05_richness/
      ├── 13_plot_fg_richness_maps.R
      └── 14_plot_delta_richness_maps.R

flowchart TD
  A[01 Build FG matrices] --> B[02 Calculate FG TBI]
  B --> C[03 Generate TBI rasters]
  C --> D[04 Plot FG TBI changes]
  D --> E[05 Summarize FG TBI changes]
  E --> F[06 Build log-mass matrices]
  F --> G[07 Calculate body-mass discontinuities]
  G --> H[08 Plot body-mass aggregations]
  H --> I[09 Calculate resilience metrics]
  I --> J[10 Plot resilience metrics]
  J --> K[11 Baseline resilience maps]
  I --> L[12 Run NMDS 5 regions]
  L --> M[15 Plot NMDS graphs]
  I --> N[13 FG richness maps]
  N --> O[14 Delta richness maps]
  I --> P[16 FG persistence]

