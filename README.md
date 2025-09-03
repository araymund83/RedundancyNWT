## RedundancyNWT
This repository contains scripts and outputs for the functional resilience and turnover analysis of boreal landbird communities under projected climate and habitat change in the Taiga Plains in the Northwest Territories in Canada. Proposed in the manuscript entitled:Northern landbird communities show geographically constrained resilience to projected climate and forest changes.
The workflow integrates functional group (FG) binary matrices, temporal beta diversity (TBI), body-mass discontinuity analyses, redundancy metrics and Non-Metric Multidimensional Scaling (NMDS) analyses.

## 📂 Repository Structure

```text
inputs/     # Raw or input data
tables/     # Processed tabular outputs (e.g., FG binary matrices, log mass matrices)
outputs/    # Intermediate outputs (TBI, discontinuities, redundancy, NMDS)
scripts/    # All R scripts used in the workflow


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


```markdown
## 🔄 Workflow (Stepwise)

| Step | Script | Purpose | Outputs |
|------|--------|---------|---------|
| **01** | `01_build_fg_binary_matrices.R` | Build FG binary matrices | `tables/fg_bin/2025/` |
| **02** | `02_calculate_fg_tbi.R` | Calculate TBI for FGs | `outputs/tbiFGBin_1km/2025/` |
| **03** | `03_generate_fg_tbi_rasters.R` | Create TBI rasters & tables (TBI, p-value, change) | `figs/TBI_FG_Ras_72spp/2025/`, `tables/tbiFGBin/2025/` |
| **04** | `04_plot_fg_tbi_changes.R` | Plot per-GCM TBI change maps | – |
| **05** | `05_summarize_fg_tbi_changes.R` | Mean ΔTBI & significance across GCMs | `figs/TBI_FG_changemean/2025/` |
| **06** | `06_build_logmass_matrices.R` | Build log10 body-mass matrices | `tables/logMassMat/` |
| **07** | `07_calculate_bodymass_discontinuities.R` | Compute mass aggregations/discontinuities | `outputs/Discontinuities/` |
| **08** | `08_plot_bodymass_aggregations.R` | Plot aggregation panels (per GCM/year) | `figs/aggregations2025_2/` |
| **09** | `09_calculate_resilience_metrics.R` | Cross-scale redundancy, within-scale redundancy, cross-scale diversity | `outputs/Redundancy/2025/` |
| **09a** | `09a_within_scale_redundancy.R` | Helper for within-scale redundancy | `outputs/Redundancy/2025/` |
| **10** | `10_plot_resilience_metrics.R` | Plot resilience metric maps | – |
| **11** | `11_plot_baseline_resilience_metrics.R` | 2011 baseline maps | `figs/Redundancy/2011_map_*.png` |
| **12** | `12_run_nmds_5regions.R` | NMDS + PERMANOVA (5 regions; all years) | `outputs/nmds_results_2025/`, `outputs/bootstrapping_2025/`, `outputs/permanova_2025/` |
| **13** | `13_plot_fg_richness_maps.R` | FG richness per year/GCM | `figs/FG_richness/2025/` |
| **14** | `14_plot_delta_richness_maps.R` | Δ richness (2011→2031, 2011→2091) | `figs/richness/2025/` |
| **15** | `15_plot_nmds_graphs_5regions.R` | Plot NMDS ordination graphs | `figs/nmds_graphs/` |
| **16** | `16_calculate_fg_persistence.R` | Calculate FG persistence | TBD |


