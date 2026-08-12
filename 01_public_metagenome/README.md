# Public fecal metagenome analysis for CRC stage-associated bacteria

This repository contains the main public-data analysis scripts used for the manuscript on Gemella morbillorum and colorectal cancer.

## Analysis overview

This analysis uses eight public colorectal cancer fecal shotgun metagenome cohorts and MetaPhlAn4 species-level relative abundance profiles.

The main comparison is:

- Early-stage CRC: Stage 0 and Stage I
- Advanced-stage CRC: Stage II, Stage III, and Stage IV

The repository includes only the main public-data analyses used in the manuscript:

1. Integrated feature analysis of Stage 0/I versus Stage II/III/IV CRC
2. Feature ranking using Extra Trees and Random Forest
3. Univariate comparison using the Mann-Whitney U test
4. Stage-wise detection rate and mean relative abundance of Gemella morbillorum

## Input data

Input files are placed in:

data/raw/metadata/metadata/
data/raw/metaphlan4_profiles/metaphlan4_profiles/

The following public cohorts are used:

- FengQ_2015
- GuptaA_2019
- VogtmannE_2016
- WirbelJ_2018
- YachidaS_2019
- YangJ_2020
- YuJ_2015
- ZellerG_2014

## Preprocessing

For machine-learning analyses, species-level relative abundance values are transformed using log1p and standardized within each cohort for each species. This cohort-standardized dataset is referred to as cohort_z_log1p.

For univariate Mann-Whitney U tests and stage-wise abundance summaries, raw species-level relative abundance values are used.

## Run

Install dependencies:

pip install -r requirements.txt

Run all analyses:

bash run_all.sh

## Scripts

1. scripts/01_main_stage_feature_selection.py

This script performs the main feature-selection analysis. Extra Trees and Random Forest are applied to cohort_z_log1p-transformed species profiles. Mann-Whitney U tests are performed using raw relative abundance values.

2. scripts/02_gm_stage_detection_abundance.py

This script calculates the detection rate and mean relative abundance of Gemella morbillorum across Healthy, Stage 0, Stage I, Stage II, Stage III, and Stage IV groups.

## Outputs

Feature-selection outputs:

results/top20_extra_trees_cohort_z_log1p_stage01_vs_stage234.csv
results/top20_random_forest_cohort_z_log1p_stage01_vs_stage234.csv
results/top20_mannwhitney_stage01_vs_stage234.csv
results/01_common_top20_features.csv

G. morbillorum stage-wise outputs:

results/02_gm_stage_summary_with_errorbars.csv
results/02_gm_stage_pairwise_tests.csv
results/02_gm_stage_pairwise_detection_fisher.csv
results/02_gm_stage_pairwise_abundance_mannwhitney.csv

Figures:

figures/top20_extra_trees_cohort_z_log1p_stage01_vs_stage234.png
figures/top20_random_forest_cohort_z_log1p_stage01_vs_stage234.png
figures/top20_mannwhitney_stage01_vs_stage234.png
figures/gm_detection_rate_by_healthy_stage.png
figures/gm_mean_abundance_by_healthy_stage.png

## Note on data redistribution

Before uploading input data files to a public GitHub repository, confirm the redistribution policy of each public dataset. If redistribution is not permitted, provide only scripts and instructions for obtaining the input data.
