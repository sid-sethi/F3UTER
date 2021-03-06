# F3UTER - Finding 3' UnTranslated Expressed Regions

This repository contains the analysis code for F3UTER.

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.


## Contents (Intructions for use)

The contents are listed in the order of performed analysis. The output of these scripts are txt files. 

### 1. [**ER_data_processing**](ER_data_processing)
Processing of ER data from Zhang et al 2020

Script | Description 
------ | -----------
[er_split_by_tissue.R](ER_data_processing/er_split_by_tissue.R) | Split ERs by tissue
[er_initial_processing.pl](ER_data_processing/er_initial_processing.pl) | Add additional gene info to ERs
[er_dataset_generate_for_analysis.R](ER_data_processing/er_dataset_generate_for_analysis.R) | Select intergenic ERs for analysis
[er_3prime_vs_5prime_analysis.R](ER_data_processing/er_3prime_vs_5prime_analysis.R) | Compare 3' intergenic ERs with 5' intergenic ERs

An example dataset of processed 3' intergenic ERs is provided within [Test_dataset](Test_dataset).

### 2. [**Training_data_generate**](Training_data_generate)

Script | Description
------ | -----------
[training_data_regions_generate.R](Training_data_generate/training_data_regions_generate.R) | Select and process regions for training

The training data can be downloaded from the F3UTER app: https://astx.shinyapps.io/F3UTER/

### 3. [**OmicFeature_generation**](OmicFeature_generation)
Calculating omic features, can be applied to training regions and ERs

Script | Description
------ | -----------
[polyA_signal.R](OmicFeature_generation/polyA_signal.R) | Scan for Poly(A) signal (number of features, n=1)
[nucleotide_frequency.R](OmicFeature_generation/nucleotide_frequency.R) | Calculate mono- and di-nucleotide frequency (n=20)
[sequence_conservation.R](OmicFeature_generation/sequence_conservation.R) | Calculate mean phastCons score (n=1)
[transposons_overlap.R](OmicFeature_generation/transposons_overlap.R) | Calculate overlap with transposons ((n=1)
[expression_features.R](OmicFeature_generation/expression_features.R) | Calculate entropy efficiency and percentage difference of expression reads (n=2)
[dna_structural_features.R](OmicFeature_generation/dna_structural_features.R) | Calculate DNA structural properties (n=16)

The training and ER data feature matrix can be downloaded from the F3UTER app: https://astx.shinyapps.io/F3UTER/

### 4. [**Uni-Multi_variate_analysis**](Uni-Multi_variate_analysis)
Univariate and multivariate analysis of omic features

Script | Description
------ | -----------
[make_ml_table.R](Uni-Multi_variate_analysis/make_ml_table.R) | Compile feature matrix
[univariate_analysis.R](Uni-Multi_variate_analysis/univariate_analysis.R) | Perform univariate analysis on features
[umap_analysis.R](Uni-Multi_variate_analysis/umap_analysis.R) | Perform UMAP analysis on features
[run_classification_models.pl](Uni-Multi_variate_analysis/run_classification_models.pl) | Wrapper for running multinomial classification models
[summarise_multinomial_models.R](Uni-Multi_variate_analysis/summarise_multinomial_models.R) | Wrapper for summarising results from multinomial models
[plot_multinomial_models.R](Uni-Multi_variate_analysis/plot_multinomial_models.R) | Plot results from multinomial models
[elasticNetMultinomialLR_model.R](Uni-Multi_variate_analysis/elasticNetMultinomialLR_model.R) | code for elastic net multinomial logistic regression
[randomForestMultinomial_model.R](Uni-Multi_variate_analysis/randomForestMultinomial_model.R) | code for multinomial random forest classification
[RF_vs_LR_analysis.R](Uni-Multi_variate_analysis/RF_vs_LR_analysis.R) | Compare multinomial logistic regression and random forest results

### 5. [**F3UTER_evaluation**](F3UTER_evaluation)
Construction and cross-validation of F3UTER

Script | Description
------ | -----------
[run_f3uter_cv.pl](F3UTER_evaluation/run_f3uter_cv.pl) | Wrapper for running F3UTER cross validation evaluation
[summarise_f3uter_cv.R](F3UTER_evaluation/summarise_f3uter_cv.R) |  Wrapper for summarising cross validation results
[f3uter_cv.R](F3UTER_evaluation/f3uter_cv.R) | Code for training F3UTER
[plot_roc_pr_curves.R](F3UTER_evaluation/plot_roc_pr_curves.R) | Plot ROC and precision-recall curves
[f3uter_trained_model.R](F3UTER_evaluation/f3uter_trained_model.R) | Save F3UTER trained model

### 6. [**Bcell_validation**](Bcell_validation)
Validation of F3UTER predictions using RNA-seq and 3'-seq data in B cells

Script | Description
------ | -----------
[generate_Bcell_predictions.R](Bcell_validation/generate_Bcell_predictions.R) | Generate 3'UTR predictions in B cell ER dataset using F3UTER
[merge_erData_with_predictionData.R](Bcell_validation/merge_erData_with_predictionData.R) | Merge ER prediction data with ER raw meta-data
[permute_random_intergenic_ERs.pl](Bcell_validation/permute_random_intergenic_ERs.pl) | Generate randomly selected intergenic ERs for permutation test
[regions_to_exclude_for_permutation.R](Bcell_validation/regions_to_exclude_for_permutation.R) | Genomic space to mask in order to produce intergenic ER space
[compare_knownThreePrime_with_polya.R](Bcell_validation/compare_knownThreePrime_with_polya.R) | Compare known 3'UTRs with poly(A) site clusters
[compare_BcellErs_with_polya.R](Bcell_validation/compare_BcellErs_with_polya.R) | Compare 3'UTR predictions in B cells with poly(A) site clusters

### 7. [**GTEx_predictions**](GTEx_predictions)
Applying F3UTER to ERs derived from GTEx tissues to predict unannotated 3'UTRs. All the ER predictions can be downloaded from the F3UTER app: https://astx.shinyapps.io/F3UTER/

Script | Description
------ | -----------
[generate_gtex_predictions.R](GTEx_predictions/generate_gtex_predictions.R) | Use F3UTER on GTEx ERs to produce predictions
[merge_gtexErData_with_predictions.R](GTEx_predictions/merge_gtexErData_with_predictions.R) | Merge ER prediction data with ER raw meta-data
[calculate_gtex_prediction_numbers.R](GTEx_predictions/calculate_gtex_prediction_numbers.R) | Calculate basic stats/numbers for prediction results across tissues
[categorise_er_tissue_specificity.R](GTEx_predictions/categorise_er_tissue_specificity.R) | Split ER predictions based on their tissue-specificity across 39 tissues
[generate_categorised_geneLists_and_tables.R](GTEx_predictions/generate_categorised_geneLists_and_tables.R) | Generate data tables and gene lists for each category to be used for downstream analysis
[calculate_cncr_score.R](GTEx_predictions/calculate_cncr_score.R) | Code to calculate CNC scores for a query region
[perform_cncr_analysis.R](GTEx_predictions/perform_cncr_analysis.R) | Perform CNCR analysis on all predictions and tissue-specific groups

### 8. [**GTEx_predictions_RBP_analysis**](GTEx_predictions_RBP_analysis)

Calculate RBP motif enrichment within ER predictions. See [**README.md**](GTEx_predictions_RBP_analysis/README.md) for details.

### 9. [**F3UTER_app**](F3UTER_app)

Source code for the F3UTER online resource.
App URL: https://astx.shinyapps.io/F3UTER/

## System requirements for data analysis

All the data analysis was performed in R version 3.6.2. For analysis or plotting, the following packages were used: ggplot2_2_3.3.2, ggridges_0.5.2, ggridges_0.5.2, ggsignif_0.6.0, rstatix_0.6.0, derfinder_1.20.0, rtracklayer_1.46.0, GenomicFeatures_1.38.2, GenomicRanges_1.38.0, ggpubr_0.4.0, stringr_1.4.0, dplyr_1.0.2, tidyverse_1.3.0

