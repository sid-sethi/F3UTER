# F3UTER - Finding 3' UnTranslated Expressed Regions

This repository contains the analysis code for F3UTER. The contents are listed in the order of performed analysis.

## Contents

### 1. [**ER_data_processing**](ER_data_processing)
Processing of ER data from Zhang et al 2020

Script | Description 
------ | -----------
[er_split_by_tissue.R](ER_data_processing/er_split_by_tissue.R) | Split ERs by tissue
[er_initial_processing.pl](ER_data_processing/er_initial_processing.pl) | Add additional gene info to ERs
[er_dataset_generate_for_analysis.R](ER_data_processing/er_dataset_generate_for_analysis.R) | Select intergenic ERs for analysis
[er_3prime_vs_5prime_analysis.R](ER_data_processing/er_3prime_vs_5prime_analysis.R) | Compare 3' intergenic ERs with 5' intergenic ERs

### 2. [**Training_data_generate**](Training_data_generate)

Script | Description
------ | -----------
[training_data_regions_generate.R](Training_data_generate/training_data_regions_generate.R) | Select and process regions for training

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
