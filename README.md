# F3UTER - Finding 3' UnTranslated Expressed Regions

This repository contains the analysis code for F3UTER

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
