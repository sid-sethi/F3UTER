# Pipeline for reproducing analysis performed on ERs from 39 GTEx tissues

This `snakemake` pipeline reproduces the analysis performed on ERs from 39 GTEx tissues, including all the numbers and the plots presented in the manuscript. Main steps include: (1) processing of ER dataset; (2) calculating omic features; (3) making unannotated 3'UTR predictions using F3UTER; (4) comparing 3'UTR predictions with latest reference gene annotations; (5) categorising ERs by tissue-specificity; (6) performing CNC score analysis.

# Getting Started

## Input

Data from a lot of publicly available resources was used in this analysis. Files which are small in size have been provided in the `/data` directory. However, other files would be required to download/processed and placed in in the `/data` directory.

- ERs in 39 GTEx tissues identified by [Zhang et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32917675/). The ERs can be downloaded from http://rytenlab.com/browser/app/vizER (`all_ERs.csv.gz`). The ERs should be split per tissue and saved in individual files. See `/test_data/ERs_by_tissue` for demo files.
- Ensembl GTF v94
- PhastCons scores from UCSC (bigwig format): https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
- Repeats data from repeatmasker.org: http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm406-dfam2.0/hg38.fa.out.gz

Process the repeatmasker file using the commands below. Please see `/test_data` for a demo file.
```bash
gunzip hg38.fa.out.gz
awk -F " " '{print $5" "$6" "$7" "$11}' hg38.fa.out | sed 1d > hg38.repeatMasker.mod.fa.out
```
- DNA structural properties: data required to calculate DNA structural properties is included in `/data` and is automatically used by the pipeline. These tables were downloaded from http://bioinformatics.psb.ugent.be/webtools/ep3/?conversion.
- GTEx expression data, summarised as mean coverage across the samples in each tissue. The mean coverage should be stored per tissue per chromosome in individual files as a Rle object with the name `tissue_coverage_w_mean_normalised`. The directory structure should look like: `/<tissue_name>/gtex_<tissue_name>_chr<chromosome_name>_mean_cov.rda`.

- Ensembl GTF v92. Place it in `/data/Ensembl_v92.gtf`.
- Ensembl GTF v104. Place it in `/data/Ensembl_v104.gtf`.
- Gencode GTF v38. Place it in `/data/Gencode_v38.gtf`.
- 3'UTRs from RefSeq curated v109 in BED format, downloaded from UCSC Table Browser. Please see `/test_data` for a demo file. Place it in `/data/refseq_curated_v109_three_prime.bed`.
- GTEx median RPKM per gene v6p: https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz. `Gunzip` and place it in `/data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct`.
- OMIM data in csv format. Please see `/test_data` for a demo file. Place it in `/data/OMIM_data.csv`.
- CNC scores for the whole genome. Available on request from [Chen et al.](https://doi.org/10.1038/s41467-021-22262-5). The data needs to be in a GRanges object called as `CNC_gr`. Place it in `/data/CNC_gr.rda`.

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (snakemake R packages) are installed via conda.

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/F3UTER.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd F3UTER/Snakemake/Gtex_analysis
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, you need to make sure that you activate the conda environment first using the command `conda activate pipeline_gtex_analysis`.

```bash
cd F3UTER/Snakemake/Gtex_analysis
conda activate pipeline_gtex_analysis
snakemake --use-conda -j <num_cores> all
```
If you provide more than one core, independent snakemake rules will be processed simultaneously. This pipeline may use up to 50 cores. It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.
