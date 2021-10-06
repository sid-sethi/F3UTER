# Pipeline for producing poly(A) site predictions within ERs in B cells using APARENT, TAPAS and GETUTR

This `snakemake` pipeline modifies ERs and applies poly(A) site prediction tools to predict novel poly(A) sites.

# Getting Started

## Input
- 3' intergenic ERs
- Aligned RNA-seq reads in BAM format. Multiple replicates were pooled together.

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- GETUTR: http://big.hanyang.ac.kr/GETUTR/download.htm
- The rest of the dependencies (snakemake, samtools and R packages) are installed via conda.

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/F3UTER.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd F3UTER/Snakemake/PAS_predictions
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, you need to make sure that you activate the conda environment first using the command `conda activate pipeline_pas_predictions`.

```bash
cd F3UTER/Snakemake/PAS_predictions
conda activate pipeline_pas_predictions
snakemake --use-conda -j <num_cores> all
```
If you provide more than one core, independent snakemake rules will be processed simultaneously. This pipeline can use up to 5 cores. It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.
