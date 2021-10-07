# Pipeline for calculating RBP density enrichment in unannotated 3'UTR predictions

This `snakemake` pipeline calculates and plots RBP density enrichment in known 3'UTRs, predicted 3'UTRs, predicted non-3'UTRs and di-nucleotide shuffled predicted 3'UTRs.

# Getting Started

## Input

- Requires results from [Gtex_analysis](https://github.com/sid-sethi/F3UTER/tree/master/Snakemake/Gtex_analysis).

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- MEME suite (v5.3.3): https://meme-suite.org/meme/
- The rest of the dependencies (snakemake and R packages) are installed via conda.

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/F3UTER.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd F3UTER/Snakemake/RBP_analysis
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, you need to make sure that you activate the conda environment first using the command `conda activate pipeline_rbp_analysis`.

```bash
cd F3UTER/Snakemake/RBP_analysis
conda activate pipeline_rbp_analysis
snakemake --use-conda -j <num_cores> all
```
If you provide more than one core, independent snakemake rules will be processed simultaneously. This pipeline can use up to 6 cores. It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.
