# Pipeline for reproducing F3UTER evaluation on drosophila known 3'UTRs

This `snakemake` pipeline extracts known 3'UTRs from dm6 GTF, calculates omic features, makes 3'UTR predictions using F3UTER and evaluates its performance.

# Getting Started

## Input

- Ensembl GTF v104 (`Drosophila_melanogaster.BDGP6.32.104.gtf`)
- PhastCons scores from UCSC (bigwig format): https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phastCons27way/dm6.27way.phastCons.bw
- Repeats data from repeatmasker.org: https://www.repeatmasker.org/genomes/dm6/RepeatMasker-rm406-dfam2.0/dm6.fa.out.gz
- Aligned RNA-seq reads in bigwig format for midgut tissue (n = 3) from the study `SRP197261`.
- Gene TPM abundance for each midgut sample (n = 3) in `SRP197261`calculated using Kallisto.

Process the repeatmasker file using the commands below.
```bash
gunzip dm6.fa.out.gz
awk -F " " '{print $5" "$6" "$7" "$11}' dm6.fa.out | sed 1,3d > dm6.repeatMasker.mod.fa.out
```

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (snakemake and R packages) are installed via conda.

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/F3UTER.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd F3UTER/Snakemake/F3UTER_drosophila_eval
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, you need to make sure that you activate the conda environment first using the command `conda activate pipeline_drosophila_eval`.

```bash
cd F3UTER/Snakemake/F3UTER_drosophila_eval
conda activate pipeline_drosophila_eval
snakemake --use-conda -j <num_cores> all
```
If you provide more than one core, independent snakemake rules will be processed simultaneously. This pipeline uses 6 cores at most. It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.
