# Pipeline for reproducing F3UTER evaluation on zebrafish known 3'UTRs

This `snakemake` pipeline extracts known 3'UTRs from GRCz10 GTF, calculates omic features, makes 3'UTR predictions using F3UTER and evaluates its performance.

# Getting Started

## Input

- Ensembl GTF v91 (`Danio_rerio.GRCz10.91.gtf`)
- Repeats data from repeatmasker.org: https://www.repeatmasker.org/genomes/danRer10/RepeatMasker-rm406-dfam2.0/danRer10.fa.out.gz
- Aligned RNA-seq reads in bigwig format for liver tissue (n = 2) from the study `SRP213938`.
- Gene TPM abundance for each liver sample (n = 2) in `SRP213938`calculated using Kallisto.

Process the repeatmasker file using the commands below.
```bash
gunzip danRer10.fa.out.gz
awk -F " " '{print $5" "$6" "$7" "$11}' danRer10.fa.out | sed 1,3d > danRer10.repeatMasker.mod.fa.out
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
cd F3UTER/Snakemake/F3UTER_zebrafish_eval
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, you need to make sure that you activate the conda environment first using the command `conda activate pipeline_zebrafish_eval`.

```bash
cd F3UTER/Snakemake/F3UTER_zebrafish_eval
conda activate pipeline_zebrafish_eval
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
