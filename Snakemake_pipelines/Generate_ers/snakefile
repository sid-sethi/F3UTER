import os
from os import path

from snakemake.utils import min_version
min_version("5.3")


configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
BIGWIGDIR = config["bigwig"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]
gtf = config["gtf"]



rule all:
    input:
        "ERs/" + sample + "_ers_raw.txt",
        "Intergenic_ERs/" + sample + "_3prime_intergenic_ers.txt",
        "Intergenic_ERs/" + sample + "_5prime_intergenic_ers.txt"

#########################################################################


rule generate_ERs:
    """
    Call ERs using ODER
    """
    input:
        gtf_file = gtf

    output:
        opt_ers = "ERs/" + sample + "_optimal_ers.rda"

    params:
        bw_dir = BIGWIGDIR,
        resPath = "ERs",
        prefix = sample,
        script = SNAKEDIR + "/scripts/oder.R"

    conda: "envs/r_v4.yml"

    shell:
        """
        Rscript {params.script} {params.bw_dir} {input.gtf_file} {params.prefix} {params.resPath}
        """



rule generate_txdb:
    """
    Generate TxDb object from provided GTF
    """
    input:
        gtf_file = gtf,
        chrom_len = SNAKEDIR + "/data/chromInfo_hg38.txt"

    output:
        txdb = "TxDb/gtf_TxDb.Rsqlite"

    params:
        resPath = "TxDb",
        script = SNAKEDIR + "/scripts/txdb_from_gtf.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.gtf_file} {input.chrom_len} {params.resPath}
        """




rule annotate_ERs:
    """
    Compare ERs with known annotation and annotate ER type
    """
    input:
        opt_ers = rules.generate_ERs.output.opt_ers,
        gtf_txdb = rules.generate_txdb.output.txdb

    output:
        ers = "ERs/" + sample + "_ers_raw.txt"

    params:
        resPath = "ERs",
        tissue = sample,
        script = SNAKEDIR + "/scripts/annotate_ERs.R",
        scriptPath = SNAKEDIR + "/scripts/"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.opt_ers} {input.gtf_txdb} {params.tissue} {params.resPath} {params.scriptPath}
        """




rule generate_intergenic_data:
    """
    Select optimal intergenic ERs (within 10 kb of a protein-coding gene)
    """
    input:
        ers = rules.annotate_ERs.output.ers,
        gtf_file = gtf

    output:
        data_three = "Intergenic_ERs/" + sample + "_3prime_intergenic_ers.txt",
        data_five = "Intergenic_ERs/" + sample + "_5prime_intergenic_ers.txt"

    params:
        resPath = "Intergenic_ERs",
        prefix = sample,
        script = SNAKEDIR + "/scripts/select_intergenic_ers.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {input.gtf_file} {params.resPath} {params.prefix}
        """
