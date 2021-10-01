import os
from os import path

from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

er_file = config["ers"]
sample = config["sample_name"]


chr = list(range(1,23))
chr.append("X")
chr.append("Y")
chr.append("MT")

# ----------------------------------------------------------------

rule all:
    input:
        "F3UTER_features/" + sample + "_nt_freq.txt",
        "F3UTER_features/" + sample + "_phastcons.txt",
        "F3UTER_features/" + sample + "_polyA_signal.txt",
        "F3UTER_features/" + sample + "_structural_feat.txt",
        "F3UTER_features/" + sample + "_repeats.txt",
        "F3UTER_features/" + sample + "_exp_feat.txt"



#########################################################################

rule cal_nt_freq:
    input:
        ers = er_file

    output:
        out = "F3UTER_features/" + sample + "_nt_freq.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/nt_freq.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix}
        """


rule cal_phastcons:
    input:
        ers = er_file,
        bw = config["phastcons"]

    output:
        out = "F3UTER_features/" + sample + "_phastcons.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/phastCons.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {input.bw}
        """



rule cal_pas:
    input:
        ers = er_file

    output:
        out = "F3UTER_features/" + sample + "_polyA_signal.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/polyA_signal.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix}
        """



rule cal_structural_feat:
    input:
        ers = er_file

    output:
        out = "F3UTER_features/" + sample + "_structural_feat.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/structural_feat.R",
        conversion_tables = SNAKEDIR + "/data/dna_structural_conversion_tables"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {params.conversion_tables}
        """



rule cal_repeats:
    input:
        ers = er_file,
        repeats_data = config["repeats"]

    output:
        out = "F3UTER_features/" + sample + "_repeats.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/repeats.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {input.repeats_data}
        """



rule generate_coverage:
    input:
        chrom_info = SNAKEDIR + "/data/hg38_v2.chrom.sizes"

    output:
        out = expand("Coverage_files/" + sample + "_{seqnames}_mean_cov.rda", seqnames=chr)

    params:
        bw_dir = config["bigwig_dir"],
        resPath = "Coverage_files",
        prefix = sample,
        script = SNAKEDIR + "/scripts/generate_coverage.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {params.bw_dir} {params.resPath} {params.prefix} {input.chrom_info}
        """



rule cal_exp_feat:
    input:
        ers = er_file,
        coverage_files = rules.generate_coverage.output.out

    output:
        out = "F3UTER_features/" + sample + "_exp_feat.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/exp_features.R",
        coverage_dir = "Coverage_files"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {params.coverage_dir}
        """
