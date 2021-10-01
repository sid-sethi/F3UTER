import os
from os import path

from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

features_dir = config["features_dir"]
sample = config["sample_name"]

# ----------------------------------------------------------------

rule all:
    input:
        "F3UTER_prediction/" + sample + "_predictions_meta.txt"



#########################################################################

rule make_ml_table:
    output:
        ml_table = "F3UTER_prediction/" + sample + "_ml_table.rds"

    params:
        feat_dir = features_dir,
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/make_ML_table.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {params.feat_dir} {params.resPath} {params.prefix}
        """


rule make_f3uter_predictions:
    input:
        ml_table = rules.make_ml_table.output.ml_table,
        prediction_model = SNAKEDIR + "/data/rf_model.rda"

    output:
        preds = "F3UTER_prediction/" + sample + "_predictions.txt"

    params:
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/f3uter_predict.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ml_table} {params.resPath} {params.prefix} {input.prediction_model}
        """



rule add_metadata_predictions:
    input:
        preds = rules.make_f3uter_predictions.output.preds,
        er_file = config["ers"]

    output:
        meta = "F3UTER_prediction/" + sample + "_predictions_meta.txt",
        bed = "F3UTER_prediction/" + sample + "_predictions.bed"

    params:
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/add_metadata_predictions.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.preds} {input.er_file} {params.resPath} {params.prefix}
        """
