args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(GenomicFeatures)
  library(derfinder)
  library(data.table)
})
# Set WD ----------------------------------------------------------------------------------------------

oc_maxgap = args[1]
gtf_txdb = args[2]
tissue = args[3]
resPath = args[4]
scriptPath = args[5]

# Load data -------------------------------------------------------------------------------------------

load(oc_maxgap)

# Functions -------------------------------------------------------------------------------------------

generate_genomic_states_list <- function(txdb_sqlite_paths){

  genomic_states_list <- list()

  for(i in seq_along(txdb_sqlite_paths)){

    txdb_sqlite_path <- txdb_sqlite_paths[i]


    print(str_c(Sys.time(), " - ", i, " - ", txdb_sqlite_path))

    txdb_ref_annot <- loadDb(file = txdb_sqlite_path)

    genomic_state_ref_annot <-
      generate_genomic_state(txdb_ref_annot,
                             output_path =
                               txdb_sqlite_path %>%
                               str_replace("\\.Rsqlite", "_genomic_state.rda")
                               )

    genomic_states_list <-
      c(genomic_states_list, genomic_state_ref_annot)

  }

  names(genomic_states_list) <-
    txdb_sqlite_paths %>%
    basename() %>%
    str_replace("\\.Rsqlite", "")

  return(genomic_states_list)

}



annotate_ERs_w_genomic_state <- function(ERs_gr, genomic_states_list){

  for(j in 1:length(genomic_states_list)){

    reference_annot_ver <- genomic_states_list[j] %>% names()
    genomic_state <- genomic_states_list[[j]]

    print(str_c(Sys.time(), " - ", j, " - ", reference_annot_ver))

    ER_annotation_count_table <-
      annotateRegions(regions = ERs_gr,
                      genomicState = genomic_state$fullGenome,
                      maxgap = -1L, minoverlap = 1L,
                      annotate = F)

    ER_annotation_count_table_w_region_annot <-
      convert_annot_count_table_to_region_annot(count_table = ER_annotation_count_table$countTable)

    stopifnot(nrow(ER_annotation_count_table_w_region_annot) == length(ERs_gr))

    elementMetadata(ERs_gr)[[str_c(reference_annot_ver, "_region_annot")]] <-
      ER_annotation_count_table_w_region_annot[["region_annot"]]



  }

  return(ERs_gr)

}



##### Second Level #####

source(str_c(scriptPath, "/generate_genomic_state.R"))
source(str_c(scriptPath, "/convert_annot_count_table_to_region_annot.R"))
source(str_c(scriptPath, "/mark_overlapping_genes_gr.R"))


# Main ------------------------------------------------------------------------------------------------

# generate a list of genomic states from different reference annotations ready to annotate the ERs
txdb_sqlite_paths <- c(gtf_txdb)
genomic_states_list <- generate_genomic_states_list(txdb_sqlite_paths)


ensembl_grch38_v92_TxDb <- loadDb(file = gtf_txdb)


for(i in 1:length(list_ERs_each_tissue_optimal_cut_off_maxgap)){

  cut_off_tissue <- names(list_ERs_each_tissue_optimal_cut_off_maxgap[i])

  ERs_one_tissue_optimal_cut_off_no_scaffold <-
    list_ERs_each_tissue_optimal_cut_off_maxgap[[i]] %>%
    keepSeqlevels(str_c("chr", c(1:22, "X", "Y")), pruning.mode = "coarse") %>%
    sortSeqlevels() %>%
    sort()

  print(str_c(Sys.time(), " - ", i, " - ", names(list_ERs_each_tissue_optimal_cut_off_maxgap[i])))

  print(str_c(Sys.time(), " - adding ER index"))

  elementMetadata(ERs_one_tissue_optimal_cut_off_no_scaffold)[["ER_index"]] <- 1:length(ERs_one_tissue_optimal_cut_off_no_scaffold)



  ERs_w_annotation <- ERs_one_tissue_optimal_cut_off_no_scaffold


  print(str_c(Sys.time(), " - annotating ERs by annotation feature"))

  ERs_w_annotation <-
    annotate_ERs_w_genomic_state(ERs_gr = ERs_w_annotation, genomic_states_list = genomic_states_list)


  print(str_c(Sys.time(), " - adding tissue"))

  elementMetadata(ERs_w_annotation)[["tissue"]] <- tissue

  ERs_w_annotation_df <- ERs_w_annotation %>% as.data.frame()


  write.table(ERs_w_annotation_df, file = str_c(resPath, "/", tissue, "_ers_raw.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

}
