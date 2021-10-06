args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(derfinder)
  library(grid)
  library(tools)
  library(AnnotationDbi)
  library(ggpubr)
  library(caret)
})

pred <- args[1]
outpath <- args[2]
tissue <- args[3]
pas_dir <- args[4]
txdb_sqlite <- args[5]
scriptPath <- args[6]


source(str_c(scriptPath, "/generate_genomic_state.R"))
source(str_c(scriptPath, "/convert_annot_count_table_to_region_annot.R"))


########################################
############### Functions ##############
########################################

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


txdb_sqlite_paths <- c(txdb_sqlite)
genomic_states_list <- generate_genomic_states_list(txdb_sqlite_paths)


cal_over <- function(er.gr, pac.gr){
  i = -1
  pr.pac = IRanges::subsetByOverlaps(er.gr, pac.gr, ignore.strand = TRUE, maxgap= i, type = "any") %>%
      length()
  return(pr.pac)
}

######################################################
######################################################


print(str_c(Sys.time(), " - ", tissue))
print(str_c(Sys.time(), " - processing PAS data"))

pas_files = list.files(pas_dir, pattern = ".bed", full.names = TRUE)

pas_gr = GRanges()
for(f in pas_files){
  print(f)
  pa = read.table(f,header=F,sep="\t")
  colnames(pa) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
  gr = makeGRangesFromDataFrame(pa, keep.extra.columns = TRUE)
  pas_gr <- c(pas_gr,gr)
}


pac <- IRanges::reduce(pas_gr)
pac = keepStandardChromosomes(pac, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(pac), "UCSC")
pac <- renameSeqlevels(pac, newStyle)
pac <- annotate_ERs_w_genomic_state(ERs_gr = pac, genomic_states_list = genomic_states_list)
pac_intergenic = pac[pac$gtf_TxDb_region_annot %in% "intergenic"]



print(str_c(Sys.time(), " - processing predictions"))


pr = read.table(pred, header=TRUE, sep="\t")


prob_cuts <- c(seq(0.1, 0.9, 0.1))

pos.df <- data.frame()
i <- 0
for(cut in prob_cuts){

  i = i + 1

  pr.pos = pr %>% dplyr::filter(predicted.prob > cut) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

  pos_overlap <- cal_over(pr.pos, pac_intergenic)
  PPV <- round((pos_overlap/length(pr.pos))*100,2)

  pos.df[i, "cutoff"] <- cut
  pos.df[i, "metric"] <- "PPV"
  pos.df[i, "value"] <- PPV
  pos.df[i, "count"] <- pos_overlap
  pos.df[i, "total"] <- length(pr.pos)

}


neg.df <- data.frame()
i <- 0
for(cut in prob_cuts){

  i = i + 1

  pr.neg = pr %>% dplyr::filter(predicted.prob <= cut) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

  neg_overlap <- cal_over(pr.neg, pac_intergenic)
  FOR <- round((neg_overlap/length(pr.neg))*100,2)

  neg.df[i, "cutoff"] <- cut
  neg.df[i, "metric"] <- "FOR"
  neg.df[i, "value"] <- FOR
  neg.df[i, "count"] <- neg_overlap
  neg.df[i, "total"] <- length(pr.neg)

}


res.df = bind_rows(pos.df, neg.df) %>%
  mutate(tissue = tissue, tool = "F3UTER")

write.table(res.df, file = str_c(outpath, "/", tissue, "_f3uter_compare_pas_multiProb.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
