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


pred_file <- args[1]
er_file <- args[2]
outpath <- args[3]
tissue <- args[4]
pas_dir <- args[5]
txdb_sqlite <- args[6]
scriptPath <- args[7]


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
  res = data.frame()
  i = -1
  pr.pac = IRanges::subsetByOverlaps(er.gr, pac.gr, ignore.strand = TRUE, maxgap= i, type = "any") %>%
    length()
  res[1, "overlap"] = pr.pac

  return(res)
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


#ER data
er = read.table(er_file, header=T, sep="\t")

## prediction data
pred = read.table(pred_file, header=FALSE, sep="\t") %>%
  mutate(ER_id = str_replace(V4, "^NM_", ""))

pred.group <- pred %>% group_by(ER_id) %>%
  summarise(pcs = str_c(V2, collapse=","), .groups="keep")

pred.er <- left_join(er, pred.group, by = "ER_id")


# calculating getutr response
# adding getutr response
pred.er <- pred.er %>%
  mutate(getutr.response = as.factor(ifelse(!is.na(pcs), "3-UTR", "Non-3-UTR")))


pr.pos = pred.er %>%
  dplyr::filter(getutr.response %in% "3-UTR") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pr.neg = pred.er %>%
  dplyr::filter(getutr.response %in% "Non-3-UTR") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


pas_regions <- pred.er %>%
  dplyr::filter(getutr.response %in% "3-UTR") %>%
  tidyr::separate_rows(pcs, sep=",", convert = TRUE) %>%
  rowwise() %>%
  mutate(pas_point = pcs, start = pas_point - 30, end = pas_point + 30) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pas_ovr_pos = IRanges::subsetByOverlaps(pas_regions, pac_intergenic, ignore.strand = TRUE, maxgap= -1, type = "any") %>%
  .$ER_id %>% unique()

pas_ovr_neg = IRanges::subsetByOverlaps(pr.neg, pac_intergenic, ignore.strand = TRUE, maxgap= -1, type = "any") %>%
  .$ER_id %>% unique()


# adding PAS ground truth
pred.er <- pred.er %>% mutate(pas.truth = as.factor(ifelse(ER_id %in% pas_ovr_pos | ER_id %in% pas_ovr_neg, "3-UTR", "Non-3-UTR")))


confusemat = confusionMatrix(data=pred.er$getutr.response,reference=pred.er$pas.truth, positive = "3-UTR")
results = t(data.frame(cbind(t(round(confusemat$byClass,3)),t(round(confusemat$overall,3))))) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("metric") %>%
  plyr::rename(c("V1" = "GETUTR"))


tp <- confusemat$table[1]
fn <- confusemat$table[2]
fp <- confusemat$table[3]
tn <- confusemat$table[4]
confusetable <- data.frame(metric = c("tp", "fn", "fp", "tn"), value = c(tp, fn, fp, tn))


pos.res3 <- data.frame(overlap = length(pas_ovr_pos),
                       group = "3-UTR",
                       total = length(pr.pos),
                       perc = round((length(pas_ovr_pos)/length(pr.pos))*100,2),
                       comparison_type = "pas_location"
)


neg.res3 <- cal_over(pr.neg, pac_intergenic) %>% mutate(group = "Non-3-UTR", total = length(pr.neg), perc = round((overlap/length(pr.neg))*100,2), comparison_type = "pas_location")

res.df = bind_rows(pos.res3, neg.res3) %>%
  mutate(tissue = tissue, tool = "GETUTR")

write.table(results, file = str_c(outpath, "/", tissue, "_getutr_performance.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table(confusetable, file = str_c(outpath, "/", tissue, "_getutr_confusemat.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table(pred.er %>% dplyr::select(ER_id:pas.truth), file = str_c(outpath, "/", tissue, "_getutr_df.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table(res.df, file = str_c(outpath, "/", tissue, "_getutr_compare_pas.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
