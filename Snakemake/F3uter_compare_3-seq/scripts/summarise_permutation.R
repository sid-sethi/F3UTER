args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(derfinder)
  library(AnnotationDbi)
  library(ggpubr)
  library(grid)
  library(tools)
})


pas_dir <- args[1]
reg_dir <- args[2] # _ShuffledRegions directory
f3uter_res <- args[3]
outpath <- args[4]
tissue <- args[5]
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

############################################
############################################


print(str_c(Sys.time(), " - ", tissue))
print(str_c(Sys.time(), " - processing PAS data"))

pas_files = list.files(pas_dir, pattern = ".bed", full.names = TRUE)

pas_gr = GRanges()
for(f in pas_files){
  print(f)
  pa = read.table(pas_files[2],header=F,sep="\t")
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



print(str_c(Sys.time(), " - compiling ", tissue, " permutation data"))

files <- list.files(path = reg_dir, full.names = TRUE)

simRes = data.frame()
i = 0
for(f in files){
  i = i + 1
  simBed = read.table(f, header = FALSE, sep = "\t")
  colnames(simBed) = c("seqnames", "start", "end")
  simBed.gr = makeGRangesFromDataFrame(simBed, keep.extra.columns = TRUE)

  overlap = IRanges::subsetByOverlaps(simBed.gr, pac_intergenic, ignore.strand = TRUE, maxgap= -1, type = "any") %>% length()

  simRes[i, "dist"] = basename(f)
  simRes[i, "number"] = overlap
  simRes[i, "group"] = "Random"
  simRes[i, "total"] = length(simBed.gr)
  simRes[i, "perc"] = round((overlap/length(simBed.gr))*100,2)
  simRes[i, "tissue"] = tissue

}


write.table(simRes, file = str_c(outpath, "/", tissue, "_simRes.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")



ntimes = 10000

f3uter.df <- read.table(f3uter_res, header=TRUE, sep="\t")
o_obs = f3uter.df %>% filter(group %in% "3-UTR") %>% .$perc

o_perm = mean(simRes$perc)
perm_high_count = simRes %>% filter(perc > o_obs) %>% nrow()
pval = ifelse(perm_high_count == 0, 0.0001,
              perm_high_count/ntimes
)
zscore = round((o_obs - o_perm)/sd(simRes$perc) , 2)
label = paste(paste("p-value < ",pval),paste("Z-score: ",zscore),paste("n perm: ",ntimes), sep="\n")

# permutation plot
p4 <- gghistogram(simRes, x = "perc", y = "..count..",
                  binwidth = 0.2,
                  fill = "#66B2FF", color = "#66B2FF",
                  title = tissue
)

png(str_c(outpath, "/", tissue, "_simRes.png"), width = 4.25, height = 3.75, res = 600, units = "in")
p4 +
  theme(plot.title = element_text(size=12, hjust=0.5),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  scale_x_continuous(name = "% of overlap with intergenic poly(A) clusters", limits = c(0,o_obs)) +
  scale_y_continuous(name = "Frequency") +
  geom_vline(xintercept = o_obs, color = "black", linetype = "dashed") +
  geom_vline(xintercept = mean(simRes$perc), color = "black", linetype = "dashed")+
  annotate("text",x=o_obs*0.91, y = max(ggplot_build(p4)$data[[1]]$count)*0.96, label = expression(atop(O[obs], "F3UTER")), parse=T, size = 3, color = "red") +
  annotate("text",x=mean(simRes$perc), y = max(ggplot_build(p4)$data[[1]]$count)*0.96, label = expression(O[perm]), parse=T, size = 3, color="red") +
  annotate("text",x=o_obs*0.52, y = max(ggplot_build(p4)$data[[1]]$count)*0.75, label = label, parse=F, size = 3)
dev.off()
