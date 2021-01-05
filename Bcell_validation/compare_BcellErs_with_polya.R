library(tidyverse)
library(rtracklayer)
library(derfinder)
library(grid)
library(tools)
library(AnnotationDbi)
library(ggpubr)

# makePredictions directory
path = "/Bcells"

if(!dir.exists(paste(path, "/Intergenic_3/Validation", sep=""))){
  system(paste("mkdir -m a=rwx ",path, "/Intergenic_3/Validation", sep=""))
}

source("/generate_genomic_state.R")
source("/convert_annot_count_table_to_region_annot.R")

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
    str_replace("/.*/", "") %>%
    str_replace("\\.Rsqlite", "")

  return(genomic_states_list)

}


txdb_sqlite_paths <- c("/gtf_TxDb_v92.Rsqlite")
genomic_states_list <- generate_genomic_states_list(txdb_sqlite_paths)


cal_over <- function(er.gr){
  res = data.frame()
  j=0
  for(i in c(-1, seq(5, 300, 5))) {
    j= j+1
    pr.pac = IRanges::subsetByOverlaps(er.gr, pac, ignore.strand = TRUE, maxgap= i, type = "any") %>%
      length()

    res[j, "dist"] = i
    res[j, "number"] = pr.pac
  }
  return(res)
}

######################################################
######################################################


# B cells
# requires tissue name hardcoded
print(str_c(Sys.time(), " - GCB1"))
## GCB1
tissue = "GCB1"

print(str_c(Sys.time(), " - processing GCB1 poly(A) data"))
# polya site data
pa = read.table("/GSE111310/GCB1/Polya_clusters/GCB1.clusters.hg38.2-0.bed",header=F,sep="\t")
colnames(pa) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
pac = makeGRangesFromDataFrame(pa, keep.extra.columns = TRUE)
pac = keepStandardChromosomes(pac, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(pac), "UCSC")
pac <- renameSeqlevels(pac, newStyle)
pac <- annotate_ERs_w_genomic_state(ERs_gr = pac, genomic_states_list = genomic_states_list)
pac = pac[pac$gtf_TxDb_v92_region_annot %in% "intergenic"]


print(str_c(Sys.time(), " - processing GCB1 ER predictions"))
pr = readRDS(str_c(path, "/PredictionMerge/", tissue, "_mergedData.rds"))

pr.pos = pr %>% dplyr::filter(predicted.prob > 0.60) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
pr.neg = pr %>% dplyr::filter(predicted.prob <= 0.60) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pos.res <- cal_over(pr.pos) %>% mutate(group = "3-UTR", total = length(pr.pos), perc = round((number/length(pr.pos))*100,2))
neg.res <- cal_over(pr.neg) %>% mutate(group = "Non-3-UTR", total = length(pr.neg), perc = round((number/length(pr.neg))*100,2))

df.gcb1 = bind_rows(pos.res, neg.res) %>%
  mutate(group2 = "B cells Rep1")


###################################################
####### compiling permutation results #############
###################################################

print(str_c(Sys.time(), " - compiling GCB1 permutation data"))

files <- list.files(path = str_c(path, "/Intergenic_3/Permutation/GCB1_ShuffledRegions"), full.names = TRUE)

simRes.gcb1 = data.frame()
i = 0
for(f in files){
  i = i + 1
  simBed = read.table(f, header = FALSE, sep = "\t")
  colnames(simBed) = c("seqnames", "start", "end")
  simBed.gr = makeGRangesFromDataFrame(simBed, keep.extra.columns = TRUE)

  overlap = IRanges::subsetByOverlaps(simBed.gr, pac, ignore.strand = TRUE, maxgap= -1, type = "any") %>% length()

  simRes.gcb1[i, "dist"] = basename(f)
  simRes.gcb1[i, "number"] = overlap
  simRes.gcb1[i, "group"] = "Random"
  simRes.gcb1[i, "total"] = length(simBed.gr)
  simRes.gcb1[i, "perc"] = round((overlap/length(simBed.gr))*100,2)
  simRes.gcb1[i, "group2"] = "B cells Rep1"

}

ntimes = 10000
o_obs = df.gcb1 %>% filter(dist == -1, group %in% "3-UTR") %>% .$perc
o_perm = mean(simRes.gcb1$perc)
perm_high_count = simRes.gcb1 %>% filter(perc > o_obs) %>% nrow()
pval = ifelse(perm_high_count == 0, 0.0001,
              perm_high_count/ntimes
)
zscore = round((o_obs - o_perm)/sd(simRes.gcb1$perc) , 2)
label = paste(paste("p-value < ",pval),paste("Z-score: ",zscore),paste("n perm: ",ntimes), sep="\n")

# permutation plot
p3 <- gghistogram(simRes.gcb1, x = "perc", y = "..count..",
                  binwidth = 0.2,
                  fill = "#66B2FF", color = "#66B2FF",
                  title = "Bcells Rep1"
)

png(str_c(path, "/Intergenic_3/Validation/", tissue, "_permutation.png"), width = 4.25, height = 3.75, res = 600, units = "in")
p3 +
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
  geom_vline(xintercept = mean(simRes.gcb1$perc), color = "black", linetype = "dashed")+
  annotate("text",x=o_obs*0.91, y = max(ggplot_build(p3)$data[[1]]$count)*0.96, label = expression(atop(O[obs], "predicted 3'UTRs")), parse=T, size = 3, color = "red") +
  annotate("text",x=mean(simRes.gcb1$perc), y = max(ggplot_build(p3)$data[[1]]$count)*0.96, label = expression(O[perm]), parse=T, size = 3, color="red") +
  annotate("text",x=o_obs*0.52, y = max(ggplot_build(p3)$data[[1]]$count)*0.75, label = label, parse=F, size = 3)
dev.off()





print(str_c(Sys.time(), " - GCB2"))
## GCB2
rm(pa, pr)
tissue = "GCB2"

print(str_c(Sys.time(), " - processing GCB2 poly(A) data"))
# polya site data
pa = read.table("/GSE111310/GCB2/Polya_clusters/GCB2.clusters.hg38.2-0.bed",header=F,sep="\t")
colnames(pa) = c("chr","start","end","id","exp_thisSample","strand", "perc_sample_support", "number_protocol_support", "avg_exp_allSamples", "cluster_annotation","pas")
pac = makeGRangesFromDataFrame(pa, keep.extra.columns = TRUE)
pac = keepStandardChromosomes(pac, species = "Homo_sapiens", pruning.mode="coarse")
newStyle <- mapSeqlevels(seqlevels(pac), "UCSC")
pac <- renameSeqlevels(pac, newStyle)
pac <- annotate_ERs_w_genomic_state(ERs_gr = pac, genomic_states_list = genomic_states_list)
pac = pac[pac$gtf_TxDb_v92_region_annot %in% "intergenic"]


print(str_c(Sys.time(), " - processing GCB2 ER predictions"))
pr = readRDS(str_c(path, "/Intergenic_3/PredictionApp/", tissue, "_appData.rds"))

pr.pos = pr %>% dplyr::filter(predicted.prob > 0.60) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
pr.neg = pr %>% dplyr::filter(predicted.prob <= 0.60) %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

pos.res <- cal_over(pr.pos) %>% mutate(group = "3-UTR", total = length(pr.pos), perc = round((number/length(pr.pos))*100,2))
neg.res <- cal_over(pr.neg) %>% mutate(group = "Non-3-UTR", total = length(pr.neg), perc = round((number/length(pr.neg))*100,2))

df.gcb2 = bind_rows(pos.res, neg.res) %>%
  mutate(group2 = "B cells Rep2")


###################################################
####### compiling permutation results #############
###################################################
print(str_c(Sys.time(), " - compiling GCB2 permutation data"))

files <- list.files(path = str_c(path, "/Intergenic_3/Permutation/GCB2_ShuffledRegions"), full.names = TRUE)

simRes.gcb2 = data.frame()
i = 0
for(f in files){
  i = i + 1
  simBed = read.table(f, header = FALSE, sep = "\t")
  colnames(simBed) = c("seqnames", "start", "end")
  simBed.gr = makeGRangesFromDataFrame(simBed, keep.extra.columns = TRUE)

  overlap = IRanges::subsetByOverlaps(simBed.gr, pac, ignore.strand = TRUE, maxgap= -1, type = "any") %>% length()

  simRes.gcb2[i, "dist"] = basename(f)
  simRes.gcb2[i, "number"] = overlap
  simRes.gcb2[i, "group"] = "Random"
  simRes.gcb2[i, "total"] = length(simBed.gr)
  simRes.gcb2[i, "perc"] = round((overlap/length(simBed.gr))*100,2)
  simRes.gcb2[i, "group2"] = "B cells Rep2"

}

ntimes = 10000
o_obs = df.gcb2 %>% filter(dist == -1, group %in% "3-UTR") %>% .$perc
o_perm = mean(simRes.gcb2$perc)
perm_high_count = simRes.gcb2 %>% filter(perc > o_obs) %>% nrow()
pval = ifelse(perm_high_count == 0, 0.0001,
              perm_high_count/ntimes
              )
zscore = round((o_obs - o_perm)/sd(simRes.gcb2$perc) , 2)
label = paste(paste("p-value < ",pval),paste("Z-score: ",zscore),paste("n perm: ",ntimes), sep="\n")

# permutation plot
p4 <- gghistogram(simRes.gcb2, x = "perc", y = "..count..",
            binwidth = 0.2,
            fill = "#66B2FF", color = "#66B2FF",
            title = "Bcells Rep2"
)

png(str_c(path, "/Intergenic_3/Validation/", tissue, "_permutation.png"), width = 4.25, height = 3.75, res = 600, units = "in")
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
geom_vline(xintercept = mean(simRes.gcb2$perc), color = "black", linetype = "dashed")+
annotate("text",x=o_obs*0.91, y = max(ggplot_build(p4)$data[[1]]$count)*0.96, label = expression(atop(O[obs], "predicted 3'UTRs")), parse=T, size = 3, color = "red") +
annotate("text",x=mean(simRes.gcb2$perc), y = max(ggplot_build(p4)$data[[1]]$count)*0.96, label = expression(O[perm]), parse=T, size = 3, color="red") +
annotate("text",x=o_obs*0.52, y = max(ggplot_build(p4)$data[[1]]$count)*0.75, label = label, parse=F, size = 3)
dev.off()




# combining GCB1 and GCB2 data
print(str_c(Sys.time(), " - plotting combined plots...."))

x = bind_rows(
  df.gcb1 %>% dplyr::filter(dist == -1) %>% mutate(dist = as.character(dist)),
  df.gcb2 %>% dplyr::filter(dist == -1)%>% mutate(dist = as.character(dist)),
  simRes.gcb1,
  simRes.gcb2
)

x.summarise = x %>% group_by(group, group2) %>%
  summarise(mean = mean(perc))

x.summarise$group <- ifelse(x.summarise$group %in% "3-UTR", "Predicted 3-UTRs", x.summarise$group)
x.summarise$group <- ifelse(x.summarise$group %in% "Non-3-UTR", "Predicted non-3-UTRs", x.summarise$group)


png(str_c(path, "/Intergenic_3/Validation/", "GCB", "_polyaOverlap.png"), width = 4.25, height = 3.75, res = 600, units = "in")
ggplot(data=x.summarise, aes(x= group2, y=mean, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  theme_bw() +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "gray40")) +
  theme(plot.title = element_text(size=11, hjust=0.5),
        legend.key.size = unit(0.35,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position= c(0.85,0.90),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name = "% of predictions") +
  scale_x_discrete(name = "Sample") +
  ggtitle("Overlap with intergenic poly(A) clusters")
dev.off()
