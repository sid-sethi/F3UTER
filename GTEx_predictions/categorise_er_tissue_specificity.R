args <- commandArgs(TRUE) # 1st  = ML output folder

library(tidyverse)
library(rtracklayer)
options(stringsAsFactors = FALSE)

dir = args[1]

#make_directory function
make_dir <- function(resPath, name){
  if(!dir.exists(paste(resPath,"/", name, sep=""))){
    system(paste("mkdir -m a=rwx ",resPath, "/", name, sep=""))
  }
}

# new output path
make_dir(str_c(dir, "/PredictionMerge"), "Tis_spec")
make_dir(str_c(dir, "/PredictionMerge"), "Brain_spec")
make_dir(str_c(dir, "/PredictionMerge"), "Common")
make_dir(str_c(dir, "/PredictionMerge"), "Ambiguous")

out.path1 = str_c(dir, "/PredictionMerge/Tis_spec")
out.path2 = str_c(dir, "/PredictionMerge/Brain_spec")
out.path3 = str_c(dir, "/PredictionMerge/Common")
out.path4 = str_c(dir, "/PredictionMerge/Ambiguous")


print(str_c(Sys.time(), " - loading ERs from all tissues"))
allErs = readRDS(str_c(dir, "/PredictionMerge_AllTissues/all_tissues_mergedData.rds"))

allErs.gr = makeGRangesFromDataFrame(allErs, keep.extra.columns = TRUE)

tissues = c(
  "adipose_subcutaneous",
  "adrenal_gland",
  "artery_aorta",
  "artery_coronary",
  "artery_tibial",
  "adipose_visceral_omentum",
  "brain_amygdala",
  "brain_anterior_cingulate_cortex_ba24",
  "brain_caudate_basal_ganglia",
  "brain_cerebellar_hemisphere",
  "brain_frontal_cortex_ba9",
  "brain_hippocampus",
  "brain_hypothalamus",
  "brain_nucleus_accumbens_basal_ganglia",
  "brain_putamen_basal_ganglia",
  "brain_spinal_cord_cervical_c_1",
  "brain_substantia_nigra",
  "colon_sigmoid",
  "colon_transverse",
  "esophagus_gastroesophageal_junction",
  "esophagus_mucosa",
  "esophagus_muscularis",
  "heart_atrial_appendage",
  "heart_left_ventricle",
  "kidney_cortex",
  "liver",
  "lung",
  "minor_salivary_gland",
  "muscle_skeletal",
  "nerve_tibial",
  "pancreas",
  "pituitary",
  "skin_not_sun_exposed_suprapubic",
  "skin_sun_exposed_lower_leg",
  "small_intestine_terminal_ileum",
  "spleen",
  "stomach",
  "thyroid",
  "whole_blood"
)


tissueMapping = data.frame(set1 = c("brain_amygdala",
                                    "brain_anterior_cingulate_cortex_ba24",
                                    "brain_caudate_basal_ganglia",
                                    "brain_cerebellar_hemisphere",
                                    "brain_frontal_cortex_ba9",
                                    "brain_hippocampus",
                                    "brain_hypothalamus",
                                    "brain_nucleus_accumbens_basal_ganglia",
                                    "brain_putamen_basal_ganglia",
                                    "brain_spinal_cord_cervical_c_1",
                                    "brain_substantia_nigra"),
                           set2 = c("amygdala",
                                    "anteriorcingulatecortexba24",
                                    "caudatebasalganglia",
                                    "brain_cerebellar_hemisphere",
                                    "frontalcortexba9",
                                    "hippocampus",
                                    "hypothalamus",
                                    "nucleusaccumbensbasalganglia",
                                    "putamenbasalganglia",
                                    "spinalcordcervicalc-1",
                                    "substantianigra"),
                           set3 = c("amygdala",
                                    "anterior_cingulate_cortex_ba24",
                                    "caudate_basal_ganglia",
                                    "cerebellar_hemisphere",
                                    "frontal_cortex_ba9",
                                    "hippocampus",
                                    "hypothalamus",
                                    "nucleus_accumbens_basal_ganglia",
                                    "putamen_basal_ganglia",
                                    "spinal_cord_cervical_c_1",
                                    "substantia_nigra")
                          )



brain_tissues <- tissueMapping$set2


for(tissue in tissues){

  print(str_c(Sys.time(), " - ", tissue))

  # other name for brain tissue
  tis <- ifelse(str_detect(tissue, "brain"), tissueMapping[tissueMapping$set1 %in% tissue, "set2"] , tissue)


  er = readRDS(str_c(dir, "/PredictionMerge/", tissue, "_mergedData.rds"))
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

 # calculating overlap hits
  hits = IRanges::findOverlaps(er.gr, allErs.gr, ignore.strand = TRUE)

  x = er.gr[queryHits(hits)] %>%
    as.data.frame() %>% mutate(index.tmp = 1:n()) %>%
    dplyr::select(seqnames:predicted.prob, tissue, annotationType_split_read_annot, associated_gene, distance_from_geneEnd, index.tmp)
  y = allErs.gr[subjectHits(hits)] %>%
    as.data.frame() %>% mutate(index.tmp = 1:n()) %>%
    dplyr::select(seqnames:predicted.prob, tissue, annotationType_split_read_annot, associated_gene, distance_from_geneEnd, index.tmp)

  xy = dplyr::left_join(x, y, by = "index.tmp")

  er.ovr = xy %>%
    dplyr::group_by(id.x) %>%
    summarise(all_tissues = str_c(unique(tissue.y), collapse = ":"),
            n = n_distinct(tissue.y))

  tis_list <- er.ovr$all_tissues %>%
    str_split(":") %>%
    IRanges::CharacterList()


  er.ovr$brain_count = tis_list %in% brain_tissues %>% sum()
  er.ovr$brain_prop = (er.ovr$brain_count/er.ovr$n)*100 %>% round(2)
  #er.ovr2 <- er.ovr %>% dplyr::select(-tis_list)
  er.ovr2 <- er.ovr

  er.ovr2.data <- dplyr::left_join(er, er.ovr2, by = c("id" = "id.x"))

  # absolute tissue-specific ERs
  tis_spec <- er.ovr2.data %>% dplyr::filter(n == 1)
  saveRDS(tis_spec, file = str_c(out.path1, "/", tissue, "_tisSpec.rds"))

  # brain-specific ERs
  brain_spec <- er.ovr2.data %>% dplyr::filter(n > 1, brain_prop > 75)
  saveRDS(brain_spec, file = str_c(out.path2, "/", tissue, "_brainSpec.rds"))

  # shared ERs
  common.ovr = IRanges::findOverlapPairs(er.gr, allErs.gr, ignore.strand = TRUE, type = "equal", maxgap = 10L) %>%
    as.data.frame() %>%
    dplyr::select(first.X.id, first.X.width, first.X.predicted.prob, first.X.tissue, first.X.annotationType_split_read_annot, first.X.associated_gene, first.X.distance_from_geneTSS, first.X.distance_from_geneEnd, second.X.id, second.X.predicted.prob, second.X.tissue) %>%
    pivot_wider(id_cols = c(first.X.id, first.X.width, first.X.predicted.prob, first.X.tissue, first.X.annotationType_split_read_annot, first.X.associated_gene, first.X.distance_from_geneTSS, first.X.distance_from_geneEnd), names_from = second.X.tissue, values_from = second.X.id) %>%
    dplyr::select(-one_of(tis)) %>%
    as.data.frame()


  # matrix of ER-Ids of shared ERs
  common.id = common.ovr[rowSums(is.na(common.ovr[ , 9:46])) <=33,] %>%
    dplyr::anti_join(brain_spec, by = c("first.X.id" = "id"))# save this
  saveRDS(common.id, file = str_c(out.path3, "/", tissue, "_commonId.rds"))

  # matrix of prediction probabilities of shared ERs
  common.prob = IRanges::findOverlapPairs(er.gr, allErs.gr, ignore.strand = TRUE, type = "equal", maxgap = 10L) %>%
    as.data.frame() %>%
    dplyr::select(first.X.id, first.X.predicted.prob, first.X.tissue, second.X.id, second.X.predicted.prob, second.X.tissue) %>%
    pivot_wider(id_cols = first.X.id, names_from = second.X.tissue, values_from = second.X.predicted.prob) %>%
    as.data.frame()
  saveRDS(common.prob, file = str_c(out.path3, "/", tissue, "_commonProb.rds"))

  # correlation between prediction probabilities of shared ERs
  common.prob.df = common.prob %>% dplyr::select(-first.X.id)
  cor = cor(common.prob.df, use = "pairwise.complete.obs") %>% as.data.frame() %>%
    rownames_to_column("query_tis") %>%
    dplyr::filter(query_tis %in% tis)

  cor.data = rbind(cor.data, cor)


  # ambiguous ERs
  ids <- c(tis_spec$id, brain_spec$id, common.id$first.X.id) %>% unique() %>% as.data.frame()
  ambi.er <- dplyr::anti_join(er, ids, by = c("id" = "."))
  saveRDS(ambi.er, file = str_c(out.path4, "/", tissue, "_ambigs.rds"))



}
