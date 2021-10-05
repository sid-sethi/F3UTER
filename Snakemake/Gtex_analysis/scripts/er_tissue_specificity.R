args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(stringsAsFactors = FALSE, warn = -1)

app_data_merged = args[1]
app_data_dir = args[2]
formatting_csv <- args[3]
outpath = args[4]


#make_directory function
make_dir <- function(resPath, name){
  if(!dir.exists(paste(resPath,"/", name, sep=""))){
    system(paste("mkdir -m a=rwx ",resPath, "/", name, sep=""))
  }
}


# new output path
make_dir(outpath, "Tis_spec")
make_dir(outpath, "Brain_spec")
make_dir(outpath, "Shared")
make_dir(outpath, "Ambiguous")
make_dir(outpath, "ER_tissue_specificity_data")


print(str_c(Sys.time(), " - loading ERs from all tissues"))
allErs = readRDS(app_data_merged)
# remove two tissues - cerebellum and cortex
allErs = allErs %>% filter(! tissue %in% c("brain_cerebellum", "cortex"))
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


cor.data = c()
numbers = data.frame()
tis_spec_df = data.frame()
brain_spec_df = data.frame()
shared_df = data.frame()
ambigs_df = data.frame()

for(tissue in tissues){

  print(str_c(Sys.time(), " - ", tissue))

  # other name for brain tissue
  tis <- ifelse(str_detect(tissue, "brain"), tissueMapping[tissueMapping$set1 %in% tissue, "set2"] , tissue)


  er = readRDS(str_c(app_data_dir, "/", tissue, "_appData.rds"))
  er.gr = makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)


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
  saveRDS(er.ovr2.data, file = str_c(outpath, "/ER_tissue_specificity_data/", tissue, "_res.rds"))


  tis_spec <- er.ovr2.data %>% dplyr::filter(n == 1)
  tis_spec_df <- rbind(tis_spec_df, tis_spec)


  brain_spec <- er.ovr2.data %>% dplyr::filter(n > 1, brain_prop > 75)
  brain_spec_df <- rbind(brain_spec_df, brain_spec)


  #commmon ERs
  common.ovr = IRanges::findOverlapPairs(er.gr, allErs.gr, ignore.strand = TRUE, type = "equal", maxgap = 10L) %>%
    as.data.frame() %>%
    dplyr::select(first.X.id, first.X.seqnames, first.X.start, first.X.end, first.X.strand, first.X.width, first.X.predicted.prob, first.X.tissue, first.X.annotationType_split_read_annot, first.X.associated_gene, first.X.distance_from_geneTSS, first.X.distance_from_geneEnd, second.X.id, second.X.predicted.prob, second.X.tissue) %>%
    pivot_wider(id_cols = c(first.X.id, first.X.seqnames, first.X.start, first.X.end, first.X.strand, first.X.width, first.X.predicted.prob, first.X.tissue, first.X.annotationType_split_read_annot, first.X.associated_gene, first.X.distance_from_geneTSS, first.X.distance_from_geneEnd), names_from = second.X.tissue, values_from = second.X.id) %>%
    dplyr::select(-one_of(tis)) %>%
    as.data.frame()


  common.id = common.ovr[rowSums(is.na(common.ovr[ , 9:46])) <=33,] %>%
    dplyr::anti_join(brain_spec, by = c("first.X.id" = "id"))# save this

  shared_spec <- common.id %>% dplyr::select(first.X.id:first.X.distance_from_geneEnd)
  headers <- colnames(shared_spec)
  colnames(shared_spec) <- str_replace(headers, "first.X.", "")

  shared_df <- rbind(shared_df, shared_spec)

  common.prob = IRanges::findOverlapPairs(er.gr, allErs.gr, ignore.strand = TRUE, type = "equal", maxgap = 10L) %>%
    as.data.frame() %>%
    dplyr::select(first.X.id, first.X.predicted.prob, first.X.tissue, second.X.id, second.X.predicted.prob, second.X.tissue) %>%
    pivot_wider(id_cols = first.X.id, names_from = second.X.tissue, values_from = second.X.predicted.prob) %>%
    as.data.frame()
  #saveRDS(common.prob, file = str_c(outpath, "/Shared/", tissue, "_commonProb.rds"))


  common.prob.df = common.prob %>% dplyr::select(-first.X.id)
  cor = cor(common.prob.df, use = "pairwise.complete.obs") %>% as.data.frame() %>%
    rownames_to_column("query_tis") %>%
    dplyr::filter(query_tis %in% tis)

  cor.data = rbind(cor.data, cor)



  # ambiguous
  ids <- c(tis_spec$id, brain_spec$id, common.id$first.X.id) %>% unique() %>% as.data.frame()
  ambi.er <- dplyr::anti_join(er, ids, by = c("id" = "."))
  ambigs_df <- rbind(ambigs_df, ambi.er)


  # all prob
  breaks = seq(0,1,0.10)

  #### tisSpec
  tis_spec_group = tis_spec %>%
    group_by(range=cut(predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(width), Gene_freq= n_distinct(associated_gene))

  cut_ranges = levels(tis_spec_group$range) %>% as.data.frame()

  tis_spec_allProb = dplyr::left_join(cut_ranges, tis_spec_group, by = c("." = "range")) %>%
    replace_na(list(frequency = 0, length=0, Gene_freq=0)) %>%
    mutate(tissue = tissue, group = "tissue-specific", total_in_group = nrow(tis_spec), perc_group = (frequency/total_in_group)*100) %>%
    plyr::rename( c("." = "range")) %>%
    arrange(as.numeric(range))

  ### brainSpec
  brain_spec_group = brain_spec %>%
    group_by(range=cut(predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(width), Gene_freq= n_distinct(associated_gene))

  cut_ranges = levels(brain_spec_group$range) %>% as.data.frame()

  brain_spec_allProb = dplyr::left_join(cut_ranges, brain_spec_group, by = c("." = "range")) %>%
    replace_na(list(frequency = 0, length=0, Gene_freq=0)) %>%
    mutate(tissue = tissue, group = "highly brain-specific", total_in_group = nrow(brain_spec), perc_group = (frequency/total_in_group)*100) %>%
    plyr::rename( c("." = "range")) %>%
    arrange(as.numeric(range))


  ### shared
  shared_group = common.id %>%
    group_by(range=cut(first.X.predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(first.X.width), Gene_freq= n_distinct(first.X.associated_gene))

  cut_ranges = levels(shared_group$range) %>% as.data.frame()

  shared_allProb = dplyr::left_join(cut_ranges, shared_group, by = c("." = "range")) %>%
    replace_na(list(frequency = 0, length=0, Gene_freq=0)) %>%
    mutate(tissue = tissue, group = "shared", total_in_group = nrow(common.id), perc_group = (frequency/total_in_group)*100) %>%
    plyr::rename( c("." = "range")) %>%
    arrange(as.numeric(range))


  #### ambiguous
  ambi_group = ambi.er %>%
    group_by(range=cut(predicted.prob, breaks= seq(0, 1, by = 0.10)) ) %>%
    summarise(frequency= n(), length = sum(width), Gene_freq= n_distinct(associated_gene))

  cut_ranges = levels(ambi_group$range) %>% as.data.frame()

  ambi_allProb = dplyr::left_join(cut_ranges, ambi_group, by = c("." = "range")) %>%
    replace_na(list(frequency = 0, length=0, Gene_freq=0)) %>%
    mutate(tissue = tissue, group = "ambiguous", total_in_group = nrow(ambi.er), perc_group = (frequency/total_in_group)*100) %>%
    plyr::rename( c("." = "range")) %>%
    arrange(as.numeric(range))


  numbers = rbind(tis_spec_allProb, brain_spec_allProb, shared_allProb, ambi_allProb) %>%
    mutate(total_in_tissue = nrow(er), perc_tissue = (frequency/total_in_tissue)*100) %>%
    rbind(numbers)



  rm(tis_spec_allProb, brain_spec_allProb, shared_allProb, ambi_allProb)
  rm(tis_spec, brain_spec, ambi.er, common.id)

}


# saving data
saveRDS(tis_spec_df, file = str_c(outpath, "/Tis_spec/all_tisSpec.rds"))
saveRDS(brain_spec_df, file = str_c(outpath, "/Brain_spec/all_brainSpec.rds"))
saveRDS(shared_df, file = str_c(outpath, "/Shared/all_sharedId.rds"))
saveRDS(ambigs_df, file = str_c(outpath, "/Ambiguous/all_ambigs.rds"))


formatting <- read_delim(formatting_csv, delim = ",") %>% mutate(tissue_color_hex = str_c("#", tissue_color_hex)) %>%
  drop_na() %>%
  dplyr::select(gtex_tissues_name_formatted_2, gtex_tissue_group, gtex_tissues_name_to_plot, tissue_color_hex, OMIM_gtex_name)


# function to fix tissue names
fix_tissue_names <- function(x){
  colnames_new <- c()
  for(name in x){
    if(name %in% formatting$OMIM_gtex_name){
      new_name = formatting[formatting$OMIM_gtex_name %in% name, "gtex_tissues_name_to_plot"]
    } else{
      new_name = name
    }
    colnames_new <- c(colnames_new, new_name)
  }
  return(colnames_new)
}


colnames = colnames(cor.data)
rownames = rownames(cor.data)

colnames(cor.data) <- fix_tissue_names(colnames)
cor.data$query_tis <- fix_tissue_names(cor.data$query_tis)

saveRDS(cor.data, file = str_c(outpath, "/Shared/all_corData.rds"))


numbers$tissue <- str_replace(numbers$tissue, "brain_", "")
saveRDS(numbers, file = str_c(outpath, "/all_numbers.rds"))
