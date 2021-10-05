args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

meancov_dir <- args[1]
outpath <- args[2]
min_meancov <- 5

files <- list.files(meancov_dir, pattern = "_meancov.txt", full.names=TRUE)


brain_tissues <- c( "brain_amygdala",
                    "brain_anterior_cingulate_cortex_ba24",
                    "brain_caudate_basal_ganglia",
                    "brain_cerebellar_hemisphere",
                    "brain_frontal_cortex_ba9",
                    "brain_hippocampus",
                    "brain_hypothalamus",
                    "brain_nucleus_accumbens_basal_ganglia",
                    "brain_putamen_basal_ganglia",
                    "brain_spinal_cord_cervical_c_1",
                    "brain_substantia_nigra")


readFiles_c <- function(x){

  tisName = basename(x) %>% stringr::str_replace(., "_meancov.txt", "")

  read.table(x, header=T) %>%
    dplyr::select(-id) %>%
    plyr::rename(c("meanCov" = tisName))
}


readFiles_r <- function(x){

  tisName = basename(x) %>% stringr::str_replace(., "_meancov.txt", "")

  read.table(x, header=T) %>% mutate(tissue = tisName)
}


########### rowwise operations ##########

data_r = files %>%
  map_dfr(readFiles_r)

data_r <- data_r %>% mutate(is_expressed = ifelse(meanCov >= min_meancov, "yes", "no"))

exp.df.long <- data_r %>% dplyr::filter(is_expressed %in% "yes")

exp.group <- exp.df.long %>%
  dplyr::group_by(id) %>%
  summarise(all_tissues = str_c(unique(tissue), collapse = ":"),
          n = n_distinct(tissue),
          .groups = "keep"
        )

tis_list <- exp.group$all_tissues %>%
  str_split(":") %>%
  IRanges::CharacterList()


exp.group$brain_count = tis_list %in% brain_tissues %>% sum()
exp.group$brain_prop = (exp.group$brain_count/exp.group$n)*100 %>% round(2)

tis_spec <- exp.group %>% dplyr::filter(n == 1) %>%
  dplyr::select(id) %>%
  mutate(outcome = "Absolute tissue-specific") %>%
  as.data.frame()

brain_spec <- exp.group %>% dplyr::filter(n > 1, brain_prop > 75) %>%
  dplyr::select(id) %>%
  mutate(outcome = "Highly brain-specific") %>%
  as.data.frame()



########### colwise operations ##########

data_c = files %>%
  map_dfc(readFiles_c)

data_c[data_c < min_meancov] <- "not_expressed"

id = read.table(files[1], header=T) %>%
  dplyr::select(id)

data_id = bind_cols(id, data_c)

all_not_exp <- data_id %>%
  dplyr::filter_at(vars(-id), all_vars(. %in% "not_expressed")) %>%
  dplyr::select(id) %>%
  mutate(outcome = "Not expressed") %>%
  as.data.frame()


data_id$total_exp_tissues = rowSums(!data_id[,-1] == "not_expressed")
shared <- data_id %>% dplyr::filter(total_exp_tissues >= 5) %>%
  dplyr::anti_join(brain_spec, by = "id") %>%
    dplyr::select(id) %>%
    mutate(outcome = "Shared") %>%
    as.data.frame()


outcome = bind_rows(tis_spec, brain_spec, shared, all_not_exp)

res <- dplyr::left_join(id, outcome, by = "id") %>%
  tidyr::replace_na(list(outcome = "Ambiguous"))


write.table(res, file = str_c(outpath, "/utr_tissue_specificity_results.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
