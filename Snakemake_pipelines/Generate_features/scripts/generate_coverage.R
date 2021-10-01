args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(derfinder)
  library(rtracklayer)
  library(GenomicRanges)
  library(tidyverse)
  library(stringr)
})

options(stringsAsFactors=F)

# Set variables ----------------------------------------------------------------------------------------------

bw_dir <- args[1]
outpath <- args[2]
tissue <- args[3]
chromInfo <- args[4]

bw_paths <- list.files(path = bw_dir, pattern = "\\.bw", recursive = TRUE, full.names = TRUE)


#-------------------------------------------------------------

metadata_df = data.frame()
i = 0
for(f in bw_paths){
  cat(str_c(Sys.time(), ": calculating auc for ", f, "\n"))
  i = i + 1
  auc = derfinder::getTotalMapped(f)
  metadata_df[i, "bigWig_paths"] = f
  metadata_df[i, "auc"] = auc
}


chrominfo_no_scaffold <- read.table(chromInfo, header=FALSE) %>%
    filter(V1 %in% c(1:22, "X", "Y", "MT"))

gtex_tissues <- tissue



# Main ------------------------------------------------------------------------------------------------

for(i in seq_along(gtex_tissues)){

  gtex_tissue_to_filter <- gtex_tissues[i]
  print(str_c(Sys.time(), " - ", i, " - ", gtex_tissue_to_filter))

  for(j in 1:nrow(chrominfo_no_scaffold)){

    chr_to_load <- chrominfo_no_scaffold[["V1"]][j]
    chr_length <- chrominfo_no_scaffold[["V2"]][j]

    print(str_c(Sys.time(), " - loading chromosome: ", chr_to_load))


    # normalise regions to minimum total auc
    tissue_coverage_w_mean_normalised <-
      derfinder::loadCoverage(files = metadata_df$bigWig_paths,
                   totalMapped = metadata_df$auc, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                   targetSize = min(metadata_df$auc), # target by auc
                   chr = chr_to_load,
                   chrlen = chr_length,
                   inputType = "BigWig",
                   returnMean = T,
                   returnCoverage = F,
                   verbose = F,
                   cutoff = NULL) # setting cutoff as null here and instead will be applied downstream in findRegions()

    save(tissue_coverage_w_mean_normalised,
          file = str_c(outpath, "/", gtex_tissue_to_filter, "_", chr_to_load, "_mean_cov.rda"))


    # to save memory after each run
    rm(tissue_coverage_w_mean_normalised)

  }

}
