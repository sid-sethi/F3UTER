args <- commandArgs(TRUE)
suppressPackageStartupMessages({
	library("tidyverse")
	library("ODER")
	library("derfinder")
})

bw_dir <- args[1]
gtf <- args[2]
sample <- args[3]
resPath <- args[4]

bw_paths <- list.files(path = bw_dir, pattern = "\\.bw", recursive = TRUE, full.names = TRUE)


bw_df = data.frame()
i = 0
for(f in bw_paths){
  cat(str_c(Sys.time(), ": calculating auc for ", f, "\n"))
  i = i + 1
  auc = derfinder::getTotalMapped(f)
  bw_df[i, "bw_path"] = f
  bw_df[i, "auc"] = auc
}


opt_ers <- ODER(
  bw_paths = bw_df$bw_path, auc_raw = bw_df$auc,
  auc_target = min(bw_df$auc),
  #chrs = c("chr20", "chr21"),
  chrs = c(str_c("chr", seq(1:22)), "chrX", "chrY", "chrM"),
  genome = "hg38", 
  mccs = seq(6, 30, 2), 
  mrgs = seq(20, 140, 10), # 20 - 100
  gtf = gtf, 
  ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, 
  bw_chr = "nochr"
)


x <- plot_ers(opt_ers$deltas, opt_ers$opt_mcc_mrg)
png(str_c(resPath,"/",sample,"_ers_cutoff.png"),bg="transparent",units="in",width = 5.25, height= 6.75 ,res=600)
plot(x)
dev.off()


ers.gr <- opt_ers$opt_ers
mcc <- opt_ers$opt_mcc_mrg[1]
mrg <- opt_ers$opt_mcc_mrg[2]



###### For plugging into local "annotate_ERs" code ########

list_ERs_each_tissue_optimal_cut_off_maxgap <- list()
list_name <- str_c(sample, "-cutoff:", mcc, "-maxgap:", mrg)
list_ERs_each_tissue_optimal_cut_off_maxgap[[1]] <- ers.gr
names(list_ERs_each_tissue_optimal_cut_off_maxgap) <- c(list_name)
  

save(list_ERs_each_tissue_optimal_cut_off_maxgap, file = str_c(resPath, "/", sample, "_optimal_ers.rda"))
  


