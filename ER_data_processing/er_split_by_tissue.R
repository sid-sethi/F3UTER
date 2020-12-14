library(dplyr)

#output path
path = "/ERs/ERs_by_tissue"

load("/ERs/ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific.rda") # dataset from Zhang et al.2020

tissues = unique(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific$tissue)

for(i in 1:length(tissues)){
  x = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% dplyr::filter(tissue %in% tissues[i])
  t = tissues[i]
  file = paste(path, "/", t, ".txt", sep="")
  write.table(x, file = file, col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")
}
