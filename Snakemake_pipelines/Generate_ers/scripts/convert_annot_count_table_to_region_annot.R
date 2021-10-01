#' Converting count table output of derfinder to region annotation
#'
#' \code{convert_annot_count_table_to_region_annot} takes as input the derfinder output and converts to useful region annotation - "intron", "exon", "intergenic"... etc
#'
#' @param count_table count table output from \code{\link{derfinder::annotateRegions}}
#'
#' @return df with the notation "intron", "exon", "intergenic"
#'
#' @export
convert_annot_count_table_to_region_annot <- function(count_table){

  count_table_tmp <-
  count_table %>%
    mutate(exon = ifelse(exon > 0, 1, 0),
           intergenic = ifelse(intergenic > 0, 1, 0),
           intron = ifelse(intron > 0, 1, 0),
           annot_code_tmp = str_c(exon, intergenic, intron)) %>%
    as_tibble()

  count_table[["region_annot"]] <-
    count_table_tmp[["annot_code_tmp"]] %>%
    map(convert_annot_code_to_region_annot) %>%
    unlist()

  count_table_w_region_annot <-
    count_table %>%
    mutate(ER_index = row_number())

  return(count_table_w_region_annot)

}

convert_annot_code_to_region_annot <- function(annot_code){

  region_annot <-
  switch(annot_code,

         `100` = "exon",
         `110` = "exon, intergenic",
         `101` = "exon, intron",
         `111` = "exon, intergenic, intron",
         `010` = "intergenic",
         `011` = "intron, intergenic",
         `001` = "intron",
         `000` = "none"

         )

  return(region_annot)

}
