#' Will produce a column that indicates T/F whether the gr 1 overlaps with elements in gr 2
#' allows you to ensure if the two grs tested are identical, you are not just picking up identical ranges overlapping
#'
#' @export
mark_overlapping_genes_gr <- function(gr_1, gr_2, identical_gr = F, maxgap = -1L, minoverlap = 1L, ignore.strand = T, ...){

  gr_1_overlapping_counts <-
    gr_1 %>%
    countOverlaps(gr_2, maxgap = maxgap, minoverlap = minoverlap, ignore.strand = ignore.strand, ...)

  if(identical_gr == T){

    stopifnot(all(gr_1_overlapping_counts >= 1))

    gr_1$overlap_gr2 <-
      gr_1_overlapping_counts > 1

  }else{

    gr_1$overlap_gr2 <-
      gr_1_overlapping_counts > 0

  }

  num_overlapping_genes <- sum(gr_1$overlap_gr2)

  print(str_c(num_overlapping_genes, "/", length(gr_1), " (propor: ", round(num_overlapping_genes/length(gr_1), digits = 2), ")", " overlapping ranges.."))

  return(gr_1)

}
