#' Generating a genomic state object from Txdb
#'
#' \code{generate_genomic_state} takes txdb object and makes a genomic state
#'
#' @param ensembl_grch38_v87_TxDb Txdb object
#' @param output_path path to output - if it already exists, then function will simply load
#' @param chrs_to_keep chromosomes to keep in genomic state
#'
#' @return genomicstate
#'
#' @export
generate_genomic_state <- function(ensembl_grch38_v87_TxDb, output_path, chrs_to_keep = str_c("chr", c(1:22, "X", "Y", "M"))){

  if(!file.exists(output_path)){

    genomic_state_ensembl_grch38_v87 <-
      makeGenomicState(txdb = ensembl_grch38_v87_TxDb, chrs = chrs_to_keep)

    save(genomic_state_ensembl_grch38_v87,
         file = output_path)

  }else {

    load(file = output_path)

  }

  return(genomic_state_ensembl_grch38_v87)

}
