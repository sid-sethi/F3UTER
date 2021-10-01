args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(rtracklayer)
})
options(stringsAsFactors=F)


er_file = args[1]
gtf_file <- args[2]
outpath = args[3]
tissue = args[4]



gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
seqlevelsStyle(gtf) <- "UCSC"
gtf.df=as.data.frame(gtf,stringsAsFactor=F)

gtf.gene = gtf.df %>% dplyr::filter(type %in% "gene")
gtf.gene.gr <- makeGRangesFromDataFrame(gtf.gene, keep.extra.columns = TRUE)


er <- read.table(er_file, header=TRUE, sep="\t") %>%
  dplyr::filter(gtf_TxDb_region_annot %in% "intergenic")


##### add nearest any gene  ############
er.gr <- makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

nearest_hit <- nearest(er.gr, gtf.gene.gr,
               select=c("arbitrary"), ignore.strand=FALSE)

er$nearest_any_gene_v94_name <- gtf.gene.gr[nearest_hit]$gene_id


###### calculating associated gene and adding meta data

er_combined <- er %>%
  mutate(associated_gene = nearest_any_gene_v94_name)

er_combined$nearest_any_gene_v94_strand <- dplyr::left_join(er_combined, gtf.gene, by = c("nearest_any_gene_v94_name" = "gene_id")) %>% .$strand.y
er_combined$nearest_any_gene_v94_biotype <- dplyr::left_join(er_combined, gtf.gene, by = c("nearest_any_gene_v94_name" = "gene_id")) %>% .$gene_biotype

# selecting only the required columns in gtf
gtf.gene.small <- gtf.gene %>%
  dplyr::select(seqnames, start, end, gene_id, strand, gene_biotype)


# join the ER details with the gene details to calculate 5' or 3' association of ER
er_combined <- er_combined %>%
  dplyr::left_join(gtf.gene.small, by = c("associated_gene" = "gene_id")) %>%
  mutate(diff_ends = end.y - end.x,
         prime_5_3 = ifelse(diff_ends <= 0, "downstream", "upstream"),
         prime_5_3_corrected = ifelse(strand.y == "-",
                                      ifelse(prime_5_3 == "upstream", 3, 5),
                                      ifelse(prime_5_3 == "upstream", 5, 3)),
         prime_5_3_corrected_chr = ifelse(prime_5_3_corrected == 3, "3'", "5'"),
         distance_from_geneTSS = ifelse(strand.y == "+",
                                        end.x - start.y,
                                        start.x - end.y
                                        ),
         distance_from_geneEnd = ifelse(strand.y == "+",
                                        start.x - end.y,
                                        end.x - start.y
                                        )
         ) %>%
         plyr::rename(c("seqnames.x"="seqnames", "start.x"="start", "end.x"="end", "strand.x"="strand", "strand.y"="associated_gene_strand", "gene_biotype"="associated_gene_biotype")) %>%
         dplyr::select(-c(seqnames.y, start.y, end.y, diff_ends, prime_5_3, prime_5_3_corrected_chr)) %>%
         mutate(ER_id = paste(tissue, seqnames, start, end, strand, sep = ":"))




er.5prime = er_combined %>% dplyr::filter(
                                    prime_5_3_corrected %in% 5,
                                    distance_from_geneTSS < 10000,
                                    distance_from_geneTSS > -10000,
                                    associated_gene_biotype %in% "protein_coding",
                                    width >= 40, width <=2000
                                  )

er.3prime = er_combined %>% dplyr::filter(
                                    prime_5_3_corrected %in% 3,
                                    distance_from_geneEnd < 10000,
                                    distance_from_geneEnd > -10000,
                                    associated_gene_biotype %in% "protein_coding",
                                    width >= 40, width <=2000
                                  )


write.table(er.5prime,
            file = str_c(outpath, "/", tissue, "_5prime_intergenic_ers.txt"),
            row.names = FALSE, quote = FALSE, sep="\t")



write.table(er.3prime,
            file = str_c(outpath, "/", tissue, "_3prime_intergenic_ers.txt", sep=""),
            row.names = FALSE, quote = FALSE, sep="\t")
