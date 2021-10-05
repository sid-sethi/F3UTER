args<-commandArgs(TRUE)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
})

options(stringsAsFactors=F)

er_file = args[1]
gtf_file = args[2]
expressed_transcript_file = args[3]
outpath = args[4]
tissue <- args[5]


#make_directory function
make_dir <- function(resPath, name){
  if(!dir.exists(paste(resPath,"/", name, sep=""))){
    system(paste("mkdir -m a=rwx ",resPath, "/", name, sep=""))
  }
}

make_dir(outpath, "5prime")
make_dir(outpath, "3prime")
make_dir(outpath, "Stats")
make_dir(outpath, "5prime_wo_length_filter")
make_dir(outpath, "3prime_wo_length_filter")


gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
seqlevelsStyle(gtf) <- "UCSC"
gtf.df=as.data.frame(gtf,stringsAsFactor=F)

gtf.gene = gtf.df %>% dplyr::filter(type %in% "gene")
gtf.gene.gr <- makeGRangesFromDataFrame(gtf.gene, keep.extra.columns = TRUE)

# expressed gtf - extracting only expressed transcripts
expressed_genes <- read.table(expressed_transcript_file, header=TRUE, sep="\t")

gtf.exp.gr <- dplyr::semi_join(gtf.gene, expressed_genes, by = c("gene_id" = "Name")) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)




er <- read.table(er_file, header=TRUE, sep="\t") %>%
  dplyr::filter(ensembl_grch38_v92_region_annot %in% "intergenic")


##### add nearest any gene and nearest expressed gene ############
er.gr <- makeGRangesFromDataFrame(er, keep.extra.columns = TRUE)

nearest_hit <- nearest(er.gr, gtf.gene.gr,
               select=c("arbitrary"), ignore.strand=FALSE)

er$nearest_any_gene_v94_name <- gtf.gene.gr[nearest_hit]$gene_id

exp_hit <- nearest(er.gr, gtf.exp.gr,
            select=c("arbitrary"), ignore.strand=FALSE)

er$nearest_expressed_gene_v94_name <- gtf.exp.gr[exp_hit]$gene_id



###### calculating associated gene and adding meta data

# get only intergenic regions connected to ONE gene through split read
er_sr <- er %>%
  dplyr::filter(!is.na(uniq_genes_split_read_annot), !str_detect(uniq_genes_split_read_annot, ",")) %>%
  mutate(associated_gene = uniq_genes_split_read_annot)

er_nosr <- er %>%
  dplyr::filter(annotationType_split_read_annot != "partially annotated split read") %>%
  mutate(associated_gene = nearest_expressed_gene_v94_name)

er_combined <- bind_rows(er_sr, er_nosr)


er_combined$nearest_any_gene_v94_strand <- dplyr::left_join(er_combined, gtf.gene, by = c("nearest_any_gene_v94_name" = "gene_id")) %>% .$strand.y
er_combined$nearest_any_gene_v94_biotype <- dplyr::left_join(er_combined, gtf.gene, by = c("nearest_any_gene_v94_name" = "gene_id")) %>% .$gene_biotype

# selecting only the required columns in gtf
gtf.gene.small <- gtf.gene %>%
  dplyr::select(seqnames, start, end, gene_id, strand, gene_biotype)


# join the ER details with the gene details
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
         mutate(ER_id = paste(tissue, seqnames, start, end, strand, sep = ":"),
            is_associated_equal_nearest = ifelse(associated_gene %in% nearest_any_gene_v94_name, "yes", "no"))




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
            file = str_c(outpath, "/5prime/", tissue, ".txt"),
            row.names = FALSE, quote = FALSE, sep="\t")



write.table(er.3prime,
            file = str_c(outpath, "/3prime/", tissue, ".txt", sep=""),
            row.names = FALSE, quote = FALSE, sep="\t")



## ER dataset without length filter
er.5prime_v2 = er_combined %>% dplyr::filter(
                                    prime_5_3_corrected %in% 5,
                                    distance_from_geneTSS < 10000,
                                    distance_from_geneTSS > -10000,
                                    associated_gene_biotype %in% "protein_coding"
                                  )

er.3prime_v2 = er_combined %>% dplyr::filter(
                                    prime_5_3_corrected %in% 3,
                                    distance_from_geneEnd < 10000,
                                    distance_from_geneEnd > -10000,
                                    associated_gene_biotype %in% "protein_coding"
                                  )

write.table(er.5prime_v2,
            file = str_c(outpath, "/5prime_wo_length_filter/", tissue, ".txt"),
            row.names = FALSE, quote = FALSE, sep="\t")


write.table(er.3prime_v2,
            file = str_c(outpath, "/3prime_wo_length_filter/", tissue, ".txt", sep=""),
            row.names = FALSE, quote = FALSE, sep="\t")




all_data <- rbind(er.5prime, er.3prime)

x <- all_data %>% dplyr::filter(is_associated_equal_nearest %in% "yes") %>% nrow()
y <- all_data %>% dplyr::filter(annotationType_split_read_annot != "partially annotated split read", is_associated_equal_nearest %in% "yes") %>% nrow()
z <- all_data %>% dplyr::filter(annotationType_split_read_annot != "partially annotated split read", prime_5_3_corrected %in% 3, is_associated_equal_nearest %in% "yes") %>% nrow()


res <- data.frame(
  type = c("all_data", "all_proximity_data", "3prime_proximity_data"),
  is_associated_gene_also_nearest = c(x, y, z),
  total = c(
    nrow(all_data),
    all_data %>% dplyr::filter(annotationType_split_read_annot != "partially annotated split read") %>% nrow(),
    all_data %>% dplyr::filter(annotationType_split_read_annot != "partially annotated split read", prime_5_3_corrected %in% 3) %>% nrow()
  )
)

res$tissue <- tissue

write.table(res,
            file = str_c(outpath, "/Stats/", tissue, ".txt", sep=""),
            row.names = FALSE, quote = FALSE, sep="\t")
