args <- commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

options(warn=-1)

gtf_file <- args[1]
outpath <- args[2]
sample <- args[3]
tpm_dir <-args[4]

####### Getting expressed transcripts ######

files <- list.files(path = tpm_dir, pattern = "_abundance.tsv", full.names = TRUE)

tpm_data = files %>%
  map_dfc(read.table, header=TRUE) %>%
  rename_at(1, ~"transcript_id" ) %>%
  dplyr::select(transcript_id, starts_with("tpm"))


expressed_transcripts <- tpm_data %>%
  mutate(rowMeans = rowMeans(.[,-1], na.rm=TRUE)) %>%
  dplyr::filter(rowMeans > 0.1) %>%
  mutate(transcript_id = str_split_fixed(transcript_id, "\\.", n=2)[,1])


######################################
############ GTF processing ##########
######################################

gtf <- rtracklayer::import(gtf_file)
gtf = keepStandardChromosomes(gtf, species = "Drosophila_melanogaster", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)

# expressed gtf - extracting only expressed transcripts
gtf.exp <- dplyr::semi_join(gtf.df, expressed_transcripts, by = "transcript_id")

# Select protein-coding genes and transcripts
gtf.df.pc = gtf.exp %>% dplyr::filter(gene_biotype %in% "protein_coding", transcript_biotype %in% "protein_coding")


##############################################################
################# Extract 3 prime UTRs #######################
##############################################################

utr_all = gtf.df.pc %>% dplyr::filter(type %in% "three_prime_utr")

# Collapsing the 3'UTRs among the transcripts for each gene
utr_all_grList = makeGRangesListFromDataFrame(utr_all, split.field = "gene_id", names.field = "transcript_id")

utr_all_collapse = IRanges::reduce(utr_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"), class = "three_prime")


# removing regions shorter than 40 bp
utr = utr_all_collapse %>% dplyr::filter(width >=40)



#######################################################################
################# Extract Internal coding exons #######################
#######################################################################

# Extract all the coding exons
cds_all = gtf.df.pc %>% dplyr::filter(type %in% "CDS") %>%
  mutate(CDS_id = paste(gene_id,seqnames, start, end, strand, sep=":"))

# number of coding exons per transcript
cds_count = table(cds_all$transcript_id) %>% as.data.frame()


# Only the transcripts with >2 coding exons will contain internal exons, so I remove the transcripts with < 3 coding exons
trans_to_remove = cds_count %>% dplyr::filter(Freq < 3) %>% plyr::rename(c("Var1" = "transcript_id"))
cds.filt = dplyr::anti_join(cds_all,trans_to_remove, by = "transcript_id")
cds.filt$exon_number = as.numeric(cds.filt$exon_number)

# extract INTERNAL CODING EXONS (ICEs): remove the first and the last coding exon
internal.cds = cds.filt %>% dplyr::group_by(transcript_id) %>%
  dplyr::filter(exon_number < max(exon_number), exon_number > min(exon_number)) %>%
  as.data.frame()

# Collapsing the ICEs amongst the transcripts for each gene
internal_cds_grList = makeGRangesListFromDataFrame(internal.cds, split.field = "gene_id", names.field = "transcript_id")

internal_cds_collapse = IRanges::reduce(internal_cds_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"), class = "ICE")

# removing regions shorter than 40 bp
ice = internal_cds_collapse %>% dplyr::filter(width >=40)



#######################################################################
################# Extract five prime  #################################
#######################################################################

# extract 5' UTRs
five_prime = gtf.df.pc %>% dplyr::filter(type %in% "five_prime_utr")

# Collapsing the 5'UTRs among the transcripts for each gene
five_prime_grList = makeGRangesListFromDataFrame(five_prime, split.field = "gene_id", names.field = "transcript_id")

five_prime_collapse = IRanges::reduce(five_prime_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"), class = "five_prime")

five = five_prime_collapse %>% dplyr::filter(width >=40)



####################################################################
#################  Extract ncRNAs  #################################
####################################################################

ncRNA = c(
  "miRNA",
  "rRNA",
  "snRNA",
  "snoRNA",
  "ncRNA"
)

# Select ncRNA
ncrna.gtf = gtf.exp %>% dplyr::filter(gene_biotype %in% ncRNA, type %in% "exon")


# Collapsing among the transcripts for each gene
ncrna_all_grList = makeGRangesListFromDataFrame(ncrna.gtf, split.field = "gene_id", names.field = "transcript_id")

ncrna_all_collapse = IRanges::reduce(ncrna_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep=":"), class = "ncrna")


ncrna = ncrna_all_collapse %>% dplyr::filter(width >=40)



####################################################################
#################  Extract pseudogenes  #################################
####################################################################

pseudogene = c(
  "pseudogene"
)

# Select pseudoGenes
pseudo.gtf = gtf.exp %>% dplyr::filter(gene_biotype %in% pseudogene, type %in% "exon")


# Collapsing among the transcripts for each gene
pseudo_all_grList = makeGRangesListFromDataFrame(pseudo.gtf, split.field = "gene_id", names.field = "transcript_id")

pseudo_all_collapse = IRanges::reduce(pseudo_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep=":"), class = "pseudogene")


pseudogene = pseudo_all_collapse %>% dplyr::filter(width >=40)


data <- rbind(utr, ice, five, ncrna, pseudogene) %>% dplyr::distinct(id, .keep_all=TRUE)
write.table(data, file = str_c(outpath, "/", sample, "_knownRegions.txt"), sep="\t", row.names=FALSE, quote=FALSE)
