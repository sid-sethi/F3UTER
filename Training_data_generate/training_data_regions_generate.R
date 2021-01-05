library(tidyverse)
library(rtracklayer)

######################################
############ GTF processing ##########
######################################

gtf <- rtracklayer::import("/Annotation/Homo_sapiens.GRCh38.94.gtf")
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)

# Select protein-coding genes and transcripts
gtf.df.pc = gtf.df %>% dplyr::filter(gene_biotype %in% "protein_coding", transcript_biotype %in% "protein_coding")
# Number of protein-coding genes
length(unique(gtf.df.pc$gene_id))
# Number of protein-coding transcripts
length(unique(gtf.df.pc$transcript_id))


# Select TSL level 1
gtf.df.pc$transcript_support_level = gsub("\\s*\\([^\\)]+\\)","",as.numeric(gtf.df.pc$transcript_support_level))
gtf.df.pc.tsl1 = gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)
# Number of protein-coding genes with TSL-1
length(unique(gtf.df.pc.tsl1$gene_id))
# Number of protein-coding transcripts with TSL-1
length(unique(gtf.df.pc.tsl1$transcript_id))


# Removing PAR regions (PAR associated genes retrieved from \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2435358/}{PMID:18660847})
par = read.table("/PAR_genes.txt", header=F, sep="\t")
colnames(par) = c("gene_id", "gene_name")
# removing the PAR genes
gtf.df.pc.tsl1 = dplyr::anti_join(gtf.df.pc.tsl1, par, by = "gene_id")



##############################################################
################# Extract 3 prime UTRs #######################
##############################################################

utr_all = gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "three_prime_utr")

# Collapsing the 3'UTRs among the transcripts for each gene
utr_all_grList = makeGRangesListFromDataFrame(utr_all, split.field = "gene_id", names.field = "transcript_id")

utr_all_collapse = IRanges::reduce(utr_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), three_prime_utr_id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"))

# Total number of 3'UTRs extracted
nrow(utr_all_collapse)
# the plot below shows the distribution of 3'UTRs retrieved per gene
hist(table(utr_all_collapse$group_name))

# removing regions shorter than 40 bp
utr = utr_all_collapse %>% dplyr::filter(width >=40)

# number of UTRs removed using length filter
nrow(utr_all_collapse) - nrow(utr)



#######################################################################
################# Extract Internal coding exons #######################
#######################################################################

# Extract all the coding exons
cds_all = gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "CDS") %>%
  mutate(CDS_id = paste(gene_id,seqnames, start, end, strand, sep=":"))

# number of coding exons per transcript
cds_count = table(cds_all$transcript_id) %>% as.data.frame()

png("3utr_exons_hist.png",type = "cairo", bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
ggplot(cds_count, aes(x=Freq)) +
  geom_histogram(aes(y= (..count..)/sum(..count..)), alpha=0.3, position="identity", lwd=0.3, binwidth=1, color = "blue", fill="blue")+
  theme_bw() +
  #theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(plot.title = element_text(size=11, hjust=0.5),
        #legend.key.size = unit(0.55,"cm"),
        legend.key.size = unit(0.30,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position= c(0.77,0.84),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=9),    #20 before
        axis.text.y = element_text(size=9),
        #panel.border = element_rect(colour="BLACK",size=0.4),
        panel.border = element_blank(),
        axis.title.x = element_text(size=11, vjust = 0.1),
        axis.title.y = element_text(size=11,angle = 90, vjust = 0.6),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  theme(axis.line = element_line(color = "black", size=0.4)) +
  scale_y_continuous(name="% of transcripts (protein-coding + TSL-1)", labels=percent) +
  scale_x_continuous(name = "Coding-exons frequency per transcript")
dev.off()



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
  mutate(elements_collapsed = lengths(revmap), cds_id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"))

# Total number of ICEs extracted
nrow(internal_cds_collapse)

# removing regions shorter than 40 bp
ice = internal_cds_collapse %>% dplyr::filter(width >=40)

# number of ICEs removed using length filter
nrow(internal_cds_collapse) - nrow(ice)

save(utr, ice, file = "/GS.RData")



#######################################################################
################# Extract five prime  #################################
#######################################################################

# extract 5' UTRs
five_prime = gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "five_prime_utr")

# Collapsing the 5'UTRs among the transcripts for each gene
five_prime_grList = makeGRangesListFromDataFrame(five_prime, split.field = "gene_id", names.field = "transcript_id")

five_prime_collapse = IRanges::reduce(five_prime_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), five_prime_utr_id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"))

five = five_prime_collapse %>% dplyr::filter(width >=40)

# saving the data
save(five, file = "/five_prime.RData")



#######################################################################
#################  Extract lncRNAs  ###################################
#######################################################################

lncRNA = c(
  "non_coding",
  "3prime_overlapping_ncRNA",
  "antisense",
  "lincRNA",
  "sense_intronic",
  "sense_overlapping",
  "macro_lncRNA"
)


# Select lncRNA
lncrna.gtf = gtf.df.pc.tsl1 %>% dplyr::filter(gene_biotype %in% lncRNA, type %in% "exon")

# Number of lncRNAs
length(unique(lncrna.gtf$gene_id))
# Number of lncRNA transcripts
length(unique(lncrna.gtf$transcript_id))

# Collapsing the transcripts for each gene
lncrna_all_grList = makeGRangesListFromDataFrame(lncrna.gtf, split.field = "gene_id", names.field = "transcript_id")

lncrna_all_collapse = IRanges::reduce(lncrna_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), lncrna_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep=":"))

# Total number of lncRNAs extracted
nrow(lncrna_all_collapse)

lncrna = lncrna_all_collapse %>% dplyr::filter(width >=40)

# number of lncRNAs removed using length filter
nrow(lncrna_all_collapse) - nrow(lncrna)

# saving the data
save(lncrna, file = "/lncRNA.RData")

# Number of lncRNA exons captured for training (TSL == 1, collapsed, >=40 bp in length)
nrow(lncrna)



####################################################################
#################  Extract ncRNAs  #################################
####################################################################

ncRNA = c(
  "miRNA",
  "misc_RNA",
  "rRNA",
  "snRNA",
  "snoRNA",
  "vaultRNA"
)

# Select ncRNA
# Comments: Transcript support level on ncRNA gene biotype has not been performed by Ensembl i.e. TSL category == NA. Therefore, all transcripts were taken into account.
ncrna.gtf = gtf.df %>% dplyr::filter(gene_biotype %in% ncRNA, type %in% "exon")

# Number of ncRNAs
length(unique(ncrna.gtf$gene_id))
# Number of ncRNA transcripts
length(unique(ncrna.gtf$transcript_id))

# Collapsing among the transcripts for each gene
ncrna_all_grList = makeGRangesListFromDataFrame(ncrna.gtf, split.field = "gene_id", names.field = "transcript_id")

ncrna_all_collapse = IRanges::reduce(ncrna_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), ncrna_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep=":"))

# Total number of ncRNA extracted
nrow(ncrna_all_collapse)

ncrna = ncrna_all_collapse %>% dplyr::filter(width >=40)

# number of lncRNAs removed using length filter
nrow(ncrna_all_collapse) - nrow(ncrna)

# saving the data
save(ncrna, file = "/ncRNA.RData")

# Number of ncRNA exons captured for training (collapsed, >=40 bp in length)
nrow(ncrna)


####################################################################
#################  Extract pseudogenes  #################################
####################################################################

pseudogene = c(
  "pseudogene",
  "processed_pseudogene",
  "unprocessed_pseudogene",
  "transcribed_processed_pseudogene",
  "transcribed_unitary_pseudogene",
  "transcribed_unprocessed_pseudogene",
  "translated_processed_pseudogene",
  "unitary_pseudogene",
  "unprocessed_pseudogene",
  "TR_V_pseudogene",
  "TR_J_pseudogene",
  "rRNA_pseudogene",
  "polymorphic_pseudogene",
  "IG_V_pseudogene",
  "IG_pseudogene",
  "IG_J_pseudogene",
  "IG_C_pseudogene"
)

# Select pseudoGenes
pseudo.gtf = gtf.df.pc.tsl1 %>% dplyr::filter(gene_biotype %in% pseudogene, type %in% "exon")

# Number of pseudoGenes
length(unique(pseudo.gtf$gene_id))
# Number of ncRNA transcripts
length(unique(pseudo.gtf$transcript_id))

# Collapsing among the transcripts for each gene
pseudo_all_grList = makeGRangesListFromDataFrame(pseudo.gtf, split.field = "gene_id", names.field = "transcript_id")

pseudo_all_collapse = IRanges::reduce(pseudo_all_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), pseudoGene_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep=":"))

# Total number of pseudoGenes extracted
nrow(pseudo_all_collapse)

pseudogene = pseudo_all_collapse %>% dplyr::filter(width >=40)

# number of gene removed using length filter
nrow(pseudo_all_collapse) - nrow(pseudogene)

# saving the data
save(pseudogene, file = "/pseudogene.RData")

#Number of pseudoGene exons captured for training (TSL ==1, collapsed, >=40 bp in length)
nrow(pseudogene)
