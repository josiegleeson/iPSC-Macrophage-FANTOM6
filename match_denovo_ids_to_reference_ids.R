

# run gffcompare first:
# gffcompare bambu.gtf -r reference.gtf
# outputs files required: gffcmp.annotated.gtf, gffcmp.tracking

# execute after running gffcompare
library(tidyverse)
library(data.table)
library(reshape2)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# args
combined_gtf <- "gffcmp.annotated.gtf"
tracking_file <- "gffcmp.tracking"
reference_gtf_file <- "gencode.v39.primary_assembly.annotation.gtf"

# import tracking from gffcompare
tracking <- fread(tracking_file, header=F, sep="\t")

# name columns
colnames(tracking) <- c("uniq_id","xloc","ref_id","class_code", "q1")

# separate columns into gene ID and transcript IDs
tracking <- tracking %>% 
  separate(q1, into=c("q_frag", "q1_gene_id", "q1_transcript_id", NULL), sep="[:|]", extra="drop") %>% dplyr::select(-q_frag) %>% 
  separate(ref_id, into=c("ref_gene_id", "ref_transcript_id", NULL), sep="\\|", extra="drop")

# gene lookup database
genedb <- tracking %>% 
  group_by(xloc, q1_gene_id) %>% 
  arrange(desc(ref_gene_id)) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  mutate(ref_gene_id = case_when(
    ref_gene_id == "-" ~ paste0(q1_gene_id),
    TRUE ~ paste0(ref_gene_id)
  )) %>% 
  dplyr::select(xloc, q1_gene_id, ref_gene_id)

# transcript lookup database
txdb <- tracking %>%
  mutate(ref_transcript_id = case_when(
    class_code == "=" ~ paste0(ref_transcript_id),
    class_code == "c" ~ paste0(ref_transcript_id, "_c"),
    TRUE ~ paste0(q1_transcript_id)
  ))

# import combined gtf from gffcompare
gtf <- import(combined_gtf)

# create dataframe from gtf
fromGTFinfo <- data.frame(q1_gene_id = mcols(gtf)$gene_id,
                          q1_transcript_id = mcols(gtf)$transcript_id,
                          xloc = mcols(gtf)$xloc)

# merge gene names
fromGTFinfo_genes <- merge(fromGTFinfo, genedb, by=c("q1_gene_id"), all.x=T, all.y=F)

# arrange by order of gene IDs 
fromGTFinfo_genes <- fromGTFinfo_genes %>%
  arrange(match(q1_gene_id, fromGTFinfo$q1_gene_id))

# add new gene ID to gtf
mcols(gtf)$gene_id <- fromGTFinfo_genes$ref_gene_id

# merge transcript names
fromGTFinfo_txs <- merge(fromGTFinfo, txdb, by="q1_transcript_id", all.x=T, all.y=F)

# arrange by order of transcript IDs
fromGTFinfo_txs <- fromGTFinfo_txs %>%
  arrange(match(q1_transcript_id, fromGTFinfo$q1_transcript_id))

# add new transcript ID to gtf
mcols(gtf)$transcript_id <- fromGTFinfo_txs$ref_transcript_id

# find lines containing novel genes or transcripts and subset 
bambu_index <- grepl("Bambu", mcols(gtf)$transcript_id) | grepl("Bambu", mcols(gtf)$gene_id)
bambu_gr <- gtf[bambu_index]

# function to replace 'Bambu' with 'denovo' in a vector (so that bambu can be re-run)
replace_bambu <- function(x) {
  if (is.character(x)) {
    return(gsub("Bambu", "denovo", x))
  }
  x
}

# apply the replacement function to all metadata columns
mcols(bambu_gr) <- lapply(mcols(bambu_gr), replace_bambu)

# change source column in gtf (so that bambu can be re-run)
bambu_gr$source <- c("denovo")

# import reference GTF
refgtf <- import(reference_gtf_file)

# export combined GTF of novel features and original reference
keep <- c("source", "type", "score", "phase", "transcript_id", "gene_id", "exon_number")

# filter for desired columns
mcols(bambu_gr) <- mcols(bambu_gr)[,keep]

# filter for denovo
condition <- with(mcols(bambu_gr), {
  (grepl("^denovo", gene_id))  
})

# filter for known genes
filtered_bambu_gr <- bambu_gr[!condition]

# filter for desired columns
mcols(refgtf) <- mcols(refgtf)[,keep]

# combine the bambu gtf with the reference gtf
combined <- c(filtered_bambu_gr, refgtf)

# set gtf column order
new_order <- c("source", "type", "score", "phase", "gene_id", "transcript_id",  "exon_number")

# reorder the columns
mcols(combined) <- mcols(combined)[, new_order]

# export the new GTF to supply for bambu
export(combined, "denovo_and_gencode_v39_combined.gtf", format = "gtf")


