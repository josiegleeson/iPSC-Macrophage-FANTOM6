library(stringr)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(rtracklayer)

# import bambu gtf
gr <- rtracklayer::import("Extended_novomix_annotation.gtf")

# import transcripts that had expression over 0.1
filt <- fread("filtered0.1_ids.csv")

# convert granges to df and filter for transcripts that met expression thresholds
df <- gr %>% as_tibble() %>% dplyr::filter(transcript_id %in% filt$transcript_id)

# separate out novel transcripts to rename
novels <- df %>% filter(!str_starts(transcript_id, "ENST"))

# sort and get unique IDs
txids <- data.frame(transcript_id = sort(unique(novels$transcript_id)))

# renumber all transcripts and make them start with BambuTx
txids$new <- paste0("BambuTx", 1:nrow(txids))

# merge dfs
merged_df <- merge(novels, txids, by="transcript_id")

# remove old ID
merged_df$transcript_id <- NULL

# set new ID
merged_df$transcript_id <- df3$new
merged_df$new <- NULL

# write out so that old IDs can be compared to new ones if needed
write.csv(txids, "txids_lookup.csv")

# convert the renamed transcripts df to a granges
newgr <- makeGRangesFromDataFrame(merged_df,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqnames",
                                  start.field="start",
                                  end.field="end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)

# remove scaffolds etc
newgr <- newgr[grep("chr", seqnames(newgr))]

# keep rows without _ or . in the chromosome names
keep_rows <- !grepl("[._]", seqnames(newgr))

# filter the GRanges to apply above
newgr <- newgr[keep_rows]

# filter based on strand
okstrand <- c("+", "-")
newgr <- newgr[strand(newgr) %in% okstrand]

# separate out known transcripts and then convert to new granges
known <- df %>% filter(str_starts(transcript_id, "ENST"))
knowngr <- makeGRangesFromDataFrame(known,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqnames",
                                  start.field="start",
                                  end.field="end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)

# combine the renamed novel transcripts with known ones
comb <- c(newgr, knowngr)

# separate transcript and exon features for sorting
bambu_exons <- comb[comb$type == "exon"]
bambu_transcripts <- comb[comb$type == "transcript"]

# sort by chr and locations
bambu_exons <- sortSeqlevels(bambu_exons)
bambu_transcripts <- sortSeqlevels(bambu_transcripts)
bambu_exons <- sort(bambu_exons)
bambu_transcripts <- sort(bambu_transcripts)

# recombine transcript and exon features
bambu_export <- c(bambu_transcripts, bambu_exons)

# export the final gtf (this was supplied to GenomeProt)
export(bambu_export, "renamed_novomix_annotation.gtf", format="gtf")

# read in transcript counts to rename with new transcript IDs
counts <- fread("bambu_cpm.csv")

# separate known and novel
counts_novel <- counts %>% dplyr::filter(!str_starts(transcript_id, "ENST"))
counts_known <- counts %>% dplyr::filter(str_starts(transcript_id, "ENST"))

# merge new IDs
counts_novel_merged <- merge(counts_novel, txids, by="transcript_id", all.x=T, all.y=F)
counts_novel_merged$transcript_id <- counts_novel_merged$new
counts_novel_merged$new <- NULL

# recombine with known
counts_export <- rbind(counts_novel_merged,counts_known)

# export final counts
write.csv(counts_export, "renamed_bambu_cpm.csv")






