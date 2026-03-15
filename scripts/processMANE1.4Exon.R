## Download and process MANE and MANE Clinical file.
## Use R 4.3.3 in Ubuntu; makeTxDbFromGFF in txdbmaker package in current Bioconductor.

library(R.utils)
library(GenomicFeatures)

mane_url = "https://ftp.ncbi.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gff.gz"
download.file(mane_url, destfile = "MANE.GRCh38.v1.4.refseq_genomic.gff.gz")

mane_sum_url = "https://ftp.ncbi.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz"
download.file(mane_sum_url, destfile = "MANE.GRCh38.v1.4.summary.txt.gz")

## Extract exon and add symbol.
mane_db = makeTxDbFromGFF("MANE.GRCh38.v1.4.refseq_genomic.gff.gz")
saveDb(mane_db, file = "mane1.4_db")
mane1.4_exon = exonsBy(mane_db) # in GRlist, "seqnames" for chromosome
mane1.4_exon_df = as.data.frame(mane1.4_exon)
nrow(mane1.4_exon_df) # 204569
mane1.4_exon_df[, 1:2] = NULL # "group"
mane1.4_exon_df$exon_name = gsub("-\\d*", "", mane1.4_exon_df$exon_name)
head(mane1.4_exon_df)

## The text file for gene symbols.
gunzip("MANE.GRCh38.v1.4.summary.txt.gz")
mane_sum = read.table("MANE.GRCh38.v1.4.summary.txt", 
                      header = TRUE, fill = TRUE, quote = "", 
                      comment.char = "", sep = "\t", encoding="UTF-8")
colnames(mane_sum)[1] = "NCBI_GeneID"
mane_sum_simplifed = mane_sum[, c("NCBI_GeneID", "symbol", "RefSeq_nuc", 
                                  "chr_start", "chr_end", "chr_strand")]
mane_sum_simplifed$NCBI_GeneID = gsub("GeneID:", "", 
                                 mane_sum_simplifed$NCBI_GeneID, fixed = TRUE)

## Merge the two tables for gene symbol, extra columns for sanity check.
mane1.4_exon_df = merge(mane1.4_exon_df, mane_sum_simplifed, 
                        by.x = "exon_name", by.y = "RefSeq_nuc")
mane1.4_exon_df = mane1.4_exon_df[, 1:12]
cols = c("seqnames", "start", "end", "strand", "exon_name", "exon_rank", "symbol")
mane1.4_exon_df = mane1.4_exon_df[, cols]

## MANE includes transcripts in *_fix and *_alt chromosomes.
# not duplicate entry in *_fix|alt gene with others, as shown below,
# therefore *_fix|alt was removed from the name.
fa_name = grepl("_fix|alt", mane1.4_exon_df$seqnames)
sum(fa_name) # 732
fa_symbol = mane1.4_exon_df[fa_name, ]$symbol
not_fa_symbol = mane1.4_exon_df[-which(fa_name), ]$symbol
length(fa_symbol) # 732
fa_gene = unique(fa_symbol)
length(fa_gene) # 62
length(not_fa_symbol) # 203837
not_fa_symbol = unique(not_fa_symbol)
length(not_fa_symbol) # 19276
sum(fa_symbol %in% not_fa_symbol) # 0, no duplicate

## Remove *_fix|alt. The file preserves the exon order for a gene. 
mane1.4_exon_df$seqnames = gsub("_.*", "", mane1.4_exon_df$seqnames)
length(unique(mane1.4_exon_df$seqnames)) # 24

colnames(mane1.4_exon_df) = c("chrom", "start", "end", "strand", 
                              "transcript", "exon", "symbol")
write.table(mane1.4_exon_df, "mane1.4_exon", row.names = FALSE, 
            quote = FALSE, sep = "\t")

