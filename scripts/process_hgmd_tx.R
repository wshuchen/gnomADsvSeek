## HGMD transcript
library(txdbmaker)
library(GenomicFeatures)

hgmd = makeFDbPackageFromUCSC(genome = "hg38", tablename = "ncbiRefSeqHgmd")
saveDb(hgmd, "hgmd_db")

# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: UCSC
# Genome: hg38
# Organism: Homo sapiens
# Taxonomy ID: 9606
# UCSC Table: ncbiRefSeqHgmd
# UCSC Track: RefSeq HGMD
# Resource URL: https://genome.ucsc.edu/
# Type of Gene ID: no gene ids
# Full dataset: yes
# Nb of transcripts: 15691
# Db created by: txdbmaker package from Bioconductor
# Creation time: 2026-09-11 21:21:26 -0400 (Fri, 11 Sep 2026)
# txdbmaker version at creation time: 1.8.0
# RSQLite version at creation time: 3.53.3
# DBSCHEMAVERSION: 1.2

## basic stats
# total
mane = loadDb("mane1.4_db")
tx_hgmd = transcripts(hgmd)
length(tx_hgmd) # 15691
tx_mane = transcripts(mane)
length(tx_mane) # 19404

# NR transcript
nr_hgmd = tx_hgmd$tx_name[grepl("NR_", tx_hgmd$tx_name)]
length(nr_hgmd) # 14
nr_mane = tx_mane$tx_name[grepl("NR_", tx_mane$tx_name)]
length(nr_mane) # 50
sum(nr_hgmd %in% nr_mane) # 6

# NM transcript
nm_hgmd = tx_hgmd$tx_name[grepl("NM_", tx_hgmd$tx_name)]
length(nm_hgmd) # 15677
nm_mane = tx_mane$tx_name[grepl("NM_", tx_mane$tx_name)]
length(nm_mane) # 19354
sum(nm_hgmd %in% nm_mane) # 11996
nm_hgmd_not_in_mane = unique(nm_hgmd[!nm_hgmd %in% nm_mane])
length(nm_hgmd_not_in_mane) # 3234

## process the data using functions from recessing the mane transcripts.
db = "hgmd_db"
txdb <- loadDb(db) 
cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
utr3 <- threeUTRsByTranscript(txdb, use.names = TRUE)

exons = exonsBy(txdb, by = "tx", use.names = TRUE)

## Prepare CDS and UTR data and combine them. 
## tx is a transcript ID, e.g., M_001005484.2.
process_feature = function(tx, feature) {
    if (feature == "cds") data = cds
    if (feature == "utr5") data = utr5
    if (feature == "utr3") data = utr3
    
    gr = unlist(data[names(data) == tx])
    mcols(gr) = mcols(gr)[, 2:3]
    mcols(gr)[, 1] = names(gr)[1]
    names(mcols(gr)) = c("transcript", "exon")
    mcols(gr)$feature = feature
    mcols(gr) = mcols(gr)[, c(2:3, 1)]
    gr
}

## Add CDS range in the form of A-B to the data (c. omitted for simplicity).
## tx_x is unlisted transcript data (CDS + UTR).
get_crange = function(tx_x) {
    cumsum_cds = cumsum(width(tx_x)[mcols(tx_x)$feature=="cds"])
    cranges = paste0("1-", cumsum_cds[1])
    if (length(cumsum_cds) > 1) {
        cumsum_cds = c(0, cumsum_cds)
        for (i in 2:(length(cumsum_cds)-1)) {
            cranges = c(cranges, paste0(cumsum_cds[i]+1, "-", cumsum_cds[i+1]))
        }
    }
    cranges
}

process_tx = function(tx) {
    # No all transcripts have UTR.
    if (tx %in% names(utr5) && !tx %in% names(utr3)) {
        features = c("cds", "utr5")
    }
    else if (!tx %in% names(utr5) && tx %in% names(utr3)) {
        features = c("cds", "utr3")
    }
    else if (!tx %in% names(utr5) && !tx %in% names(utr3)) {
        features = "cds"
    } else {
        features = c("cds", "utr5", "utr3")
    }
    
    tx_data = lapply(features, function(x){process_feature(tx, x)})
    tx_data = sort(do.call(c, tx_data))
    if (as.vector(strand(tx_data))[1] == "-") {tx_data = rev(tx_data)}
    mcols(tx_data)$CDS = ifelse(mcols(tx_data)$feature == "cds", "coding",
                                "non-coding")
    cranges = get_crange(tx_data)
    mcols(tx_data)$CDS[mcols(tx_data)$CDS != "non-coding"] = cranges
    mcols(tx_data) = mcols(tx_data)[, c(1, 4, 2, 3)]
    names(tx_data) = NULL
    tx_data
}

process_tx_data = function(tx_names) {
    tx_data = lapply(tx_names, function(x){process_tx(x)})
    tx_data = GRangesList(tx_data)
    names(tx_data) = tx_names
    return(tx_data) # 
}

tx_names = names(cds)
tx_data = process_tx_data(tx_names)

saveRDS(tx_data, "hgmd_transcript_grl.rds")

## Turn the gr list into a data frame. 
tx = as.data.frame(tx_data)
tx$transcript = tx$group_name
tx = tx[, 3:length(tx)] # remove "group", "group_name"
colnames(tx)[1] = "chrom"

## It was easier to prepare a transcript-symbol table from downloaded 
## UCSC refseq-HGMD table (as "refseqHGMD.txt").
hgmd_table = read.table("refseqHGMD.txt", header = TRUE, sep = "\t")
colnames(hgmd_table)
# [1] "bin"          "name"         "chrom"        "strand"       "txStart"      "txEnd"       
# [7] "cdsStart"     "cdsEnd"       "exonCount"    "exonStarts"   "exonEnds"     "score"       
# [13] "name2"        "cdsStartStat" "cdsEndStat"   "exonFrames" 
hgmd_summary = hgmd_table[c("chrom", "txStart", "txEnd", 
                            "cdsStart", "cdsEnd", "strand", "exonCount", "name", "name2")]
colnames(hgmd_summary) = c("chrom", "tx_start", "tx_end", 
                           "cds_start", "cds_end", "strand", "exon_count", "transcript", "symbol")
write.table(hgmd_summary, "hgmd_summary", row.names = FALSE, quote = FALSE)

## Add gene symbol to transcript data frame. 
transcript_symbol = hgmd_summary[, c("transcript", "symbol")]
tx = merge(tx, transcript_symbol, by = "transcript", sort = FALSE)
tx = tx[, c(2:length(tx), 1)]
tx = tx[!duplicated(tx), ]

## Save as RDS. 
saveRDS(tx, "hgmd_transcript.rds")

## Append the HGMD transcripts not in MANE to the MANE table.
mane_tx = readRDS("rdsData/mane1.4_transcript.rds")
length(unique(mane_tx$transcript)) # 19354
hgmd_tx_not_in_mane = hgmd_tx[hgmd_tx$transcript %in% nm_hgmd_not_in_mane, ]
length(unique(hgmd_tx_not_in_mane$transcript)) # 3234
# Found out that the *_alt|fix entries in HGMD are quite messy, better without them.
hgmd_tx_not_in_mane = hgmd_tx_not_in_mane[!grepl("_alt|fix", hgmd_tx_not_in_mane$chrom), ]
length(unique(hgmd_tx_not_in_mane$transcript)) # 3234
mane_hgmd_tx = do.call(rbind, list(mane_tx, hgmd_tx_not_in_mane))
length(unique(mane_hgmd_tx$transcript)) # 22588

## Check _*_alt, _*_fix from chromosome names. 
alt_fix = grepl("_alt|fix", mane_hgmd_tx$chrom)
alt_fix_gene = unique(mane_hgmd_tx[alt_fix, ]$symbol)
length(alt_fix_gene) # 62
not_alt_fix_gene = unique(mane_hgmd_tx[!alt_fix, ]$symbol)
length(not_alt_fix_gene) # 19253
sum(alt_fix_gene %in% not_alt_fix_gene) # 4

alt_fix_tx = unique(mane_hgmd_tx[alt_fix, ]$transcript)
length(alt_fix_tx) # 62
not_alt_fix_tx = unique(mane_hgmd_tx[!alt_fix, ]$transcript)
length(not_alt_fix_tx) # 22526
sum(alt_fix_tx %in% not_alt_fix_tx) # 0 -> no duplicate

saveRDS(mane_hgmd_tx, "mane_hgmd_transcript.rds")












