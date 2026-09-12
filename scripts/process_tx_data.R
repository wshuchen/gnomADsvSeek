### Prepare transcript data with CDS and UTR.

library(GenomicFeatures)

db = "mane1.4_db"
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

saveRDS(tx_data, "mane1.4_transcript_grl.rds")

## Turn the gr list into a data frame. 
tx = as.data.frame(tx_data)
tx$transcript = tx$group_name
tx = tx[, 3:length(tx)] # remove "group", "group_name"
colnames(tx)[1] = "chrom"

## Prepare a data frame of refseq summary. 
mane_sum_url = "https://ftp.ncbi.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz"
download.file(mane_sum_url, destfile = "MANE.GRCh38.v1.4.summary.txt.gz")
gunzip("MANE.GRCh38.v1.4.summary.txt.gz")
mane_sum = read.table("MANE.GRCh38.v1.4.summary.txt", 
                      header = TRUE, fill = TRUE, quote = "", 
                      comment.char = "", sep = "\t", encoding="UTF-8")

mane_sum = mane_sum[, -grep("Ensembl", colnames(mane_sum))]
mane_sum = mane_sum[, !colnames(mane_sum) %in% c("HGNC_ID", "name", "MANE_status")]
colnames(mane_sum) = c("GeneID", "symbol", "transcript", 
                       "protein", "chrom", "start", "end", "strand")
mane_sum$GeneID = gsub("GeneID:", "", mane_sum$GeneID, fixed = TRUE)
mane_sum$chrom = gsub(".*\\.", "", mane_sum$chrom)
mane_sum$chrom = paste0("chr", mane_sum$chrom)
mane_sum = mane_sum[mane_sum$transcript %in% tx$transcript, ]

write.table(mane_sum, "mane1.4_summary", row.names = FALSE, quote = FALSE)

## Add gene symbol to transcript data frame. 
transcript_symbol = mane_sum[, c("transcript", "symbol")]

tx = merge(tx, transcript_symbol, by = "transcript", sort = FALSE)
tx = tx[, c(2:length(tx), 1)]

## Save as RDS. 
saveRDS(tx, "mane1.4_transcript.rds")












