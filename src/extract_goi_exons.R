
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(rtracklayer)

# catch data
gff_file <- snakemake@input[["gff"]]
transcript_id <- snakemake@params[["goi"]]
goi_exons_file <- snakemake@output[["goi_exons"]]

# read gff and subset
gff <- import.gff(gff_file, feature.type = "exon")
goi_exons <- gff[gff$transcript_id %in% c(transcript_id)]

# write a csv of positions
goi <- as.data.table(goi_exons)[, .(seqnames, start, end)]
fwrite(goi, goi_exons_file)

# log
sessionInfo()
