log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)
library(scales)

###########
# GLOBALS #
###########

hist_before_file <- snakemake@input[["hist_before"]]
hist_after_file <- snakemake@input[["hist_after"]]
peak_file <- snakemake@input[["peaks"]]
plot_file <- snakemake@output[["plot"]]

# dev
# hist_before_file <- "output/030_norm/hist.txt"
# hist_after_file <- "output/030_norm/hist_out.txt"
# peak_file <- "output/030_norm/peaks.txt"
# plot_file <- "output/030_norm/kmer_plot.pdf"

########
# MAIN #
########

# read data
before_data <- fread(hist_before_file, colClasses = "numeric")
after_data <- fread(hist_after_file, colClasses = "numeric")
peaks <- fread(paste("grep '^[^#]'", peak_file))

# combine
combined_data <- rbindlist(list(
    Raw = before_data,
    Normalised = after_data),
    idcol = "type")

# arrange plot
combined_data[, type := factor(type, levels = c("Raw", "Normalised"))]

# hlines
mincov <- peaks[which.max(V4), V1]
p1 <- peaks[which.max(V4), V2]
maxcov <- peaks[which.max(V4), V3]

# plot title
gt <- paste0(
    p1, "× 31-mer coverage. ",
    "Main peak: ", mincov, "×–", maxcov, "×"
)

# plot
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
kmer_plot <- ggplot(combined_data, aes(x = `#Depth`, y = Unique_Kmers,
                                       colour = type)) +
    theme_minimal() +
    theme(legend.position = c(5/6, 2/4)) +
    geom_vline(xintercept = c(mincov, p1, maxcov),
               colour = Set1[9]) +
    geom_path(alpha = 0.75) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = NULL)) +
    scale_y_continuous(
        trans = "log10",
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("31-mer depth") + ylab("Number of unique 31-mers") + ggtitle(gt)

# write output
ggsave(filename = plot_file,
       plot = kmer_plot,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
