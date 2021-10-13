# plotting rarefaction curves

# streamlined script for analysis thus far ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(vegan)
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

# set seed
set.seed(42)

# figure path
path_fig <- 'sequencing_16S/plots/analyses'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('sequencing_16S/data/output/run_2/ps_prev_filtered.rds')

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 50947. Woof.
hist(sample_sums(ps))

sort(sample_sums(ps))
# 4th lowest sample has


num_reads <- tibble(reads = sample_sums(ps), sample = names(sample_sums(ps)), x = 'blah')

# remove samples with too few reads
label_reads <- filter(num_reads, reads >= quantile(reads, 0.95) | reads <= quantile(reads, 0.05))

pos = position_jitter(width = 0.2, seed = 1)

ggplot(num_reads, aes(x, reads)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(position = position_jitter(width = 0.2), size = 3, fill = 'white', shape = 21, data = filter(num_reads, ! sample %in% label_reads$sample)) +
  geom_point(position = position_jitter(width = 0.2, seed = 1), size = 3, fill = 'white', shape = 21, data = label_reads) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^0, 10^6),
                minor_breaks = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(x = '',
       y = 'Number of reads',
       title = 'Number of reads per sample') +
  ggrepel::geom_label_repel(aes(label = sample), data = label_reads, position = position_jitter(width = 0.2, seed = 1))

ggsave(file.path(path_fig, 'read_distribution.pdf'), last_plot(), width = 5, height = 9)

# plot differences in sampling effort

# not going to rarefy those samples yet

# can plot rarefaction curves
# code from https://www.fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
# written by Gavin Simpson himself (an author of vegan)

# check rarefaction curves ####
ps_otu_table <- data.frame(otu_table(ps))
raremax <- min(sample_sums(ps))
col <- c('black', 'blue', 'yellow', 'red', 'orange', 'grey', 'hotpink', 'purple', 'green')
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, label = FALSE))

# save plot
pdf(file.path(path_fig, 'rarefaction_curve.pdf'), width = 9, height = 7)
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, label = FALSE))
dev.off()