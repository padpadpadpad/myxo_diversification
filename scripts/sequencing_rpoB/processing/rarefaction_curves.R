#-------------------------------------------------------#
# Code to make rarefaction curves of myxo sequencing ####
#-------------------------------------------------------#

# load in packages
librarian::shelf(phyloseq, reshape2, tidyverse, parallel)

# custom rarefaction curve
# taken from this GitHub issue: https://github.com/joey711/phyloseq/issues/143
calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=TRUE) {
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- reshape2::melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# set seed
set.seed(42)

# load in dataset
ps <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/clustered/ps_otu_asv.rds'))

# Setting up and registering the cluster
cl <- makeCluster(detectCores(all.tests=TRUE))
registerDoParallel(cl)

to_do <- seq(1, max(sample_sums(ps)), length.out = 100)

# calculate rarefaction curves
# did ten at each number of reads in to-do
rarefaction_curves <- calculate_rarefaction_curves(ps, c('Observed'), rep(to_do, each = 10), parallel=TRUE)

# summarise the rarefaction curves
rarefaction_curves <- group_by(rarefaction_curves, Depth, Sample, Measure) %>%
  summarise(diversity_mean = mean(Alpha_diversity),
            sd = sd(Alpha_diversity),
            .groups = 'drop') %>%
  janitor::clean_names()

# join the rarefaction curves with the metadata for each sample
# read in group clusters from the 16S analysis
clusters <- read.csv('data/sequencing_16s/sample_cluster_assignments.csv') %>%
  mutate(across(where(is.numeric), as.character))

# we will use medoid clustering
clusters <- clusters %>%
  dplyr::select(., sample, medoid_nbclust) %>%
  mutate(clusters, clust = case_when(medoid_nbclust == '1' ~ 'terrestrial',
                                     medoid_nbclust == '2' ~ 'freshwater',
                                     medoid_nbclust == '3' ~ 'marine_mud'))

# from the PCoA plot of the myxobacteria, it can be seen that sample s46 actually is not terrestrial (as it was for the 16s)
# it is instead marine mud and the misidentification likely happened during the DNA extraction / sequencing

# change this cluster assignment
clusters <- mutate(clusters, clust = ifelse(sample == 'sample_s46', 'marine_mud', clust))

# join together
rarefaction_curves <- left_join(rarefaction_curves, select(clusters, sample, clust))

# plot them
ggplot(filter(rarefaction_curves, measure == 'Observed')) +
  geom_line(aes(depth, diversity_mean, group = sample, col = clust), key_glyph = 'point') +
  geom_ribbon(aes(x = depth, ymin = diversity_mean - sd, ymax = diversity_mean + sd, group = sample), alpha = 0.3) +
  theme_bw(base_size = 14) +
  labs(x = 'Sequencing depth (Number of reads)',
       y = 'Estimated species richness') +
  scale_color_manual('Biome', values = c('#5ECAE2', '#3911EE', '#53C20A'), labels = c('Freshwater', 'Marine', 'Land')) +
  theme(legend.position = c(0.85, 0.2)) +
  xlim(c(0, 9e5)) +
  NULL

ggsave('plots/manuscript_plots/rarefaction_rpoB.png', last_plot(), width = 7, height = 5)
