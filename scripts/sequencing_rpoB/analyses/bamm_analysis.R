# bamm analysis

# Analysis of the Bamm (Bayesian Analysis of Macroevolutionary Mixtures) run looking for variation in speciation, extinction, and diversification rates on the phylogenetic tree.
# follows lots of the work here: http://bamm-project.org/postprocess.html#bammtools

#--------------------------#
# what this script does ####
#--------------------------#

# checks convergence of Bamm run
# looks at the number of rate shifts in the tree
# plots the tree with most likely rate shift configuration
# runs a phylogenetic regression on tip diversification rate across biome preference

#------------------------------#
# load in packages and data ####
#------------------------------#

# load packages
librarian::shelf(caper, ggtree, ggnewscale, RColorBrewer, patchwork, ape, phytools, BAMMtools, coda, MetBrewer, nlme, emmeans, tidyverse, magick)

# read in habitat preference
d_habpref <- read.csv('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv')

# read in phyloseq object and grab tax table
d_taxa <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds') %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
d_meta <- left_join(dplyr::select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), dplyr::select(d_taxa, otu:family))

# load in colours
cols_hab <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds')

# load in tree
#tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre')
tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')
tree2 <- tree
is.ultrametric(tree)


# load in bamm run
mcmcout <- read.csv('data/sequencing_rpoB/bamm/revision/asv/bamm_asv_SF0.5_mcmc_out.txt')

# load in bamm event data
edata <- getEventData(tree, eventdata = 'data/sequencing_rpoB/bamm/revision/asv/bamm_asv_SF0.5_event_data.txt', burnin = 0.3)

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

#-----------------------------------#
# assess convergence of BAMM run ####
#-----------------------------------#

# look at the max generation of the MCMC simulation
max(mcmcout$generation)

# discard some runs as burnin. We will discard the first 30% of samples
burnstart <- floor(0.3 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# calculate effective sample size
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# in general, we want the effective sample size to be at least 200 for these

# visualise prior and posterior simultaneously to look at convergence
d_prior <- plotPrior(mcmcout, expectedNumberOfShifts=500, burnin = 0.3) %>%
  data.frame() %>%
  janitor::clean_names() %>%
  pivot_longer(cols = contains('probs'), names_to = 'type', values_to = 'prob', names_pattern = "(.*)_probs")

#------------------------------------------------------#
# Look at number of rate shifts and model selection ####
#------------------------------------------------------#

# look at probability of different numbers of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs

# get the same thing from the event data
shift_probs <- summary(edata)

# plot the probabilities
plot_rateshifts <- ggplot(shift_probs, aes(shifts, prob)) +
  geom_col(col = 'black', fill = 'light grey') +
  theme_bw(base_size = 14) +
  labs(x = 'Number of shifts',
       y = 'Probability')

# calculate 95% CIs for the number of shifts
n_shifts_ci <- tibble(mean_shifts = mean(postburn$N_shifts),
                      lower_ci = quantile(postburn$N_shifts, 0.025),
                      upper_ci = quantile(postburn$N_shifts, 0.975))
n_shifts_ci
# average number of shifts is 31, lower CI 23, upper CI 40.

# calculate Bayes factors for each number of shifts in rate
# use 0 burnin to sample 0 shifts
mcmc_file = 'data/sequencing_rpoB/bamm/revision/asv/bamm_asv_SF0.5_mcmc_out.txt'
bayes_factors <- computeBayesFactors(mcmc_file, expectedNumberOfShifts=500, burnin=0)

# grab the columns for pairwise comparisons between 0 shifts and number of shifts
d_bayes_factors <- bayes_factors[,1] %>%
  data.frame() %>%
  rownames_to_column(var = 'n_shifts') %>%
  rename(., bayes_factor = `.`)

# we can rank bayes factors and then find the the difference between these
d_bayes_factors <- arrange(d_bayes_factors, desc(bayes_factor)) %>%
  mutate(diff = c(0, abs(diff(bayes_factor))),
         cum_diff = cumsum(diff))

head(d_bayes_factors)
# models withj 28-32 shifts within 20 shifts is the best supported by a long way!
# Bayes factors greater than 20 generally imply strong evidence for one model over another; values greater than 50 are very strong evidence in favour of the numerator model. There is no definitive Bayes factor criterion for “significance”, but many researchers consider values greater than 12 to be consistent with at least some effect.

#---------------------------------#
# Plot results from bamm model ####
#---------------------------------#

# calculate credible shift set - the distinct set of shift configurations that account for 95% of the probability of the data
d_css <- credibleShiftSet(edata, expectedNumberOfShifts = 500, threshold = 5, set.limit = 0.95)

# number of distinct configurations in the data
d_css$number.distinct

# view more information about the credible set
summary(d_css)
# even the single best shift configuration has a very very low posterior probablity

# extract the shift configuration that has the maximum marginal probability and plot
# calculate max shift credibility
msc_set <- maximumShiftCredibility(edata, maximize='product')

# grab the best configuration and plot it
msc_config <- subsetEventData(edata, index = msc_set$sampleindex)
plot.bammdata(msc_config, lwd=2)
addBAMMshifts(msc_config, cex = 2)

# plot this configuration using ggtree

# get mean phylorates that underly the colorised plot produced by plot.bammdata
# from here: https://groups.google.com/g/bamm-project/c/W6s38xzm6OU/m/LALF47xVS54J
#mbt <- getMeanBranchLengthTree(edata, rate = "speciation")

# get the mean branch lengths from the best tree configuration as identified from maximumShiftCredibility 
mbt <- getMeanBranchLengthTree(msc_config, rate = "ndr")

# get shift nodes from "best model"
shiftnodes <- getShiftNodesFromIndex(edata, index = msc_set$sampleindex)

# get tree
tree_bamm <- mbt$phy

tree_bamm$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# get the edge lengths in a dataframe
d_tree_bamm <- data.frame(tree_bamm$edge, edge_num=1:nrow(tree_bamm$edge), edge_length = tree_bamm$edge.length)
colnames(d_tree_bamm)=c("parent", "node", "edge_num", 'edge_length')

# transform these to log for the colour scale
d_tree_bamm <- mutate(d_tree_bamm, log_edge_length = log(edge_length))

# constrained families
constrained_families <- c('Myxococcaceae', 'Vulgatibacteraceae', 'Anaeromyxobacteraceae', 'Polyangiaceae', 'Sandaracinaceae', 'Nannocystaceae', 'Haliangiaceae')

# reorder d_meta so the tip labels link to the order of the tips in the tree
d_meta <- tibble(tip_label = tree_bamm$tip.label) %>%
  left_join(., dplyr::rename(d_meta, tip_label = otu))

# find the mrca of each of the constrained families
d_meta2 <- filter(d_meta, family %in% constrained_families) %>%
  dplyr::select(family, tip_label) %>%
  group_by(family) %>%
  nest() %>%
  mutate(mrca = NA)

for(i in 1:nrow(d_meta2)){
  d_meta2$mrca[i] <- findMRCA(tree_bamm, tips = d_meta2$data[[i]]$tip_label)
}

d_meta2 <- dplyr::select(d_meta2, family2 = family, mrca) %>%
  mutate(blank_label = '')

# add colour for the different families
cols <- c(colorRampPalette(brewer.pal(11, "Spectral"))(nrow(d_meta2)))
names(cols) <- sort(d_meta2$family2)

# remove any tip labels not in our tree from d_meta
d_meta <- filter(d_meta, tip_label %in% tree_bamm$tip.label)

# make a separate column for the rare states to make their size bigger!
d_meta <- mutate(d_meta, rare = ifelse(habitat_preference %in% c('marine mud generalist'), 'rare', 'common'))

# plot tree using ggtree
# first colour branches and add rate shifts
p1 <- ggtree(tree, layout = 'circular', aes(col = log_edge_length)) %<+% d_tree_bamm +
  scale_color_gradientn('Net diversification\n(branch colours)', colors = met.brewer(name='Hiroshige', direction=-1, override.order = F), breaks=c(min(d_tree_bamm$log_edge_length, na.rm = TRUE) + abs(min(d_tree_bamm$edge_length, na.rm = TRUE))*0.2, max(d_tree_bamm$log_edge_length, na.rm = TRUE) * 0.95), labels=c("Slow","Fast")) +
  geom_point2(aes(subset=(node %in% shiftnodes)), color="black",size=5)+
  NULL

# next add tip points
p2 <- p1 %<+% d_meta +
  new_scale_color() +
  geom_tippoint(aes(x=x+x*0.04, col = habitat_preference, size = rare), position = position_jitter(width = 0.025, height = 0)) +
  scale_color_manual('Biome preference\n(tip points)', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  scale_size_manual(values = c(0.6, 3)) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = 'none')

tree_plot <- p2 +
  new_scale_color() +
  geom_cladelab(data = d_meta2,
                mapping = aes(node = mrca,
                              color = family2,
                              label = blank_label),
                offset = castor::get_all_distances_to_root(tree) %>% max() * 0.08,
                barsize = 2) +
  scale_color_manual('Family (outer bar)', values = cols) +
  guides(color = guide_legend(override.aes = list(size = 0.1, shape = 1)))

tree_plot

# save plot out
ggsave('plots/sequencing_rpoB/bamm_tree.png', tree_plot, height = 9, width = 12)

#----------------------------------------#
# look at rate variation through time ####
#----------------------------------------#

# write function to get rate through time into the correct format
get_rate_through_time_df <- function(ephy, ...){
  rtt <- getRateThroughTimeMatrix(ephy, ...)
  
  # get dataframe of each part
  # first speciation
  rtt_sp <- rtt$lambda %>%
    data.frame() %>%
    mutate(sample = 1:n()) %>%
    pivot_longer(starts_with('X'), names_to = 'time_point', values_to = 'speciation') %>%
    mutate(time_point = parse_number(time_point)) %>%
    group_by(sample) %>%
    mutate(time = unname(rtt$times)) %>%
    ungroup()
  
  # second extinction
  rtt_ex <- rtt$mu %>%
    data.frame() %>%
    mutate(sample = 1:n()) %>%
    pivot_longer(starts_with('X'), names_to = 'time_point', values_to = 'extinction') %>%
    mutate(time_point = parse_number(time_point)) %>%
    group_by(sample) %>%
    mutate(time = unname(rtt$times)) %>%
    ungroup()
  
  rtt_comb <- left_join(rtt_sp, rtt_ex)
  
  rtt_comb <- mutate(rtt_comb, net_div = speciation - extinction) %>%
    pivot_longer(cols = c(speciation, extinction, net_div), names_to = 'process', values_to = 'rate')
  
  return(rtt_comb)
}

# get rate through time estimates for the whole tree
rtt_all <- get_rate_through_time_df(ephy = edata)

# create means
rtt_combine_means <- group_by(rtt_all, time_point, process, time) %>%
  summarise(ave_rate = mean(rate), .groups = 'drop',
            lower_ci = quantile(rate, 0.025),
            upper_ci = quantile(rate, 0.975))

# get rate through time estimates for each shift node
rtt_shift <- tibble(shift_node = shiftnodes,
                    n = 1:length(shiftnodes)) %>%
  nest(data = shift_node) %>%
  mutate(temp = purrr::map(data, possibly(~get_rate_through_time_df(ephy = edata, node = .x$shift_node, nodetype = 'include'), otherwise = NA_real_)),
         is_tib = purrr::map_dbl(temp, is_tibble))

rtt_shift2 <- filter(rtt_shift, is_tib == 1) %>%
  dplyr::select(data, temp) %>%
  unnest(data) %>%
  unnest(temp)

# create means
rtt_shift_means <- group_by(rtt_shift2, time_point, process, time, shift_node) %>%
  summarise(ave_rate = mean(rate), .groups = 'drop',
            lower_ci = quantile(rate, 0.025),
            upper_ci = quantile(rate, 0.975))

# create a plot
ggplot() +
  geom_line(aes(time, ave_rate, group = shift_node), rtt_shift_means, col = 'dark grey') +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci, group = shift_node), alpha = 0.01, rtt_shift_means, fill = 'dark grey') +
  geom_line(aes(time, ave_rate), rtt_combine_means) +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci), alpha = 0.1, rtt_combine_means) +
  facet_wrap(~process, scales = 'free') +
  theme_bw(base_size = 14) +
  labs(x = 'Relative time',
       y = 'rate')

# save plot out
ggsave('plots/sequencing_rpoB/bamm_rate_through_time.png', last_plot(), height = 5, width = 12)

# create a plot for just net diversification
p_rtt <- ggplot() +
  geom_line(aes(time, ave_rate, group = shift_node), filter(rtt_shift_means, process == 'net_div'), col = 'dark grey') +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci, group = shift_node), alpha = 0.05, filter(rtt_shift_means, process == 'net_div'), fill = 'dark grey') +
  geom_line(aes(time, ave_rate), filter(rtt_combine_means, process == 'net_div')) +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci), alpha = 0.1, filter(rtt_combine_means, process == 'net_div')) +
  theme_bw(base_size = 10) +
  labs(x = 'Relative time',
       y = 'Net diversification rate')

#--------------------------------------------#
# look at tip-specific evolutionary rates ####
#--------------------------------------------#

# grab out tip rates 
tip_rates <- data.frame(tip_label2 = edata$tip.label,
                        speciation = edata$meanTipLambda,
                        extinction = edata$meanTipMu) %>%
  mutate(tip_label = purrr::map_chr(tip_label2, strsplit_mod)) %>%
  left_join(., dplyr::select(d_meta, tip_label, habitat_preference)) %>%
  filter(!is.na(habitat_preference)) %>%
  group_by(habitat_preference) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(hab_pref_axis = gsub(':', '/ ', habitat_preference),
         hab_pref_axis = gsub('_', ' ', hab_pref_axis),
         net_diversification = speciation - extinction)

# set up correlation matrix for the tree
cor_lambda <- corPagel(value = 1, phy = tree2, form = ~tip_label2)

# fit phylogenetic generalised linear model
mod <- gls(net_diversification ~ habitat_preference, data = tip_rates, correlation = cor_lambda)

# do contrasts between habitat preferences
contrasts <- emmeans(mod, pairwise ~ habitat_preference)

p1 <- ggplot(tip_rates, aes(forcats::fct_reorder(hab_pref_axis, n), speciation)) +
  MicrobioUoE::geom_pretty_boxplot(col='black', fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.2)) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = scales::label_wrap(13)) +
  labs(x = 'Habitat preference',
       y = 'Speciation rate')

p2 <- ggplot(tip_rates, aes(forcats::fct_reorder(hab_pref_axis, n), extinction)) +
  MicrobioUoE::geom_pretty_boxplot(col='black', fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.2)) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = scales::label_wrap(13)) +
  labs(x = 'Habitat preference',
       y = 'Extinction rate')

p3 <- ggplot(tip_rates, aes(forcats::fct_reorder(hab_pref_axis, n), net_diversification)) +
  MicrobioUoE::geom_pretty_boxplot(col='black', fill = 'black') +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.2)) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = scales::label_wrap(13)) +
  labs(x = 'Habitat preference',
       y = 'Net diversification rate')

p1 + p2 + p3

# save plot out
ggsave('plots/sequencing_rpoB/bamm_tip_rates.png', last_plot(), height = 5, width = 17)

# make plot of just net diversification rate
p_tiprates <- mutate(tip_rates, hab_pref_axis = case_when(hab_pref_axis == 'marine mud generalist' ~ 'marine generalist',
                                                          hab_pref_axis == 'freshwater + terrestrial generalist' ~ 'freshwater + land generalist',
                                                          hab_pref_axis == 'marine mud specialist' ~ 'marine specialist',
                                                          hab_pref_axis == 'terrestrial specialist' ~ 'land specialist',
                                                          hab_pref_axis == 'freshwater specialist' ~ 'freshwater specialist')) %>%
  ggplot(., aes(forcats::fct_reorder(hab_pref_axis, n), net_diversification)) +
  MicrobioUoE::geom_pretty_boxplot(aes(col = habitat_preference, fill = habitat_preference), show.legend = FALSE) +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.2), size = 0.8, stroke = 0.6) +
  theme_bw(base_size = 10) +
  scale_x_discrete(labels = scales::label_wrap(13)) +
  labs(x = 'Biome preference',
       y = 'Net diversification rate') +
  scale_color_manual('Biome preference', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  scale_fill_manual('Biome preference', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  ylim(c(0, 10))

#----------------------#
# Assemble Figure 4 ####
#----------------------#

legend <- cowplot::get_legend(tree_plot)

layout <- c(
  'AAAAAB
   CCDDEE'
)

tree_plot2 <- tree_plot + theme(legend.position = 'none')
ggsave('plots/sequencing_rpoB/bamm_tree.png', tree_plot2, height = 8, width = 8)

tree_plot2 <- image_read('plots/sequencing_rpoB/bamm_tree.png', density = 300)
tree_plot2 <- image_trim(tree_plot2)

# make plot
image_ggplot(tree_plot2) + legend + p_rtt + plot_rateshifts + p_tiprates + plot_layout(design = layout, heights = c(0.8, 0.2)) + plot_annotation(tag_levels = list(c('(a)', '', '(b)', '(c)', '(d)'))) &
  theme(plot.tag = element_text(size = 14))

ggsave('plots/manuscript_plots/Figure_3.png', last_plot(), width = 12, height = 12)