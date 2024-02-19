# look at bamm analysis across the five sampling fractions for ASV

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
tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')
tree2 <- tree
is.ultrametric(tree)

tree_97.7 <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_97.7/myxo_97.7_treepl_cv_node_labels.tre')

tree_95 <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_95/myxo_95_treepl_cv_node_labels.tre')

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

#----------------------------------#
# check convergence of each run ####
#----------------------------------#

# list all mcmcout files
mcmc_files <- list.files('data/sequencing_rpoB/bamm/revision/asv', pattern = 'mcmc', full.names = TRUE)
mcmc_files_97.7 <- list.files('data/sequencing_rpoB/bamm/revision/myxo97.7', pattern = 'mcmc', full.names = TRUE)
mcmc_files_95 <- list.files('data/sequencing_rpoB/bamm/revision/myxo95', pattern = 'mcmc', full.names = TRUE)

# write a function to check convergence
check_convergence <- function(mcmc_file, burnin = 0.25){
  # read in file
  temp <- read.csv(mcmc_file)
  
  burnstart <- floor(burnin * nrow(temp))
  postburn <- temp[burnstart:nrow(temp), ]
  
  # calculate effective sample size
  temp <- tibble(neff_nshifts = coda::effectiveSize(postburn$N_shifts),
         neff_loglik = coda::effectiveSize(postburn$logLik))
  
  return(temp)
}

# run function on mcmcfiles and purrr
convergence_asv <- purrr::map_dfr(mcmc_files, check_convergence, burnin = 0.25, .id = 'run') %>%
  mutate(file = basename(mcmc_files),
         samp_frac = parse_number(file))
# pretty good, use a burn in of 0.25 throughout

convergence_asv_97.7 <- purrr::map_dfr(mcmc_files_97.7, check_convergence, burnin = 0.6, .id = 'run') %>%
  mutate(file = basename(mcmc_files_97.7),
         samp_frac = parse_number(file))

convergence_asv_95 <- purrr::map_dfr(mcmc_files_95, check_convergence, burnin = 0.15, .id = 'run') %>%
  mutate(file = basename(mcmc_files_95),
         samp_frac = parse_number(file))

#-------------------------------------------------------#
# look at number of shifts at each sampling fraction ####
#-------------------------------------------------------#

get_shift_cis <- function(mcmc_file, tree, burnin = 0.25){
  
  # read in file
  temp <- read.csv(mcmc_file)
  
  burnstart <- floor(burnin * nrow(temp))
  postburn <- temp[burnstart:nrow(temp), ]
  
  # calculate 95% CIs for the number of shifts
  n_shifts_ci <- tibble(mean_shifts = mean(postburn$N_shifts),
                        lower_ci = quantile(postburn$N_shifts, 0.025),
                        upper_ci = quantile(postburn$N_shifts, 0.975))
  return(n_shifts_ci)
}

n_shifts_ci_asv <- purrr::map_dfr(mcmc_files, get_shift_cis, tree, burnin = 0.25, .id = 'run') %>%
  mutate(file = basename(mcmc_files),
         samp_frac = parse_number(file),
         type = 'ASV')

n_shifts_ci_97.7 <- purrr::map_dfr(mcmc_files_97.7, get_shift_cis, tree, burnin = 0.25, .id = 'run') %>%
  mutate(file = basename(mcmc_files_97.7),
         samp_frac = gsub('97.7', '', file) |> parse_number(),
         type = '97.7% OTU')

n_shifts_ci_95 <- purrr::map_dfr(mcmc_files_95, get_shift_cis, tree, burnin = 0.25, .id = 'run') %>%
  mutate(file = basename(mcmc_files_95),
         samp_frac = gsub('95', '', file) |> parse_number(),
         type = '95% OTU')

# combine these together
n_shifts_ci <- bind_rows(n_shifts_ci_asv, n_shifts_ci_97.7, n_shifts_ci_95) %>%
  mutate(samp_frac = as.character(samp_frac))

# create plot to look at number of shifts across levels of sampling fraction
ggplot(aes(x = samp_frac, y = mean_shifts, ymin = lower_ci, ymax = upper_ci), data = n_shifts_ci) +
  geom_point(size = 3) +
  geom_linerange() +
  theme_bw(base_size = 14) +
  labs(x = 'Sampling fraction', y = 'Number of shifts') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~type, labeller = labeller(type = MicrobioUoE::letter_facets))

# save the plot out
ggsave('plots/sequencing_rpoB/bamm_shifts_CIs.png', width = 10, height = 4)

#-------------------------------------------------#
# grab shift nodes from each sampling fraction ####
#-------------------------------------------------#

# use the maximum marginal probability set

# write a function
get_shift_nodes <- function(event_file, tree, burnin = 0.25){
  
  # read in file
  temp <- getEventData(tree, eventdata = event_file, burnin = burnin)
  
  # extract the shift configuration that has the maximum marginal probability and plot
  # calculate max shift credibility
  msc_set <- maximumShiftCredibility(temp, maximize='product')
  
  # get shift nodes from "best model"
  shiftnodes <- getShiftNodesFromIndex(temp, index = msc_set$sampleindex)
  
  return(tibble(shiftsnodes = shiftnodes))
}

# list event files
event_files_asv <- list.files('data/sequencing_rpoB/bamm/revision/asv', pattern = 'event', full.names = TRUE)
event_files_97.7 <- list.files('data/sequencing_rpoB/bamm/revision/myxo97.7', pattern = 'event', full.names = TRUE)
event_files_95 <- list.files('data/sequencing_rpoB/bamm/revision/myxo95', pattern = 'event', full.names = TRUE)

# get shift nodes for each sampling fraction
shift_nodes_asv <- tibble(file = basename(event_files_asv),
                          samp_frac = parse_number(file)) %>%
  mutate(shiftsnodes = purrr::map(event_files_asv, get_shift_nodes, tree2, burnin = 0.25)) %>%
  unnest(cols = shiftsnodes) %>%
  count(shiftsnodes) %>%
  arrange(desc(n))

shift_nodes_97.7 <- tibble(file = basename(event_files_97.7),
                          samp_frac = parse_number(file)) %>%
  mutate(shiftsnodes = purrr::map(event_files_97.7, get_shift_nodes, tree_97.7, burnin = 0.25)) %>%
  unnest(cols = shiftsnodes) %>%
  count(shiftsnodes) %>%
  arrange(desc(n))

shift_nodes_95 <- tibble(file = basename(event_files_95),
                          samp_frac = parse_number(file)) %>%
  mutate(shiftsnodes = purrr::map(event_files_95, get_shift_nodes, tree_95, burnin = 0.15)) %>%
  unnest(cols = shiftsnodes) %>%
  count(shiftsnodes) %>%
  arrange(desc(n))

#-------------------------------------------------------#
# make plot of shift nodes across sampling fractions ####
#-------------------------------------------------------#

make_phylorate_panel_plot <- function(event_files, tree, burnin = 0.25){
  
  temp_file_names <- tempfile(letters[1:length(event_files)])
  temp_file_names <- paste(temp_file_names, '.png', sep = '')
  
  for(i in 1:length(event_files)){
    
    temp <- getEventData(tree, eventdata = event_files[i], burnin = burnin)
    
    # extract the shift configuration that has the maximum marginal probability and plot
    # calculate max shift credibility
    msc_set <- maximumShiftCredibility(temp, maximize='product')
    
    # grab the best configuration and plot it
    msc_config <- subsetEventData(temp, index = msc_set$sampleindex)
    
    # get shift nodes from "best model"
    shiftnodes <- getShiftNodesFromIndex(temp, index = msc_set$sampleindex)
    
    # get the mean branch lengths from the best tree configuration as identified from maximumShiftCredibility 
    mbt <- getMeanBranchLengthTree(msc_config, rate = "ndr")
    
    # get tree
    tree_bamm <- mbt$phy
    
    # get the edge lengths in a dataframe
    d_tree_bamm <- data.frame(tree_bamm$edge, edge_num=1:nrow(tree_bamm$edge), edge_length = tree_bamm$edge.length)
    colnames(d_tree_bamm)=c("parent", "node", "edge_num", 'edge_length')
    
    # transform these to log for the colour scale
    d_tree_bamm <- mutate(d_tree_bamm, log_edge_length = log(edge_length))
    
    # plot tree
    p1 <- ggtree(tree, layout = 'circular', col = 'light grey') %<+% d_tree_bamm +
      geom_point2(aes(subset=(node %in% shiftnodes)), color="black",size=3) +
      theme(legend.position = 'none') +
      NULL
    
    ggsave(temp_file_names[i], p1, width = 6, height = 6)
    
    }
  
  p1 <- image_trim(image_read(temp_file_names[1]))
  p2 <- image_trim(image_read(temp_file_names[2]))
  p3 <- image_trim(image_read(temp_file_names[3]))
  p4 <- image_trim(image_read(temp_file_names[4]))
  p5 <- image_trim(image_read(temp_file_names[5]))
    
  image_ggplot(p1) + image_ggplot(p2) + image_ggplot(p3) + image_ggplot(p4) + image_ggplot(p5) + plot_annotation(tag_levels = list(c('(a) SF = 0.0625', '(b) SF = 0.125', '(c) SF = 0.25', '(d) SF = 0.5', '(e) SF = 1'))) &
    theme(plot.tag = element_text(size = 14),
          plot.tag.position = 'top')
}

make_phylorate_panel_plot(event_files_asv, tree2, burnin = 0.25)
ggsave('plots/sequencing_rpoB/bamm_asv_sampfracs.png', width = 10, height = 8)

make_phylorate_panel_plot(event_files_97.7, tree_97.7, burnin = 0.25)
ggsave('plots/sequencing_rpoB/bamm_97.7_sampfracs.png', width = 10, height = 8)

make_phylorate_panel_plot(event_files_95, tree_95, burnin = 0.25)
ggsave('plots/sequencing_rpoB/bamm_95_sampfracs.png', width = 10, height = 8)
