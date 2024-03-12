#---------------------------------------------------------#
# script to run models of discrete character evolution ####
#---------------------------------------------------------#

# what this script does ####

# 1. Fits Mk models for discrete character evolution
# 2. Does model selection on these models
# 3. does bootstraps on the best model sampling 80% of tips and refitting the best model

# load packages ####
librarian::shelf(here, tidyverse, ggtree, ggnewscale, RColorBrewer, patchwork, phytools, diversitree, tidygraph, igraph, ggraph, GGally, ggrepel, flextable, padpadpadpad/MicrobioUoE, cli)

# read in datasets ####

# read in habitat preference
d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv'))

# read in phyloseq object and grab tax table
d_taxa <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds')) %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')

# create d_meta
# use habitat_preference3 which collapses into marine mud generalists
d_meta <- left_join(select(d_habpref, otu, habitat_preference = habitat_preference3, num_present), select(d_taxa, otu:family))

# read in tree
tree <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

tree <- force.ultrametric(tree)

# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

source('scripts/sequencing_rpoB/analyses/diversitree_helper_functions.R')

# setup for analyses ####

# reorder metadata to match tip labels of tree
d_meta <- tibble(tip_label = tree$tip.label) %>%
  left_join(., rename(d_meta, tip_label = otu))

# make sure order of habitat preference is the same in the trait vector as the tip labels
sum(d_meta$tip_label == tree$tip.label) == length(tree$tip.label)
# SUCCESS if TRUE

# use habitat preference 3 - which collapses marine mud generalists into a single state
unique(d_meta$habitat_preference)

hab_pref <- setNames(d_meta$habitat_preference, d_meta$tip_label)
hab_pref_num <- as.numeric(as.factor(hab_pref))
hab_pref_num <- setNames(hab_pref_num, d_meta$tip_label)

# coding from numeric to character
coding <- tibble(hab_pref = unname(hab_pref), hab_pref_num = unname(hab_pref_num)) %>%
  distinct() %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub(':', '/ ', hab_pref),
  hab_pref_axis = gsub('_', ' ', hab_pref_axis),
  # rename the columns for easy naming
  initials = c('FG', 'FS', 'MG', 'MS', 'TS'))

coding

# save out coding
write.csv(coding, 'data/sequencing_rpoB/processed/transition_rates/coding.csv', row.names = FALSE)

#-------------------------------------------------------------#
# run diversitree analysis of discrete character evolution ####
#-------------------------------------------------------------#

# all of these methods needs a likelihood function, we can build a Mkn model
lik_ard <- make.mkn(tree, hab_pref_num, k = max(hab_pref_num))

argnames(lik_ard)

# make symmetric rates model
lik_sym <- constrain(lik_ard, 
                     q12~q21, q13~q31, q14~q41, q15~q51,
                     q23~q32, q24~q42, q25~q52,
                     q34~q43, q35~q53,
                     q45~q54)

argnames(lik_sym)

# make equal rates model
lik_er <- constrain(lik_ard,
                    q13~q12, q14~q12, q15~q12,
                    q21~q12, q23~q12, q24~q12, q25~q12,
                    q31~q12, q32~q12, q34~q12, q35~q12,
                    q41~q12, q42~q12, q43~q12, q45~q12,
                    q51~q12, q52~q12, q53~q12, q54~q12)

argnames(lik_er)

lik_transient <- constrain(lik_ard,
                           q24~0, q25~0, q42~0, q52~0, q14~0, q41~0, q45~0, q54~0)

# need to pass start values to it - can grab these from the ape::ace, but we will just pass an average rate to the model
inits_ard <- rep(1, length(argnames(lik_ard)))
inits_sym <- rep(1, length(argnames(lik_sym)))
inits_er <- rep(1, length(argnames(lik_er)))
inits_trans <- rep(1, length(argnames(lik_transient)))

# run first four models
mod_er <- find.mle(lik_er, inits_er, method = 'subplex', control = list(maxit = 50000))
mod_ard <- find.mle(lik_ard, inits_ard, method = 'subplex', control = list(maxit = 50000))
mod_trans <- find.mle(lik_transient, inits_trans, method = 'subplex', control = list(maxit = 50000))
mod_sym <- find.mle(lik_sym, inits_sym, method = 'subplex', control = list(maxit = 50000))

# compare models
AIC(mod_sym, mod_ard, mod_er, mod_trans) %>% arrange(AIC)

# sort parameter estimates
mod_ard$par %>% sort()

# get parameters close to 0.
mod_ard$par[mod_ard$par < 1e-03] %>% names(.)

# make custom matrix model
lik_custom1 <- constrain(lik_ard, 
                     q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0)

# make start parameters
inits_custom1 <- rep(1, length(argnames(lik_custom1)))

# refit model
mod_custom1 <- find.mle(lik_custom1, inits_custom1, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1) %>% arrange(AIC)

# filter for only estimated parameters
mod_custom1$par.full[names(mod_custom1$par.full) %in% argnames(lik_custom1)] 

# sort parameter estimates
mod_custom1$par %>% sort()

# make custom matrix model again
lik_custom2 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0)

# make start parameters
inits_custom2 <- rep(1, length(argnames(lik_custom2)))

# run model
mod_custom2 <- find.mle(lik_custom2, inits_custom2, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2) %>% arrange(AIC)

# sort parameter estimates
mod_custom2$par %>% sort()

# make custom matrix model
lik_custom3 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0,
                         q54~0)

# make start parameters
inits_custom3 <- rep(1, length(argnames(lik_custom3)))

# run model
mod_custom3 <- find.mle(lik_custom3, inits_custom3, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3) %>% arrange(AIC)

# sort parameter estimates
mod_custom3$par %>% sort()

# make custom matrix model
lik_custom4 <- constrain(lik_ard, 
                         q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                         q45~0,
                         q54~0,
                         q43~0)

# make start parameters
inits_custom4 <- rep(1, length(argnames(lik_custom4)))

# run model
mod_custom4 <- find.mle(lik_custom4, inits_custom4, method = 'subplex', control = list(maxit = 50000))

# do AIC comparison
AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>% arrange(AIC)

# sort parameter estimates
mod_custom4$par %>% sort()

# lets look at model weights
d_aic <- AIC(mod_sym, mod_ard, mod_er, mod_trans, mod_custom1, mod_custom2, mod_custom3, mod_custom4) %>%
  data.frame() %>%
  mutate(log_lik = c(mod_sym$lnLik, mod_ard$lnLik, mod_er$lnLik, mod_trans$lnLik, mod_custom1$lnLik, mod_custom2$lnLik, mod_custom3$lnLik, mod_custom4$lnLik),
         ntips = 2621) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = 'model') %>%
  mutate(aicc = -2*log_lik + 2*df*(ntips/(ntips - df - 1)),
         aic_weight = round(MuMIn::Weights(aic), 2),
         aicc_weight = round(MuMIn::Weights(aicc), 2)) %>%
  arrange(aic)
# mod custom 3 is favoured pretty highly

# make this into a table
d_table <- select(d_aic, model, df, log_lik, aic, aic_weight) %>%
  mutate(model = case_when(model == 'mod_er' ~ 'ER',
                           model == 'mod_ard' ~ 'ARD',
                           model == 'mod_sym' ~ 'SYM',
                           model == 'mod_trans' ~ 'Stepwise',
                           model == 'mod_custom1' ~ 'Simplified ARD 1',
                           model == 'mod_custom2' ~ 'Simplified ARD 2',
                           model == 'mod_custom3' ~ 'Simplified ARD 3',
                           model == 'mod_custom4' ~ 'Simplified ARD 4'))

# make table
table <- d_table %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    df = 'd.f.',
                    loglik = 'Log Likelihood',
                    aic = 'AIC',
                    aic_weight = "AIC weight") %>%
  italic(j = 2, part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

save_as_image(table, 'plots/manuscript_plots/markov_model_table.png')
save_as_docx(table, path = 'plots/manuscript_plots/markov_model_table.docx')

# save out best markov model
saveRDS(mod_custom3, 'data/sequencing_rpoB/processed/transition_rates/mod_custom_3.rds')

# plot_transition_matrix
diversitree_df <- get_diversitree_df(mod_custom3, coding$hab_pref_num, coding$hab_pref)

diversitree_df %>%
  left_join(., select(coding, state_1 = hab_pref, state_1_num = hab_pref_num, state_1_label = hab_pref_axis)) %>%
  left_join(., select(coding, state_2 = hab_pref, state_2_num = hab_pref_num, state_2_label = hab_pref_axis)) %>%
  mutate(transition_rate = round(transition_rate, 2)) %>%
  ggplot(., aes(forcats::fct_reorder(state_2_label, state_2_num), forcats::fct_reorder(state_1_label, desc(state_1_num)))) +
  geom_tile(aes(alpha = transition_rate, col = free_param), width = 0.9, height = 0.9, size = 1.1) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
  legend.position = 'none',
  axis.text.x.top = element_text(angle = 90, vjust = 0.5),
  plot.title.position = "plot") +
  scale_alpha_continuous(range = c(0, 0.6)) +
  geom_text(aes(label = transition_rate), size = MicrobioUoE::pts(10)) +
  scale_x_discrete(position = 'top', labels = scales::label_wrap(13)) +
  scale_y_discrete(position = 'left', labels = scales::label_wrap(13)) +
  labs(y = 'From',
  x = 'To',
  title = paste('all rates different with', length(mod_custom3$par), 'free parameters', sep = ' ')) +
  coord_fixed() +
  scale_color_manual(values = c('red', 'black'))

ggsave('plots/sequencing_rpoB/transition_matrix.png', last_plot(), height = 5, width = 7)

#------------------------------------#
# plot best transition rate model ####
#------------------------------------#

label_wrap2 <- function(x, width){
  unlist(lapply(strwrap(x, width = width, simplify = FALSE), 
                paste0, collapse = "\n"))
}

# calculate proportion of number of species each state has
d_habpref_summary <- group_by(d_meta, habitat_preference) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n/sum(n))

# set colours
cols_hab <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds'))

# make very simple plot to grab legend from
p_legend <- ggplot(d_habpref_summary, aes(habitat_preference, prop, col = habitat_preference)) +
  geom_point() +
  scale_color_manual('Biome preference', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=5)))
p_legend <- cowplot::get_legend(p_legend)

# turn transition matrix into network to plot
d_network <- as_tbl_graph(select(diversitree_df, state_1, state_2, transition_rate)) %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Biome preference', values = cols_hab) +
  scale_fill_manual('Biome preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_best <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white'))

p_best + (wrap_elements(p_legend)/plot_spacer()) + plot_layout(widths = c(0.8, 0.2))

# save out model
ggsave(here('plots/sequencing_rpoB/transition_plot_diversitree.png'), last_plot(), height = 5.5, width = 7)

# plot other models based on AIC weights

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_custom2, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Biome preference', values = cols_hab) +
  scale_fill_manual('Biome preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_2 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(b) Simplified ARD model 2',
       subtitle = 'AIC model weight = 0.24')

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_trans, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Biome preference', values = cols_hab) +
  scale_fill_manual('Biome preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_3 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(c) Stepwise model',
       subtitle = 'AIC model weight = 0.14')

# turn transition matrix into network to plot
d_network <- get_diversitree_df(mod_custom1, coding$hab_pref_num, coding$hab_pref) %>%
  select(., state_1, state_2, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(!is.na(transition_rate) & transition_rate > 0) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = round(transition_rate, 2)),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                strength = 1.3) + 
  geom_node_point(aes(size = prop,
                      col = name)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_size(range = c(5,20)) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_color_manual('Biome preference', values = cols_hab) +
  scale_fill_manual('Biome preference', values = cols_hab)

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_4 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', colour = 'white')) +
  labs(title = '(d) Simplified ARD model 1',
       subtitle = 'AIC model weight = 0.09')

# make layout
design <- "
  123
  453
"

p_best <- p_best + labs(title = '(a) Simplified ARD model 3',
                        subtitle = 'AIC model weight = 0.53')

p_best + p_2 + p_legend + p_3 + p_4 + plot_layout(design = design, widths = c(0.4, 0.4, 0.2))

# save out model
ggsave(here('plots/sequencing_rpoB/transition_plot_multiple.png'), last_plot(), height = 9, width = 12)
ggsave(here('plots/manuscript_plots/Figure_4.png'), last_plot(), height = 9, width = 12)

#------------------------------------------------#
# look at whether states are a source or sink ####
#------------------------------------------------#

d_source_sink_rate <- select(diversitree_df, away = state_1, into = state_2, transition_rate) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  group_by(habitat_preference, direction) %>%
  summarise(total_rate = sum(transition_rate), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_rate) %>%
  mutate(source_sink1 = away / into,
         source_sink2 = into - away,
         habitat_preference = gsub('terrestrial', 'land', habitat_preference))

table_rate <- select(d_source_sink_rate, habitat_preference, away, into, source_sink1) %>%
  mutate(across(away:source_sink1, ~round(.x, 2)),
         habitat_preference = gsub(':', ' + ', habitat_preference),
         habitat_preference = gsub('_', ' ', habitat_preference)) %>%
  arrange(desc(source_sink1)) %>%
  flextable(.) %>%
  set_header_labels(habitat_preference = 'biome preference',
                    source_sink1 = 'source sink ratio') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

# save out
save_as_image(table_rate, here('plots/manuscript_plots/source_sink_rate.png'), zoom = 3, webshot = 'webshot2')
save_as_docx(table_rate, path = here('plots/manuscript_plots/source_sink_rate.docx'))

# bootstrap ####

# bootstrap model by sampling 80% of the tips and refitting the best model

num_boots <- 1000

# set up empty tibble to store results
d_boots <- tibble(boot = 1:num_boots, 
                  rates = list(NA))

to_sample <- length(tree$tip.label) * 0.8
to_sample <- round(to_sample)

# run for loop to do bootstrapping
for(i in 1:num_boots){
  
  cat(i)
  
  # subset tree
  temp_tree <- keep.tip(tree, sample(tree$tip.label, to_sample, replace = FALSE))
  
  # subset habitat preference vector
  temp_hab_pref <- hab_pref_num[names(hab_pref_num) %in% temp_tree$tip.label]
  
  #sum(names(temp_hab_pref) == temp_tree$tip.label) == length(temp_tree$tip.label)
  
  # make model
  temp_ard <- make.mkn(temp_tree, temp_hab_pref, k = max(temp_hab_pref))
  
  # make custom matrix model
  temp_custom3 <- constrain(temp_ard, 
                            q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                            q45~0,
                            q54~0)
  
  # make start parameters
  inits_custom3 <- rep(1, length(argnames(temp_custom3)))
  
  # run model
  temp_mod_custom3 <- find.mle(temp_custom3, inits_custom3, method = 'nlminb', control = list(maxit = 50000))
  
  # save out parameters
  temp_output <- data.frame(temp_mod_custom3$par) %>%
    rownames_to_column(var = 'transition') %>%
    rename(rate = 2)
  
  d_boots$rates[[i]] <- temp_output
  
}

# create confidence intervals for each transition
d_boots2 <- unnest(d_boots, rates) %>%
  group_by(transition) %>%
  tidybayes::mean_qi(rate)

# make dataframe of model with the whole dataset
d_custom3 <- data.frame(mod_custom3$par) %>%
  rownames_to_column(var = 'transition') %>%
  rename(rate = 2)

# combine with diversitree df to get proper names of states
d_boots3 <- rename(d_boots2, param = transition) %>%
  select(-starts_with("boot")) %>%
  left_join(., select(diversitree_df, param, state_1, state_2, transition_rate)) %>%
  mutate(across(starts_with('state'), function(x) case_when(x == 'marine mud specialist' ~ 'marine specialist',
                                                            x == 'freshwater + terrestrial generalist' ~ 'freshwater + land generalist',
                                                            x == 'freshwater specialist' ~ 'freshwater specialist',
                                                            x == 'marine mud generalist' ~ 'marine generalist',
                                                            x == 'terrestrial specialist' ~ 'land specialist')),
         transition = paste(state_1, state_2, sep = ' -> '))

# make plot
ggplot(d_boots3) +
  geom_linerange(aes(x = forcats::fct_reorder(transition, rate), ymin = .lower, ymax = .upper)) +
  geom_point(aes(forcats::fct_reorder(transition, rate), transition_rate), shape = 21, size = 5, fill = 'red') +
  geom_point(aes(forcats::fct_reorder(transition, rate), rate), shape = 21, size = 3, fill = 'white') +
  theme_bw() +
  scale_x_discrete(labels = scales::label_wrap(25), guide = guide_axis(n.dodge = 2)) +
  labs(x = 'Transition',
       y = 'Transition rate')

# make network plot and table for source sink dynamics underneath
d_network <- rename(d_boots2, param = transition, transition_rate = rate) %>%
  select(-starts_with("boot")) %>%
  left_join(., select(diversitree_df, param, state_1, state_2)) %>%
  mutate(across(transition_rate:.upper, ~round(.x, 2)),
         rate = paste(transition_rate, ' (', .lower, '-', .upper, ')', sep = '')) %>%
  select(., state_1, state_2, rate, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  activate(nodes) %>%
  left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = rate),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                label_size = MicrobioUoE::pts(8),
                strength = 1.3) + 
  geom_node_point(aes(col = name, size = prop)) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_edge_width(range = c(0.5, 2), guide = 'none') +
  scale_color_manual('Biome preference', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  scale_size(range = c(5,20), guide = 'none') +
  guides(edge_alpha = 'none',
         color = guide_legend(override.aes = list(size = 3)))

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_2 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

# make table for source sink values
d_source_sink_v2 <- unnest(d_boots, rates) %>%
  left_join(., select(diversitree_df, transition = param, state_1, state_2)) %>%
  select(., boot, away = state_1, into = state_2, rate) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  mutate(habitat_preference = case_when(habitat_preference == 'marine mud specialist' ~ 'marine specialist',
                                        habitat_preference == 'freshwater + terrestrial generalist' ~ 'freshwater + land generalist',
                                        habitat_preference == 'freshwater specialist' ~ 'freshwater specialist',
                                        habitat_preference == 'marine mud generalist' ~ 'marine generalist',
                                        habitat_preference == 'terrestrial specialist' ~ 'land specialist')) %>%
  group_by(habitat_preference, direction, boot) %>%
  summarise(total_rate = sum(rate), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_rate) %>%
  mutate(source_sink1 = away / into) %>%
  group_by(habitat_preference) %>%
  ggdist::mean_qi()

table_2 <- d_source_sink_v2 %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  mutate(away = paste(away, ' (', away.lower, '-', away.upper, ')', sep = ''),
         into = paste(into, ' (', into.lower, '-', into.upper, ')', sep = ''),
         source_sink1 = paste(source_sink1, ' (', source_sink1.lower, '-', source_sink1.upper, ')', sep = '')) %>%
  arrange(desc(parse_number(source_sink1))) %>%
  select(habitat_preference, away, into, source_sink1) %>%
  flextable(.) %>%
  set_header_labels(habitat_preference = 'biome preference',
                    source_sink1 = 'source sink ratio') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

p_2 + gen_grob(table_2) + plot_layout(ncol = 1, heights = c(0.8, 0.2))

ggsave('plots/manuscript_plots/bootstrap_transitions.png', height = 5, width = 7)

#-------------------------------------------------------------#
# bootstrap tree to keep the same number of tips per group ####
#-------------------------------------------------------------#

# set number of bootstraps
num_boots <- 1000

# set up empty tibble to store results
d_boots_v2 <- tibble(boot = 1:num_boots, 
                     rates = list(NA))

# calculate number of tips per group
d_meta %>% 
  group_by(habitat_preference) %>%
  tally()

# there are 123 marine generalists so only keep 123 tips from each group
num_to_keep <- 123
  
# run for loop to do bootstrapping
for(i in 1:num_boots){
  
  # print iteration - bit shit but needs must
  cat(i)
  
  # choose tips to keep
  d_temp <- group_by(d_meta, habitat_preference) %>%
    sample_n(size = num_to_keep, replace = FALSE)
  
  # subset tree
  temp_tree <- keep.tip(tree, d_temp$tip_label)
  
  # reorder metadata to match tip labels of tree
  d_temp <- tibble(tip_label = temp_tree$tip.label) %>%
    left_join(., d_meta, by = 'tip_label')
  
  sum(d_temp$tip_label == temp_tree$tip.label) == length(temp_tree$tip.label)
  # SUCCESS if TRUE
  
  # create habitat preference vector
  temp_habpref <- setNames(d_temp$habitat_preference, d_temp$tip_label)
  temp_habpref <- as.numeric(as.factor(temp_habpref))
  temp_habpref <- setNames(temp_habpref, d_temp$tip_label)
  
  # make model
  temp_ard <- make.mkn(temp_tree, temp_habpref, k = max(temp_habpref))
  
  # make custom matrix model
  temp_custom3 <- constrain(temp_ard, 
                            q14~0, q24~0, q25~0, q35~0, q41~0, q42~0, q53~0,
                            q45~0,
                            q54~0)
  
  # make start parameters
  inits_custom3 <- rep(1, length(argnames(temp_custom3)))
  
  # run model
  temp_mod_custom3 <- find.mle(temp_custom3, inits_custom3, method = 'nlminb', control = list(maxit = 50000))
  
  # save out temp parameters
  temp_output <- data.frame(temp_mod_custom3$par) %>%
    rownames_to_column(var = 'transition') %>%
    rename(rate = 2)
  
  d_boots_v2$rates[[i]] <- temp_output
  
}

# create confidence intervals for each transition
d_boots2 <- unnest(d_boots_v2, rates) %>%
  group_by(transition) %>%
  tidybayes::mean_qi(rate)

# make dataframe of model with the whole dataset
d_custom3 <- data.frame(mod_custom3$par) %>%
  rownames_to_column(var = 'transition') %>%
  rename(rate = 2)

# combine with diversitree df to get proper names of states
d_boots3 <- rename(d_boots2, param = transition) %>%
  select(-starts_with("boot")) %>%
  left_join(., select(diversitree_df, param, state_1, state_2, transition_rate)) %>%
  mutate(across(starts_with('state'), function(x) case_when(x == 'marine mud specialist' ~ 'marine specialist',
                                                            x == 'freshwater + terrestrial generalist' ~ 'freshwater + land generalist',
                                                            x == 'freshwater specialist' ~ 'freshwater specialist',
                                                            x == 'marine mud generalist' ~ 'marine generalist',
                                                            x == 'terrestrial specialist' ~ 'land specialist')),
         transition = paste(state_1, state_2, sep = ' -> '))

# make plot
ggplot(d_boots3) +
  geom_linerange(aes(x = forcats::fct_reorder(transition, rate), ymin = .lower, ymax = .upper)) +
  geom_point(aes(forcats::fct_reorder(transition, rate), transition_rate), shape = 21, size = 5, fill = 'red') +
  geom_point(aes(forcats::fct_reorder(transition, rate), rate), shape = 21, size = 3, fill = 'white') +
  theme_bw() +
  scale_x_discrete(labels = scales::label_wrap(25), guide = guide_axis(n.dodge = 2)) +
  labs(x = 'Transition',
       y = 'Transition rate')

# make network plot and table for source sink dynamics underneath
d_network <- rename(d_boots2, param = transition, transition_rate = rate) %>%
  select(-starts_with("boot")) %>%
  left_join(., select(diversitree_df, param, state_1, state_2)) %>%
  mutate(across(transition_rate:.upper, ~round(.x, 2)),
         rate = paste(transition_rate, ' (', .lower, '-', .upper, ')', sep = '')) %>%
  select(., state_1, state_2, rate, transition_rate) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  activate(nodes) %>%
  #left_join(., select(d_habpref_summary, name = habitat_preference, prop, n)) %>%
  left_join(., tibble(name = names(cols_hab))) %>%
  mutate(order = c(2, 1, 5, 4, 3)) %>%
  arrange(order)

p <- ggraph(d_network, layout = 'linear', circular = TRUE) + 
  geom_edge_fan(aes(alpha = transition_rate, 
                    width = transition_rate,
                    label = rate),
                arrow = arrow(length = unit(4, 'mm')),
                end_cap = circle(10, 'mm'),
                start_cap = circle(10, 'mm'),
                angle_calc = 'along',
                label_dodge = unit(2.5, 'mm'),
                label_size = MicrobioUoE::pts(8),
                strength = 1.3) + 
  geom_node_point(aes(col = name), size = 15) +
  theme_void() +
  #geom_node_label(aes(label = label, x=xmin), repel = TRUE) +
  scale_edge_width(range = c(0.5, 2), guide = 'none') +
  scale_color_manual('Biome preference', values = cols_hab, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  guides(edge_alpha = 'none',
         color = guide_legend(override.aes = list(size = 3)))

# grab data for points
point_data <- p$data %>%
  select(x, y, name) %>%
  mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
         nudge_y = ifelse(y < 0, -0.2, 0.2))

p_2 <- p + 
  #geom_label(aes(nudge_x + x, nudge_y+y, label = label_wrap2(name, 15)), point_data, size = MicrobioUoE::pts(12)) +
  #coord_cartesian(clip = "off") +
  xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
  ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))


# make table for source sink values
d_source_sink_v2 <- unnest(d_boots_v2, rates) %>%
  left_join(., select(diversitree_df, transition = param, state_1, state_2)) %>%
  select(., boot, away = state_1, into = state_2, rate) %>%
  pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
  mutate(habitat_preference = case_when(habitat_preference == 'marine mud specialist' ~ 'marine specialist',
                                        habitat_preference == 'freshwater + terrestrial generalist' ~ 'freshwater + land generalist',
                                        habitat_preference == 'freshwater specialist' ~ 'freshwater specialist',
                                        habitat_preference == 'marine mud generalist' ~ 'marine generalist',
                                        habitat_preference == 'terrestrial specialist' ~ 'land specialist')) %>%
  group_by(habitat_preference, direction, boot) %>%
  summarise(total_rate = sum(rate), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = total_rate) %>%
  mutate(source_sink1 = away / into) %>%
  group_by(habitat_preference) %>%
  ggdist::mean_qi()

table_2 <- d_source_sink_v2 %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  mutate(away = paste(away, ' (', away.lower, '-', away.upper, ')', sep = ''),
         into = paste(into, ' (', into.lower, '-', into.upper, ')', sep = ''),
         source_sink1 = paste(source_sink1, ' (', source_sink1.lower, '-', source_sink1.upper, ')', sep = '')) %>%
  arrange(desc(parse_number(source_sink1))) %>%
  select(habitat_preference, away, into, source_sink1) %>%
  flextable(.) %>%
  set_header_labels(habitat_preference = 'biome preference',
                    source_sink1 = 'source sink ratio') %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', part = 'body', j = 1) %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()

p_2 + gen_grob(table_2) + plot_layout(ncol = 1, heights = c(0.8, 0.2))

ggsave('plots/manuscript_plots/bootstrap_transitions_v2.png', height = 5, width = 7)