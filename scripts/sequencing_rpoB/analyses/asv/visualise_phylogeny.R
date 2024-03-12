# plot phylogenetic tree and visualise biome preference across the tree

#--------------------------#
# what this script does ####
#--------------------------#

# reads in the best ultrametric ASV Myxobacteria tree
# plots the tree
# looks at the distribution of biome preference


#------------------------------#
# load in packages and data ####
#------------------------------#

# load packages - use librarian::shelf()
# need to ensure BiocManager and Biobase are installed - BiocManager::install("Biobase")
# can install packages from GitHub, Bioconductor, and CRAN all at once
librarian::shelf(here, tidyverse, ggtree, ggnewscale, RColorBrewer, patchwork, ape, phytools, MetBrewer, flextable, officer, magick)

# load in the tree
tree <- ape::read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')
#tree <- phytools::force.ultrametric(tree)

#tree <- ape::write.tree(tree, 'data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

# replace NA in node label with ""
tree$node.label[is.na(tree$node.label)] <- ''

is.ultrametric(tree)

# alter tip labels to remove family as they will not link to the distance matrix
# write function to remove family labels
strsplit_mod <- function(x)(strsplit(x, split = '_') %>% unlist() %>% .[1:2] %>% paste0(., collapse = '_'))

tree$tip.label <- purrr::map_chr(tree$tip.label, strsplit_mod)

# read in habitat preference
d_habpref <- read.csv('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv.csv')

# read in phyloseq object and grab tax table
d_taxa <- readRDS('data/sequencing_rpoB/phyloseq/myxococcus/prevalence_filtered/ps_otu_asv_filt.rds') %>%
  phyloseq::tax_table() %>%
  data.frame() %>%
  janitor::clean_names() %>%
  rownames_to_column('otu')


#----------------------#
# do some wrangling ####
#----------------------#

# create d_meta
d_meta <- left_join(select(d_habpref, otu, habitat_preference, num_present), select(d_taxa, otu:family)) %>%
  mutate(family2 = ifelse(is.na(family), 'uncertain', family))

# set habitat colours
cols_hab <- met.brewer('Austria', n = 7) %>% as.character()
cols_labels <- c('marine_mud', 'freshwater', 'terrestrial', 'freshwater:terrestrial', 'freshwater:marine_mud:terrestrial', 'marine_mud:terrestrial', 'freshwater:marine_mud')

# create a named vector of the colours to match up with d_meta
names(cols_hab) <- cols_labels

# create new habitat labels
hab_labels <- c('marine mud specialist', 'freshwater specialist', 'terrestrial specialist', 'freshwater + terrestrial generalist', 'full generalist', 'marine mud + terrestrial generalist', 'freshwater + marine mud generalist')
names(hab_labels) <- cols_labels

# sort the new labels to be in the correct order
hab_labels <- hab_labels[sort(names(hab_labels))]

#------------------#
# plot the tree ####
#------------------#

# constrained families
constrained_families <- c('Myxococcaceae', 'Vulgatibacteraceae', 'Anaeromyxobacteraceae', 'Polyangiaceae', 'Sandaracinaceae', 'Nannocystaceae', 'Haliangiaceae')

# find proportion of tips assigned to these families
d_common <- d_meta %>%
  group_by(family2) %>%
  tally() %>%
  ungroup() %>%
  mutate(., prop = n / sum(n)) %>%
  filter(., family2 %in% constrained_families)

sum(d_common$prop)

# add in column for colouring branches according to family - add in uncertain (family = NA), and unconstrained values (family != NA but not one of the 7 families used to constrain the tree)
d_meta <- mutate(d_meta, family_common = ifelse(family2 %in% c(constrained_families, 'uncertain'), family2, 'unconstrained')) %>%
  rename(tip_label = otu)

# group tip labels together in terms of their order
to_group <- split(d_meta$tip_label, d_meta$family_common)
tree2 <- groupOTU(tree, to_group)

# find the mrca of each of the clades
# this is used to colour the outerbars of the phylogeny
d_meta2 <- filter(d_meta, family_common %in% constrained_families) %>%
  select(family, tip_label) %>%
  group_by(family) %>%
  nest() %>%
  ungroup()

# run for loop to find MRCA of each set of tip labels
for(i in 1:nrow(d_meta2)){
  d_meta2$mrca[i] <- findMRCA(tree2, tips = d_meta2$data[[i]]$tip_label)
}

# add blank label so family names are not printed onto the graph
d_meta2 <- select(d_meta2, family2 = family, mrca) %>%
  mutate(blank_label = '')

# create colour palettes

# for different families - add in black and grey for uncertain and uncommon respectively
cols <- c(colorRampPalette(brewer.pal(11, "Spectral"))(nrow(d_common)), 'black', 'grey')
names(cols) <- c(sort(d_common$family2), 'uncertain', 'unconstrained')

# make a separate column for the rare states to make their size bigger!
d_meta <- mutate(d_meta, rare = ifelse(habitat_preference %in% c('freshwater:marine_mud:terrestrial', 'marine_mud:terrestrial', "freshwater:marine_mud"), 'rare', 'common'))

castor::get_all_distances_to_root(tree2) %>% max()

# first only plot taxonomy, also make non circular so we can see groupings properly
tree_plot <- ggtree(tree2, aes(col = group)) %<+% filter(d_meta) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +
  geom_cladelab(data = d_meta2,
                mapping = aes(node = mrca,
                              color = family2,
                              label = blank_label),
                offset = castor::get_all_distances_to_root(tree2) %>% max() * 0.03,
                barsize = 2,
                show.legend = FALSE)

tree_plot

# we can remake the tree in a circular fashion for all the other plots from now on
tree_plot <- ggtree(tree2, layout = 'circular', aes(col = group)) %<+% filter(d_meta) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +
  geom_cladelab(data = d_meta2,
                mapping = aes(node = mrca,
                              color = family2,
                              label = blank_label),
                offset = castor::get_all_distances_to_root(tree2) %>% max() * 0.08,
                barsize = 2,
                show.legend = FALSE)

tree_plot

#----------------------------------------------#
# plot the  tree with biome preference info ####
#----------------------------------------------#

# look at numbers in each biome preference
group_by(d_meta, habitat_preference) %>%
  summarise(n = n(),
            ave_present = median(num_present)) %>%
  arrange(desc(n)) %>%
  mutate(prop = round(n/sum(n)*100, 2))

# look at number that only turned up in single - or multiple - biomes
d_habpref %>%
  group_by(habitats_present) %>%
  tally() %>%
  mutate(prop = round(n/sum(n), 2))

d_habpref %>%
  group_by(habitats_present, habitat_preference) %>%
  tally() %>%
  group_by(habitat_preference) %>%
  mutate(prop = round(n/sum(n), 2)) %>%
  arrange(habitat_preference, habitats_present)

# remake tree and add in biome preference
# reuse tree made above
# now plot habitat preference around the tips
tree_plot2 <- tree_plot +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  new_scale_color() +
  geom_tippoint(aes(x=x+x*0.04, col = habitat_preference, size = rare), position = position_jitter(width = 0.025, height = 0)) +
  scale_size_manual(values = c(0.6, 3)) +
  scale_color_manual('Habitat preference (tip points)', values = cols_hab, labels = hab_labels) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = 'none')

tree_plot2 

# save tree out
ggsave(here('plots/sequencing_rpoB/tree_all_taxonomy.png'), tree_plot2, height = 9, width = 12)

#-------------------------------------------#
# plot tree after collapsing rare states ####
#-------------------------------------------#

# collapse full generalists, marine + terrestrial generalist, and freshwater + marine generalist into a single preference
# create new habitat preference vector
d_habpref <- mutate(d_habpref, habitat_preference2 = ifelse(habitat_preference %in% c('freshwater:marine_mud:terrestrial', 'freshwater:marine_mud', 'marine_mud:terrestrial'), 'marine_mud_generalist', habitat_preference),
                    # rename all habitat preference vectors for easy renaming
                    habitat_preference3 = case_when(habitat_preference2 == 'marine_mud_generalist' ~ "marine mud generalist",
                                                    habitat_preference2 == 'marine_mud' ~ "marine mud specialist",
                                                    habitat_preference == "terrestrial" ~ "terrestrial specialist",
                                                    habitat_preference2 == "freshwater" ~ "freshwater specialist",
                                                    habitat_preference2 == "freshwater:terrestrial" ~ "freshwater + terrestrial generalist"))

d_habpref %>% group_by(habitat_preference3) %>%
  tally() %>%
  arrange(-n)

# make table used in Figure 2
d_table <- group_by(d_habpref, habitat_preference2, habitat_preference) %>%
  tally() %>%
  mutate(habitat_preference = case_when(habitat_preference == 'freshwater' ~ 'freshwater specialist',
                                        habitat_preference == 'freshwater:marine_mud' ~ 'freshwater + marine generalist',
                                        habitat_preference == 'freshwater:marine_mud:terrestrial' ~ 'full generalist',
                                        habitat_preference == 'freshwater:terrestrial' ~ 'freshwater + land generalist',
                                        habitat_preference == 'marine_mud:terrestrial'~ 'land + marine generalist',
                                        habitat_preference == 'marine_mud'~ 'marine specialist',
                                        habitat_preference == 'terrestrial'~ 'land specialist'),
         habitat_preference2 = case_when(habitat_preference2 == 'freshwater' ~ 'freshwater specialist',
                                         habitat_preference2 == 'marine_mud'~ 'marine specialist',
                                         habitat_preference2 == 'terrestrial'~ 'land specialist',
                                         habitat_preference2 == 'freshwater:terrestrial' ~ 'freshwater + land generalist',
                                         habitat_preference2 == 'marine_mud_generalist' ~ 'marine generalist')) %>%
  ungroup()

head(d_table)

# rearrange table
d_table <- mutate(d_table, first_col = habitat_preference2) %>%
  select(first_col, everything())

black_dot <- "\U2B24" 

# make table in flextable
table <- flextable(d_table) %>%
  set_header_labels(first_col = '',
                    habitat_preference2 = 'Biome preference',
                    habitat_preference = 'Original preference',
                    n = 'Number of ASVs') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  compose(j = 'first_col', i=1,
          value = as_paragraph(as_chunk(black_dot, props = fp_text(color = '#16317d'))),
          part = "body") %>%
  compose(j = 'first_col', i=2,
          value = as_paragraph(as_chunk(black_dot, props = fp_text(color = '#ffcd12'))),
          part = "body") %>%
  compose(j = 'first_col', i=3,
          value = as_paragraph(as_chunk(black_dot, props = fp_text(color = '#a40000'))),
          part = "body") %>%
  compose(j = 'first_col', i=4,
          value = as_paragraph(as_chunk(black_dot, props = fp_text(color = '#721b3e'))),
          part = "body") %>%
  compose(j = 1:2, i=5,
          value = as_paragraph(as_chunk('')),
          part = "body") %>%
  compose(j = 1:2, i=6,
          value = as_paragraph(as_chunk('')),
          part = "body") %>%
  compose(j = 'first_col', i=7,
          value = as_paragraph(as_chunk(black_dot, props = fp_text(color = '#007e2f'))),
          part = "body") %>%
  hline(i = c(1,2,3,6), border = fp_border_default()) %>%
  autofit()

save_as_image(table, 'plots/manuscript_plots/hab_pref_table.png')

# remake tree with collapsed biome preference
# remake colour scheme
cols_hab2 <- cols_hab[names(cols_hab) %in% c('marine_mud', 'freshwater', 'terrestrial', 'freshwater:terrestrial')]
cols_hab2 <- c(cols_hab2, '#721b3e')
names(cols_hab2) <- c('marine mud specialist', 'freshwater specialist', 'terrestrial specialist', 'freshwater + terrestrial generalist', 'marine mud generalist')

# create d_meta
d_meta_new <- left_join(select(d_habpref, otu, habitat_preference3, num_present), select(d_taxa, otu:family)) %>%
  mutate(rare2 = ifelse(habitat_preference3 == "marine mud generalist", 'rare', "common"))

# make plot again
tree_plot3 <- tree_plot %<+% d_meta_new +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  new_scale_color() +
  geom_tippoint(aes(x=x+x*0.04, col = habitat_preference3, size = rare2), position = position_jitter(width = 0.025, height = 0)) +
  scale_color_manual('Biome preference (tip points)', values = cols_hab2, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
  scale_size_manual(values = c(0.6, 3)) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = 'none')

tree_plot3

#----------------------------------#
# make extra plots for Figure 2 ####
#----------------------------------#

# create lineage through time plot
d_ltt <-  ape::ltt.plot.coords(tree) %>%
  data.frame() %>%
  mutate(time2 = time + 1)

# create lineage through time plot
p2 <- ggplot(d_ltt, aes(time, log(N))) +
  geom_line() +
  theme_bw(base_size = 12) +
  labs(x = 'Relative time from present',
       y = 'Log number of lineages')

# make constraint tree
backbone <- read.tree(here('data/sequencing_rpoB/raxml/backbone.tre'))

backbone_tree <- ggtree(backbone) +
  geom_tippoint() +
  geom_tiplab(offset = 0.1, size = MicrobioUoE::pts(10)) +
  coord_cartesian(clip = "off") +
  theme_tree(plot.margin=margin(6, 120, 6, 6))

# make table with flex table
legend <- cowplot::get_legend(tree_plot3)

ggsave('plots/sequencing_rpoB/tree_with_constraints.png', tree_plot3 + theme(legend.position = 'none'), height = 8, width = 8)

#----------------------#
# assemble Figure 2 ####
#----------------------#

tree_plot3 <- image_read('plots/sequencing_rpoB/tree_with_constraints.png', density = 300)
tree_plot3 <- image_trim(tree_plot3)

table <- png::readPNG('plots/manuscript_plots/hab_pref_table.png', native = TRUE)

layout <- c(
  'AAAAAAAB
   AAAAAAAC
   DDDDEEEE'
)

image_ggplot(tree_plot3) + 
  legend + 
  backbone_tree + 
  p2 + 
  table + 
  plot_layout(design = layout, heights = c(0.6, 0.25, 0.15)) +
  plot_annotation(tag_levels = list(c('(a)', '','', '(b)', '(c)'))) &
  theme(plot.tag = element_text(size = 14))

# save out
ggsave('plots/manuscript_plots/Figure_2.png', last_plot(), height = 12, width = 13)

#----------------------------------------------------#
# save out objects used in other analysis scripts ####
#----------------------------------------------------#

# save out habitat preference vector
write.csv(d_habpref, 'data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv', row.names = FALSE)
saveRDS(cols_hab2, 'data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds')

#----------------------------------------------------#
# create plot looking at phylogenetic uncertainty ####
#----------------------------------------------------#

# the big whole plot is a complete mess when looking at phylogenetic uncertainty - so try plot family specific trees

# replace NA node label with ''
tree$node.label <- ifelse(tree$node.label == 'NA', '', tree$node.label)
tree$node.label <- substr(tree$node.label, 1, 4)
tree2$node.label <- ifelse(tree$node.label == 'NA', '', tree$node.label)
tree2$node.label <- substr(tree$node.label, 1, 4)

tree$node.label[1:6]

# first plot is to plot bootstrap values against distance from root
# calculate distance from root to each internal node of the tree
dist_to_root <- castor::get_all_distances_to_root(tree)
dist_to_root <- dist_to_root[(Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))]

d_dist <- data.frame(dist_to_root = dist_to_root,
                     bootstrap = as.numeric(tree$node.label))

p1 <- ggplot(d_dist, aes(dist_to_root, bootstrap)) +
  geom_point(alpha = 0.3) +
  theme_bw(base_size = 12) +
  labs(y = 'Bootstrap tbe value',
       x = 'Distance from root',
       title = '(a)')

p2 <- ggplot(d_dist, aes(bootstrap)) +
  geom_histogram(col = 'black', fill = 'light grey') +
  theme_bw(base_size = 12) +
  labs(y = 'Number of nodes',
       x = 'Bootstrap tbe value',
       title = '(b)')

p1 + p2

# save this plot out
ggsave('plots/manuscript_plots/phylogeny_checks/bootstrap_vs_distance.png', last_plot(), height = 3, width = 8)

tree_plot <- ggtree(tree2, layout = 'circular', aes(col = group))  %<+% filter(d_meta) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot + geom_label(aes(label = node), col = 'black', size = MicrobioUoE::pts(8))

# HORRIBLE

# grab legend though
legend <- cowplot::get_legend(tree_plot)

# so plot each family separately, then plot the nodes not in these for all the whole tree

# identify MRCA of all Anaeromyxobacteraceae
node_anaero <- getMRCA(tree, filter(d_meta, family == 'Anaeromyxobacteraceae')$tip_label)

tree_anaeromyxobacteraceae <- extract.clade(tree, node_anaero)

tree_anaeromyxobacteraceae_2 <- groupOTU(tree_anaeromyxobacteraceae, to_group)

tree_plot_anaero <- ggtree(tree_anaeromyxobacteraceae_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_anaeromyxobacteraceae$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_anaero + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/anaeromyxobacteraceae.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Haliangiaceae
node_haliangiaceae <- getMRCA(tree, filter(d_meta, family == 'Haliangiaceae')$tip_label)

tree_haliangiaceae <- extract.clade(tree, node_haliangiaceae)

tree_haliangiaceae_2 <- groupOTU(tree_haliangiaceae, to_group)

tree_plot_haliangiaceae <- ggtree(tree_haliangiaceae_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_haliangiaceae$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_haliangiaceae + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/haliangiaceae.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Myxococcaceae
node_myxo <- getMRCA(tree, filter(d_meta, family == 'Myxococcaceae')$tip_label)

tree_myxo <- extract.clade(tree, node_myxo)

tree_myxo_2 <- groupOTU(tree_myxo, to_group)

tree_plot_myxo <- ggtree(tree_myxo_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_myxo$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_myxo + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/myxo.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Nannocystaceae
node_nanno <- getMRCA(tree, filter(d_meta, family == 'Nannocystaceae')$tip_label)

tree_nanno <- extract.clade(tree, node_nanno)

tree_nanno_2 <- groupOTU(tree_nanno, to_group)

tree_plot_nanno <- ggtree(tree_nanno_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_nanno$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_nanno + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/nanno.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Polyangiacaeae
node_poly <- getMRCA(tree, filter(d_meta, family == 'Polyangiaceae')$tip_label)

tree_poly <- extract.clade(tree, node_poly)

tree_poly_2 <- groupOTU(tree_poly, to_group)

tree_plot_poly <- ggtree(tree_poly_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_poly$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_poly + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/poly.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Sandaracinaceae
node_sanda <- getMRCA(tree, filter(d_meta, family == 'Sandaracinaceae')$tip_label)

tree_sanda <- extract.clade(tree, node_sanda)

tree_sanda_2 <- groupOTU(tree_sanda, to_group)

tree_plot_sanda <- ggtree(tree_sanda_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_sanda$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_sanda + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/sanda.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# identify MRCA of all Vulgatibacteraceae
node_vulga <- getMRCA(tree, filter(d_meta, family == 'Vulgatibacteraceae')$tip_label)

tree_vulga <- extract.clade(tree, node_vulga)

tree_vulga_2 <- groupOTU(tree_vulga, to_group)

tree_plot_vulga <- ggtree(tree_vulga_2, aes(col = group))  %<+% filter(d_meta, tip_label %in% tree_vulga$tip.label) +
  scale_color_manual('Family (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_plot_vulga + ggnewscale::new_scale_color() +
  geom_nodelab(aes(label = label), geom = 'text', col = 'black', size = MicrobioUoE::pts(6))

ggsave('plots/manuscript_plots/phylogeny_checks/vulga.jpeg',
       last_plot(),
       height = 10,
       width = 8)

# plot tree with collapsed nodes

# do not show legend
tree_plot <- ggtree(tree2, aes(col = group), layout = 'circular')  %<+% filter(d_meta) +
  scale_color_manual('Family (branch colours)', values = cols) +
  # turn off legend
  theme(legend.position = 'none') +
  geom_text(aes(label = label), col = 'black', size = MicrobioUoE::pts(12))

tree_plot <- tree_plot %>%
  collapse(node_anaero, 'max', fill = cols[1], alpha = 0.5) %>% 
  collapse(node_haliangiaceae, 'max', fill = cols[2], alpha = 0.5) %>%
  collapse(node_myxo, 'max', fill = cols[3], alpha = 0.5) %>%
  collapse(node_nanno, 'max', fill = cols[4], alpha = 0.5) %>%
  collapse(node_poly, 'max', fill = cols[5], alpha = 0.5) %>%
  collapse(node_sanda, 'max', fill = cols[6], alpha = 0.5) %>%
  collapse(node_vulga, 'max', fill = cols[7], alpha = 0.5) 

# save out tree plot
ggsave('plots/manuscript_plots/phylogeny_checks/whole_tree.jpeg',
       tree_plot,
       height = 10,
       width = 10)

tree_plot <- image_read('plots/manuscript_plots/phylogeny_checks/whole_tree.jpeg', density = 300)
tree_plot <- image_trim(tree_plot)

image_ggplot(tree_plot) + 
  legend +
  plot_layout(widths = c(0.8, 0.2))

# save out
ggsave('plots/manuscript_plots/phylogeny_checks/whole_tree.jpeg', last_plot(), height = 10, width = 11)

