# calculate GC content for each ASV

# load packages
librarian::shelf(Biostrings, seqinr, tidyverse, here, ggtree, padpadpadpad/MicrobioUoE, phytools, nlme, emmeans, flextable, multcomp)

# load in the tree
tree <- ape::read.tree(here('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre'))

# read in colours
cols_hab <- readRDS(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_colours.rds'))

# read in sequence to otu name converter
asv_names <- read_csv(here('data/sequencing_rpoB/phyloseq/myxococcus/myxo_seq_name_conversion.csv'))
asv_names

# read habitat_preference data
d_habpref <- read.csv(here('data/sequencing_rpoB/phyloseq/myxococcus/habitat_preference/summary/habitat_preference_asv_new.csv'))

# keep only sequences present in the filtered dataset
seqs <- filter(asv_names, otu_name %in% d_habpref$otu)

# calculate G/C content
d_gc <- nest(seqs, data = -otu_name) %>%
  mutate(gc = map_dbl(data, ~GC(s2c(.x$seq)))) %>%
  select(-data)

hist(d_gc$gc)

# left join this into habitat preference
d_habpref <- left_join(d_habpref, rename(d_gc, otu = otu_name))

# plot
ggplot(d_habpref, aes(habitat_preference3, gc*100)) +
  geom_pretty_boxplot(aes(col = habitat_preference3, fill = habitat_preference3), show.legend = FALSE) +
  geom_point(shape = 21, fill = 'white', col = 'black', position = position_jitter(width = 0.15), alpha = 0.8) +
  theme_bw(base_size = 16) +
  scale_x_discrete(labels = scales::label_wrap(15)) +
  scale_color_manual(values = cols_hab) +
  scale_fill_manual(values = cols_hab) +
  ylim(c(50, 75)) +
  labs(x = 'Habitat preference',
       y = 'GC content (%)')

ggsave('plots/sequencing_rpoB/analyses/gc_content/gc_content.png', last_plot(), width = 7, height = 5)

# do phylogenetic penalised regression to look at whether there are significant differences in GC across environments after accounting for phylogenetic distance

d_habpref <- mutate(d_habpref, tip_label = otu)

# set up correlation matrix for the tree
cor_lambda <- corPagel(value = 1, phy = tree, form = ~tip_label)

d_habpref <- mutate(d_habpref, hab_pref_fac = as.factor(habitat_preference3))

# fit phylogenetic generalised linear model
mod <- gls(gc ~ habitat_preference3, data = d_habpref, correlation = cor_lambda, method = 'ML')
mod_v2 <- gls(gc ~ hab_pref_fac, data = d_habpref, correlation = cor_lambda, method = 'ML')
mod1 <- gls(gc ~ 1, data = d_habpref, correlation = cor_lambda, method = 'ML')
summary(mod)
anova(mod)
anova(mod, mod1)

# do contrasts between habitat preferences
contrasts <- emmeans(mod, pairwise ~ habitat_preference3)
contrasts2 <- glht(mod_v2, linfct=mcp(hab_pref_fac='Tukey'))

contrasts
summary(contrasts2)

# produce table
table1 <- contrasts$contrasts %>%
  data.frame() %>%
  mutate(p.value = ifelse(p.value < 0.0001, "<0.0001", as.character(round(p.value, 2)))) %>%
  flextable() %>%
  set_header_labels(contrast = 'Contrast',
                    emmean = 'Estimate',
                    SE = 's.e.',
                    df = 'd.f.',
                    t.ratio = "t-ratio",
                    p.value = "p value") %>%
  italic(j = c(3:6), part = 'header') %>% # make some column names italic
  colformat_double(j = c(2:6), digits = 2) %>% # round numbers of specific columns to 2 decimal places
  align(align = 'center', part = 'header') %>% # align column names centrally
  align(align = 'left', part = 'body') %>% # align table output to the left
  font(fontname = 'Times', part = 'all') %>% # set font name for the table
  bold(~p.value == "<0.0001", j = "p.value") %>% # make significant p values bold
  fontsize(size = 12, part = 'all') %>% # set font size for the table
  autofit() # fix any random size issues

# save out table
save_as_image(table1, 'plots/sequencing_rpoB/analyses/gc_content/contrast_table.png', webshot = 'webshot2')
