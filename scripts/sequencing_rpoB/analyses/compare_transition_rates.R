# compare the best models for the diversitree models

# load in packages
librarian::shelf(ape, diversitree, flextable, padpadpadpad/MicrobioUoE, tidyverse, patchwork)

# list files
best_asv <- 'data/sequencing_rpoB/processed/transition_rates/mod_custom_3.rds'
best_95 <- 'data/sequencing_rpoB/processed/transition_rates/mod_custom_6_otu95.rds'
best_97.7 <- 'data/sequencing_rpoB/processed/transition_rates/mod_custom_3_otu97.7.rds'
best_bootstraps <- list.files("data/sequencing_rpoB/processed/transition_rates", full.names = TRUE, pattern = 'asvboot')

# write a function to extract data we want from the diversitree models
get_transition_rate_df <- function(mod){
  temp <- readRDS(mod)
  temp <- data.frame(temp$par.full) %>%
    dplyr::rename(rate = 1) %>%
    rownames_to_column('transition')
  return(temp)
  }

# combine all the files together and run the function
d_rates <- tibble(file = c(best_asv, best_95, best_97.7, best_bootstraps),
                model = c('best_asv', 'best_95', 'best_97.7', paste('asv_bootstrap', 1:9, sep ='_'))) %>%
  mutate(data = map(file, get_transition_rate_df)) %>%
  unnest(data) %>%
  mutate(from = as.numeric(str_sub(transition, 2,2)),
         to = as.numeric(str_sub(transition, 3,3)))

# add in the coding
coding <- read.csv('data/sequencing_rpoB/processed/transition_rates/coding.csv')

d_rates <- left_join(d_rates, select(coding, to_name = hab_pref, to = hab_pref_num)) %>%
  left_join(., select(coding, from_name = hab_pref, from = hab_pref_num)) %>%
  mutate(transition2 = paste(from_name, to_name, sep = ' -> ')) %>%
  mutate(rate2 = ifelse(rate == 0, NA, rate))

# make prevalence plot for transitions
d_prevalence <- group_by(d_rates, transition2, from_name, to_name) %>%
  summarise(present = sum(!is.na(rate2)),
            total = n(),
            prop = present/total, .groups = 'drop') %>%
  arrange(prop) %>%
  mutate(order = 1:n(),
         transition2 = gsub('-> ', '->\n', transition2))

p1 <- ggplot(d_prevalence, aes(prop, forcats::fct_reorder(transition2, order))) +
  geom_col(fill = 'light grey', col = 'black') +
  labs(x = 'Proportion of models\nwith transition present', y = 'Transition') +
  theme_bw(base_size = 12)

# make plot comparing rates of transitions between models
p2 <- filter(d_rates) %>%
  mutate(transition2 = gsub('-> ', '->\n', transition2)) %>%
  select(rate2, transition2) %>%
  left_join(., select(d_prevalence, transition2, order)) %>%
  ggplot(., aes(rate2, forcats::fct_reorder(transition2, order))) +
  geom_pretty_boxplot(fill = 'light grey', col = 'light grey') +
  geom_point(position = position_jitter(height = 0.1), shape = 21, fill = 'white') +
  labs(x = 'Transition rate', y = 'Transition') +
  theme_bw(base_size = 12)

p2

p1 + p2 + theme(axis.text.y = element_blank(),
                axis.title.y = element_blank())

ggsave('plots/manuscript_plots/transition_rates_compare_plot.png', width = 8, height = 8)

# make table
d_table <- select(d_rates, model, rate2, transition2, from_name, to_name) %>%
  mutate(rate2 = ifelse(is.na(rate2), '-', round(rate2, 2))) %>%
  pivot_wider(names_from = model, values_from = rate2) %>%
  select(from_name, to_name, best_asv, contains('asv'), everything(), -transition2)

table <- flextable(d_table) %>%
  set_header_labels(from_name = 'from',
                    to_name = 'to',
                    best_asv = 'ASV tree',
                    asv_bootstrap_1 = 'ASV\nbootstrap 1',
                    asv_bootstrap_2 = 'ASV\nbootstrap 2',
                    asv_bootstrap_3 = 'ASV\nbootstrap 3',
                    asv_bootstrap_4 = 'ASV\nbootstrap 4',
                    asv_bootstrap_5 = 'ASV\nbootstrap 5',
                    asv_bootstrap_6 = 'ASV\nbootstrap 6',
                    asv_bootstrap_7 = 'ASV\nbootstrap 7',
                    asv_bootstrap_8 = 'ASV\nbootstrap 8',
                    asv_bootstrap_9 = 'ASV\nbootstrap 9',
                    best_95 = 'OTU 95%\ntree',
                    best_97.7 = 'OTU 97.7%\ntree') %>%
  align(align = 'center', part = 'all') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  hline(i = c(4,8, 12, 16)) %>%
  vline(j = c(2,3,12,13)) %>%
  autofit() %>%
  width(j = 1:2, width = 1.5) 

table  

# save out table
save_as_image(table, 'plots/manuscript_plots/transition_rates_compare.png')
