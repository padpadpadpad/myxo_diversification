# script to create a table comparing estimated transition rates from the two best models (Ubuntu and MacOS)

# load packages
library(tidyverse)
library(flextable)
library(diversitree)

# calculate all text orders of a vector of characters
get_all_combs <- function(x, vec_len = 2){
  # set command
  return(expand.grid(rep(list(x), vec_len)) %>%
           filter(apply(., 1, function(x)length(unique(x))) == vec_len) %>% 
           unite(., 'y', everything(), sep = '') %>% 
           pull(y))
}

# load in Ubuntu best model
mod_ubuntu <- readRDS('data/sequencing_rpoB/processed/diversitree_mk_model.rds')

AIC(mod_ubuntu)

# load best macOS model
mod_mac <- readRDS('data/sequencing_rpoB/processed/diversitree_mk_model_mac.rds')

AIC(mod_mac)

# load in habitat preference
hab_pref <- readRDS('data/sequencing_rpoB/processed/hab_pref_vec.rds')

# make habitat preference vector numeric
hab_pref_num <- as.numeric(as.factor(hab_pref)) -1
hab_pref_num <- setNames(hab_pref_num, names(hab_pref))
hab_pref_num2 <- hab_pref_num + 1

# coding from numeric to character
coding <- tibble(hab_pref = unname(hab_pref), hab_pref_num = unname(hab_pref_num)) %>%
  distinct() %>%
  mutate(hab_pref2 = gsub(':', '.', hab_pref),
         hab_pref_num2 = hab_pref_num + 1) %>%
  arrange(hab_pref) %>%
  mutate(hab_pref_axis = gsub(':', '/ ', hab_pref),
         hab_pref_axis = gsub('_', ' ', hab_pref_axis),
         # rename the columns for easy naming
         initials = c('F', 'FM', 'FT', 'G', 'M', 'MT', 'T'))

coding

# create q matrix data frame

# for ubuntu

# grab names of q matrices from diversitree
q_matrix_ubuntu <- names(mod_ubuntu$par)

# create dataframe with these
q_matrix_ubuntu <- tibble(q_matrix_ubuntu = q_matrix_ubuntu) %>%
  # grab first and second elements that are the habitat codings
  mutate(state_1_num = as.numeric(substr(q_matrix_ubuntu, 2,2)),
         state_2_num = as.numeric(substr(q_matrix_ubuntu, 3,3)),
         transition_rate_diversitree = unlist(mod_ubuntu$par)) %>%
  left_join(select(coding, state_1 = hab_pref, state_1_num = hab_pref_num2)) %>%
  left_join(select(coding, state_2 = hab_pref, state_2_num = hab_pref_num2)) %>%
  select(-state_1_num, -state_2_num, -q_matrix_ubuntu)

q_matrix_mac <- tibble(q_matrix_mac = names(mod_mac$par)) %>%
  # grab first and second elements that are the habitat codings
  mutate(state_1_num = as.numeric(substr(q_matrix_mac, 2,2)),
         state_2_num = as.numeric(substr(q_matrix_mac, 3,3)),
         transition_rate_diversitree = unlist(mod_mac$par)) %>%
  left_join(select(coding, state_1 = hab_pref, state_1_num = hab_pref_num2)) %>%
  left_join(select(coding, state_2 = hab_pref, state_2_num = hab_pref_num2)) %>%
  select(-state_1_num, -state_2_num, -q_matrix_mac)

nrow(q_matrix_mac)
nrow(q_matrix_ubuntu)

# create all combinations

all_combs <- expand.grid(rep(list(coding$hab_pref), 2)) %>%
  filter(Var1 != Var2) %>%
  rename(state_1=Var1, state_2 = Var2) %>%
  left_join(., rename(q_matrix_mac, transition_rate_mac = transition_rate_diversitree)) %>%
  left_join(., rename(q_matrix_ubuntu, transition_rate_ubuntu = transition_rate_diversitree)) %>%
  mutate(across(starts_with('transition'), replace_na, 0))

# how do they compare?
ggplot(all_combs, aes(transition_rate_ubuntu, transition_rate_mac)) +
  geom_point() +
  theme_bw()
# sort of good in the main

# lets visualise the differences in a table

table <- mutate(all_combs, ave = (transition_rate_ubuntu + transition_rate_mac)/2) %>%
  arrange(-ave) %>%
  select(-ave) %>%
  mutate(across(starts_with('state'), function(x) gsub('mud_and_shore', 'marine mud', x)),
         across(starts_with('state'), function(x) gsub(':', ' + ', x)),
         across(starts_with('trans'), round, 2))

table_flex <- flextable(table) %>%
  set_header_labels(state_1 = 'from',
                    state_2 = 'to',
                    transition_rate_ubuntu = 'ubuntu',
                    transition_rate_mac = 'mac') %>%
  align(align = 'center', part = 'header') %>%
  bold(part = 'header') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit()
  
# save out
save_as_image(table_flex, 'plots/sequencing_rpoB/analyses/transition_rate_comparison.png', zoom = 3, webshot = 'webshot2')
