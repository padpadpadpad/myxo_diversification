# lets look at which hidden model is best

# load packages
librarian::shelf(secsse, tidyverse)

# list all model files
files <- list.files('data/sequencing_rpoB/processed/secsse', full.names = TRUE)

# functions to clean model objects

# custom label_wrap function
label_wrap_mod <- function (x, width){
  unlist(lapply(strwrap(x, width = width, simplify = FALSE), 
                paste0, collapse = "\n"))
}

# calculate aic and get model summaries
calc_aic <- function(file, nhidden = 0){
  temp <- readRDS(file)
  
  d_temp <- data.frame(nhidden = nhidden,
                       nparams = temp$n_params,
                       loglik = temp$mod$ML,
                       conv = temp$mod$conv,
                       samp_frac = unique(temp$samp_frac),
                       run = str_extract(file, "(?<=run)\\d+(?=\\.rds)"))
  # calculate aic and specify which set of initial values are used
  d_temp <- mutate(d_temp, aic = (2*nparams) - 2*loglik,
                   initial_values = case_when(run %in% c(1, 4, 7, 10) ~ 'mk transitions',
                                              run %in% c(2, 5, 8, 11) ~ 'double speciation halve transitions',
                                              run %in% c(3, 6, 9, 12) ~ 'halve speciation double transitions'))

  return(d_temp)
}

# grab parameter values out
grab_params <- function(file){
  temp <- readRDS(file)

  d_temp_speciation <- data.frame(estimate = temp$mod$MLpars[[1]]) %>%
    rownames_to_column(var = 'param') %>%
    mutate(type = 'speciation',
           param = paste('lambda', param, sep = ''))
  d_temp_extinction <- data.frame(estimate = temp$mod$MLpars[[2]]) %>%
    rownames_to_column(var = 'param') %>%
    mutate(type = 'extinction',
           param = paste('mu', param, sep = ''))
  d_temp_transition <- data.frame(temp$mod$MLpars[[3]]) %>%
    rownames_to_column(var = 'from') %>%
    pivot_longer(names_prefix = 'X', cols = starts_with('X'), values_to = 'estimate', names_to = 'to') %>%
    mutate(param = paste('q', from, to, sep = ''),
           type = 'transition') %>%
    select(-from, -to)
  
  return(bind_rows(d_temp_speciation, d_temp_extinction, d_temp_transition) %>% mutate(run = str_extract(file, "(?<=run)\\d+(?=\\.rds)")))
  }

#-----------------------------#
# look at the musse models ####
#-----------------------------#

files_musse <- files[str_detect(files, 'musse')]

# calculate initial values
d_musse_aic <- map_df(files_musse, calc_aic) %>%
  arrange(samp_frac, aic) %>%
  mutate(model = 'musse')

# plot these
ggplot(d_musse_aic, aes(as.character(samp_frac), aic)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(aes(label = label_wrap_mod(initial_values, 20)), force_pull = 0.5, box.padding = 0.5) +
  theme_bw()

# mk transitions are the best all the time
# i.e. initial values do matter (QUITE A LOT)

# look at variation in rates for each initial rate values
d_musse_rates <- map_df(files_musse, grab_params) %>%
  left_join(., d_musse_aic) %>%
  group_by(param) %>%
  mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(estimate))

ggplot(d_musse_rates, aes(fct_reorder(param, estimate, .fun = mean), estimate)) +
  geom_point(aes(col = initial_values)) +
  facet_wrap(~type+ samp_frac, scale = 'free', labeller = labeller(.multi_line = FALSE)) +
  theme_bw(base_size = 14) +
  labs(x = 'parameter', y = 'estimate') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# read in musse model
fit_musse_no_se <- readRDS('data/sequencing_rpoB/processed/transition_rates/asv_musse_no_se.rds')
fit_musse_no_se$par %>% round(., 3) %>% sort()

# compare rates between musse and secsse
d_musse_diversitree <- data.frame(diversitree = fit_musse_no_se$par) %>%
  rownames_to_column(var = 'param')
d_musse_compare <- filter(d_musse_rates, samp_frac == 1) %>%
  filter(aic == min(aic)) %>%
  filter(!is.na(estimate) & estimate > 0) %>%
  mutate(param = gsub('A', '', param)) %>%
  left_join(., d_musse_diversitree)

ggplot(d_musse_compare, aes(diversitree, estimate)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggrepel::geom_label_repel(aes(label = param), force_pull = 0.001, box.padding = 1, point.padding = 1, force = 0.1, hjust = 0.75) +
  geom_point(aes(col = type), size = 3) +
  theme_bw(base_size = 14) +
  labs(x = 'Estimate from diversitree',
       y = 'Estimate from secsse',
       title = 'Comparison of MuSSE models fit using secsse and diversitree.')


#-----------------------------#
# look at the muhisse models ####
#-----------------------------#

files_muhisse1 <- files[str_detect(files, 'muhisseSSonly')]

# calculate initial values
d_muhisse1_aic <- map_df(files_muhisse1, calc_aic) %>%
  arrange(samp_frac, aic) %>%
  mutate(model = 'muhisse1')

# plot these
ggplot(d_muhisse1_aic, aes(as.character(samp_frac), aic)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(aes(label = label_wrap_mod(initial_values, 20)), force_pull = 0.5, box.padding = 0.5) +
  theme_bw()

# mk transitions not always the best
# i.e. initial values do matter (QUITE A LOT)

# look at variation in rates for each initial rate values
d_muhisse1_rates <- map_df(files_muhisse1, grab_params) %>%
  left_join(., d_muhisse1_aic) %>%
  group_by(param) %>%
  mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(estimate))

filter(d_muhisse1_rates, estimate >0) %>%
  ggplot(., aes(fct_reorder(param, estimate, .fun = mean), estimate)) +
  geom_point(aes(col = initial_values)) +
  facet_wrap(~type+ samp_frac, scale = 'free') +
  theme_bw(base_size = 12) +
  labs(x = 'parameter', y = 'estimate') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-----------------------------#
# look at the CTD2 models ####
#-----------------------------#

files_ctd2 <- files[str_detect(files, 'ctd2')]

# calculate initial values
d_ctd2_aic <- map_df(files_ctd2, calc_aic) %>%
  arrange(samp_frac, aic) %>%
  mutate(model = 'ctd2')

# plot these
ggplot(d_ctd2_aic, aes(as.character(samp_frac), aic)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(aes(label = label_wrap_mod(initial_values, 20)), force_pull = 0.5, box.padding = 0.5) +
  theme_bw()

# mk transitions not always the best
# i.e. initial values do matter (QUITE A LOT)

# look at variation in rates for each initial rate values
d_ctd2_rates <- map_df(files_ctd2, grab_params) %>%
  left_join(., d_ctd2_aic) %>%
  group_by(param) %>%
  mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(estimate))

filter(d_ctd2_rates, estimate >0) %>%
  ggplot(., aes(fct_reorder(param, estimate, .fun = mean), estimate)) +
  geom_point(aes(col = initial_values)) +
  facet_wrap(~type+ samp_frac, scale = 'free') +
  theme_bw(base_size = 12) +
  labs(x = 'parameter', y = 'estimate') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-----------------------------#
# look at the CTD3 models ####
#-----------------------------#

files_ctd3 <- files[str_detect(files, 'ctd3')]

# calculate initial values
d_ctd3_aic <- map_df(files_ctd3, calc_aic) %>%
  arrange(samp_frac, aic) %>%
  mutate(model = 'ctd3')

# plot these
ggplot(d_ctd3_aic, aes(as.character(samp_frac), aic)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(aes(label = label_wrap_mod(initial_values, 20)), force_pull = 0.5, box.padding = 0.5) +
  theme_bw()

# mk transitions not always the best
# i.e. initial values do matter (QUITE A LOT)

# look at variation in rates for each initial rate values
d_ctd3_rates <- map_df(files_ctd3, grab_params) %>%
  left_join(., d_ctd3_aic) %>%
  group_by(param) %>%
  mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(estimate))

filter(d_ctd3_rates, estimate >0) %>%
  ggplot(., aes(fct_reorder(param, estimate, .fun = mean), estimate)) +
  geom_point(aes(col = initial_values)) +
  facet_wrap(~type+ samp_frac, scale = 'free') +
  theme_bw(base_size = 12) +
  labs(x = 'parameter', y = 'estimate') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-----------------------------#
# look at the muhisse full models ####
#-----------------------------#

files_muhisse2 <- files[str_detect(files, 'muhissefull')]

# calculate initial values
d_muhisse2_aic <- map_df(files_muhisse2, calc_aic) %>%
  arrange(samp_frac, aic) %>%
  mutate(model = 'muhisse2')

# plot these
ggplot(d_muhisse2_aic, aes(as.character(samp_frac), aic)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(aes(label = label_wrap_mod(initial_values, 20)), force_pull = 0.5, box.padding = 0.5) +
  theme_bw()

# mk transitions not always the best
# i.e. initial values do matter (QUITE A LOT)

# look at variation in rates for each initial rate values
d_muhisse2_rates <- map_df(files_muhisse2, grab_params) %>%
  left_join(., d_muhisse2_aic) %>%
  group_by(param) %>%
  mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(estimate))

filter(d_muhisse2_rates, estimate >0) %>%
  ggplot(., aes(fct_reorder(param, estimate, .fun = mean), estimate)) +
  geom_point(aes(col = initial_values)) +
  facet_wrap(~type+ samp_frac, scale = 'free') +
  theme_bw(base_size = 12) +
  labs(x = 'parameter', y = 'estimate') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

