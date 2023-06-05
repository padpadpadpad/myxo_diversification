# check which model is best for sriswasdi

# load packages
librarian::shelf(tidyverse, flextable)
tidyverse_conflicts()

# list model output files
files <- list.files('data/sriswasdi_data/results', full.names = TRUE, recursive = TRUE)

files[1]

# write function to extract logLik from the models and create a dataframe
get_loglik <- function(secsse_file){
  temp <- readRDS(secsse_file)
  temp <- data.frame(model_name = basename(tools::file_path_sans_ext(secsse_file)),
             loglik = temp$mod$ML,
             n_params = temp$n_params)
  return(temp)
}

# read in musse models ####
fits_musse <- str_subset(files, 'musse') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 1, 2, 3, 5, 6 
# AIC = -1427.806

# read in muhisse full models ####
fits_muhisse <- str_subset(files, 'muhissefull') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 1 (and run 3)
# AIC = -2317.801

# read in muhisse 2 ####
fits_muhisse2 <- str_subset(files, 'muhisse2') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 3 and 6
# AIC = -2321.799

# read in ctd2 models
fits_ctd2 <- str_subset(files, 'ctd2') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 6: AIC = -2239.474

# read in those four models and do model comparison using AIC scores
best_ctd2 <- readRDS('data/sriswasdi_data/results/seccse_ctd2_sampfrac1_run6.rds')
best_muhisse1 <- readRDS('data/sriswasdi_data/results/seccse_muhissefull_sampfrac1_run1.rds')
best_muhisse2 <- readRDS('data/sriswasdi_data/results/seccse_muhisse2_sampfrac1_run3.rds')
best_musse <- readRDS('data/sriswasdi_data/results/seccse_musse_sampfrac1_run1.rds')

final_models <- c('data/sriswasdi_data/results/seccse_ctd2_sampfrac1_run6.rds', 
                  'data/sriswasdi_data/results/seccse_muhissefull_sampfrac1_run1.rds',
                  'data/sriswasdi_data/results/seccse_muhisse2_sampfrac1_run3.rds',
                  'data/sriswasdi_data/results/seccse_musse_sampfrac1_run1.rds')

best_ctd2$mod$ML
best_muhisse1$mod$ML
best_muhisse2$mod$ML
best_musse$mod$ML

model_table <- final_models %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         weights = MuMIn::Weights(aic) %>% round(2)) %>%
  arrange(-weights) %>%
  mutate(model = case_when(model == 'ctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhissefull' ~ 'MuHiSSE v1',
                           model == 'muhisse2' ~ 'MuHiSSE v2'))

# make table
table <- select(model_table, model, n_params, loglik, aic, weights) %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    n_params = 'Number of estimated parameters',
                    loglik = 'Log Likelihood',
                    aic = 'AIC',
                    weights = "AIC weight") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 
  
table

save_as_image(table, 'plots/manuscript_plots/sriswasdi_table.png')

# the best model is muhisse2

# look at the parameters
best_muhisse2$mod$MLpars
# BIG NUMBERS
best_muhisse2$mod$MLpars 


# look at the
