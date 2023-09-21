# lets look at which hidden model is best

# load packages
librarian::shelf(secsse, tidyverse, flextable)

# write function to extract logLik from the models and create a dataframe
get_loglik <- function(secsse_file){
  temp <- readRDS(secsse_file)
  temp <- data.frame(model_name = basename(tools::file_path_sans_ext(secsse_file)),
                     loglik = temp$mod$ML,
                     n_params = temp$n_params)
  return(temp)
}

#----------------------#
# sample fraction 1 ####
#----------------------#

# list all model files
files <- list.files('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2', full.names = TRUE)

# read in musse models
fits_musse <- str_subset(files, 'musse') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 1, 2, 3, 4, 5, 6 
# AIC = 3451.294

# read in muhisse full models
fits_muhisse <- str_subset(files, 'muhisse') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 4, 5, 6
# AIC = 2704.353

# read in ctd2
fits_ctd2 <- str_subset(files, 'ctd2') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs all of them
# AIC = 2762.069

# read in ctd3 models
fits_ctd3 <- str_subset(files, 'ctd3') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is runs 1, 2, 3, 4, 6 - AIC = 2647.654

# read in ctd4 models
fits_ctd4 <- str_subset(files, 'ctd4') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model are run 1 and run 4 - AIC = 2613.6

# read in ctd5 models
fits_ctd5 <- str_subset(files, 'ctd5') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run all of them: AIC = 3083.513

# read in those four models and do model comparison using AIC scores
best_ctd2 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd2_sampfrac1_run1_v2.rds')
best_muhisse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muhisseSSonly_sampfrac1_run4_v2.rds')
best_ctd3 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd3_sampfrac1_run1_v2.rds')
best_musse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_musse_sampfrac1_run1v2.rds')
best_ctd4 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd4_sampfrac1_run1_v2.rds')
#best_ctd5 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd4_sampfrac1_run1_v2.rds')

final_models <- c('data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd2_sampfrac1_run1_v2.rds',
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muhisseSSonly_sampfrac1_run4_v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd3_sampfrac1_run1_v2.rds', 'data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_musse_sampfrac1_run1v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_1/v2/seccse_muctd4_sampfrac1_run1_v2.rds')

best_musse$mod$ML
best_muhisse$mod$ML
best_ctd2$mod$ML
best_ctd3$mod$ML
best_ctd4$mod$ML
best_ctd5$mod$ML

best_ctd3$mod$MLpars[[1]]
best_ctd4$mod$MLpars[[1]]

model_table <- final_models %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         weights = MuMIn::Weights(aic) %>% round(2)) %>%
  arrange(aic) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

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

save_as_image(table, 'plots/manuscript_plots/secsse_table.png')

#-------------------------#
# sample fraction 0.75 ####
#-------------------------#

# list all model files
files <- list.files('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2', full.names = TRUE)

# read in musse models
fits_musse <- str_subset(files, 'musse') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 7, 8, 9, 10, 11, 12 
# AIC = 3389.113

# read in muhisse full models
fits_muhisse <- str_subset(files, 'muhisseSS') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 7, 8, 10, 11, 12 
# AIC = 2642.122

# read in ctd2
fits_ctd2 <- str_subset(files, 'ctd2') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 7, 8, 9, 10, 11, 12 
# AIC = 2699.069

# read in ctd3 models
fits_ctd3 <- str_subset(files, 'ctd3') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 9
# AIC = 2590.748

# read in ctd4 models
fits_ctd4 <- str_subset(files, 'ctd4') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 12
# AIC = 2564.462

# read in those four models and do model comparison using AIC scores
best_musse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_musse_sampfrac0.75_run10v2.rds')
best_muhisse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muhisseSSonly_sampfrac0.75_run10_v2.rds')
best_ctd2 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd2_sampfrac0.75_run7_v2.rds')
best_ctd3 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd3_sampfrac0.75_run9_v2.rds')
best_ctd4 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd4_sampfrac0.75_run12.rds')
#best_ctd5 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd3_sampfrac1_run1.rds')

final_models_0.75 <- c('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_musse_sampfrac0.75_run10v2.rds',
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muhisseSSonly_sampfrac0.75_run10_v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd2_sampfrac0.75_run7_v2.rds',
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd3_sampfrac0.75_run9_v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.75/v2/seccse_muctd4_sampfrac0.75_run12_v2.rds')

best_ctd4$mod$MLpars

#------------------------#
# sample fraction 0.5 ####
#------------------------#

# list all model files
files <- list.files('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2', full.names = TRUE)

# read in musse models
fits_musse <- str_subset(files, 'musse') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 13, 14, 15, 16, 17, 18 
# AIC = 3333.751

# read in muhisse full models
fits_muhisse <- str_subset(files, 'muhisseSS') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 14, 15, 18
# AIC = 2565.014

# read in ctd2
fits_ctd2 <- str_subset(files, 'ctd2') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best models are runs 13, 14, 15, 16, 17, 18  
# AIC = 2621.852

# read in ctd3 models
fits_ctd3 <- str_subset(files, 'ctd3') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 18
# AIC = 2517.092

# read in ctd4 models
fits_ctd4 <- str_subset(files, 'ctd4') %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik)
# best model is run 15
# AIC = 2520.147

# read in those four models and do model comparison using AIC scores
best_musse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_musse_sampfrac0.5_run13v2.rds')
best_muhisse <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muhisseSSonly_sampfrac0.5_run18_v2.rds')
best_ctd2 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd2_sampfrac0.5_run13_v2.rds')
best_ctd3 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd3_sampfrac0.5_run18_v2.rds')
best_ctd4 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd2_sampfrac0.5_run15_v2.rds')
#best_ctd5 <- readRDS('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd3_sampfrac1_run1.rds')

final_models_0.5 <- c('data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_musse_sampfrac0.5_run13v2.rds',
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muhisseSSonly_sampfrac0.5_run18_v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd2_sampfrac0.5_run13_v2.rds',
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd3_sampfrac0.5_run18_v2.rds', 
                  'data/sequencing_rpoB/processed/secsse/results/samp_frac_0.5/v2/seccse_muctd4_sampfrac0.5_run15_v2.rds')

best_ctd2$mod$MLpars

# create summary table of the sampling fractions for 0.75 and 0.5
model_table_0.75 <- final_models_0.75 %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         weights = MuMIn::Weights(aic) %>% round(2),
         sample_fraction = 0.75) %>%
  arrange(aic) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))
model_table_0.5 <- final_models_0.5 %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         weights = MuMIn::Weights(aic) %>% round(2),
         sample_fraction = 0.5) %>%
  arrange(aic) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

model_table_all <- bind_rows(model_table_0.75, model_table_0.5)

model_table_all

# make table
table <- select(model_table_all, sample_fraction, model, n_params, loglik, aic, weights) %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(sample_fraction = 'Sampled fraction',
                    model = "Model",
                    n_params = 'Number of estimated parameters',
                    loglik = 'Log Likelihood',
                    aic = 'AIC',
                    weights = "AIC weight") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  merge_v(~sample_fraction) %>%
  hline(i = c(5), border = fp_border_default()) %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  fix_border_issues() %>%
  autofit() 

table

save_as_image(table, 'plots/manuscript_plots/secsse_supp_table.png')
