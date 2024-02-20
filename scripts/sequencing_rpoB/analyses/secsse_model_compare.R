# lets look at which hidden model is best

# load packages
librarian::shelf(secsse, tidyverse, flextable, patchwork, officer)

# write function to extract logLik from the models and create a dataframe
get_loglik <- function(secsse_file){
  temp <- readRDS(secsse_file)
  temp <- data.frame(model_name = basename(tools::file_path_sans_ext(secsse_file)),
                     loglik = temp$mod$ML,
                     n_params = temp$n_params,
                     conv = temp$mod$conv)
  return(temp)
}

# write function to extract the parameters from the model
get_pars <- function(secsse_file){
  temp <- readRDS(secsse_file)
  
  temp_speciation <- data.frame(val = temp$mod$MLpars[[1]]) %>%
    rownames_to_column(var = 'param') %>%
    group_by(val) %>%
    slice_head(n=1) %>%
    mutate(param = paste('lambda', param, sep = '_'))
    
  temp_extinction <- data.frame(val = temp$mod$MLpars[[2]]) %>%
    rownames_to_column(var = 'param') %>%
    slice_head(n = 1) %>%
    mutate(param = paste('mu', param, sep = '_'))
  
  temp_transition <- data.frame(temp$mod$MLpars[[3]]) %>% 
    rownames_to_column(var = 'from') %>%
    pivot_longer(cols = -from, names_to = 'to', values_to = 'val') %>%
    mutate(to = gsub('X', '', to),
           param = ifelse(substr(to, 2, 2) == substr(from, 2, 2), paste('q_', substr(from, 1, 1), substr(to, 1, 1), sep = ''), paste('q_', substr(from, 2, 2), substr(from, 2, 2), sep = ''))) %>%
    filter(val > 0) %>%
    select(param, val) %>%
    distinct()
  
  temp <- bind_rows(temp_speciation, temp_extinction, temp_transition) %>%
    mutate(model_name = basename(tools::file_path_sans_ext(secsse_file)))

  return(temp)
  
}

#--------------------------------#
# first for the ASV best tree ####
#--------------------------------#

# read in all the files
files <- list.files('data/sequencing_rpoB/processed/secsse/results/asv', full.names = TRUE)

# read in all files
fits <- files %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         aic = round(aic, 2),
         # extract text between the first and second underscore
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         # extract text between the second and third underscore
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name))

# check the number of models per samp frac and model that have the same aic
check <- fits %>%
  group_by(samp_frac, model, aic) %>%
  tally()

# how often was the lowest AIC score found in each combination
check2 <- slice_min(check, n = 1, order_by = aic) %>%
  mutate(., best = 'yes')
# some were only found once

# filter the files to only keep the lowest aic score for each model
best_models <- fits %>%
  left_join(check2, by = c('samp_frac', 'model', 'aic')) %>%
  filter(best == 'yes') %>%
  group_by(samp_frac, model) %>%
  slice_min(order_by = aic, n = 1, with_ties = FALSE) %>%
  ungroup()

# a few of the best models were only found one time in terms of the start values, but not sure what more can be done about that, this may indicate local optima, but we have followed the best practice approaches

# find the files with the lowest aic per group and keep them to read back in
files_best <- files[str_detect(files, paste(best_models$model_name, collapse = '|'))]

# read in best files only
fits <- files_best %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         aic = round(aic, 2),
         # extract text between the first and second underscore
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         # extract text between the second and third underscore
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = parse_number(samp_frac))

# make table for sampling fraction 0.5
model_table <- fits %>%
  filter(samp_frac == '0.5') %>%
  mutate(., weights = MuMIn::Weights(aic) %>% round(2)) %>%
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
save_as_docx(table, path ='plots/manuscript_plots/secsse_table.docx')

model_table_all <- fits %>%
  filter(samp_frac != '0.5') %>%
  group_by(samp_frac) %>%
  mutate(., weights = MuMIn::Weights(aic) %>% round(2)) %>%
  arrange(samp_frac, aic) %>%
  ungroup() %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

# make table
table <- select(model_table_all, samp_frac, model, n_params, loglik, aic, weights) %>%
  mutate(samp_frac = as.character(samp_frac), 
         across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(samp_frac = 'Sampled fraction',
                    model = "Model",
                    n_params = 'Number of estimated parameters',
                    loglik = 'Log Likelihood',
                    aic = 'AIC',
                    weights = "AIC weight") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  hline(i = c(5, 10, 15), border = fp_border(color = "grey"), part = 'body') %>%
  merge_v(~samp_frac) %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  fix_border_issues(part = 'all')

print(table)
plot(table)

save_as_image(table, 'plots/manuscript_plots/secsse_supp_table.png', webshot = 'webshot')

# read in parameter values of all the best models
pars <- files_best %>%
  map(., get_pars) %>%
  list_rbind() %>%
  mutate(., model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = parse_number(samp_frac)) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

p_extinction <- filter(pars, str_detect(param, 'mu')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  ggplot(aes(samp_frac, val)) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*sd(x)/sqrt(length(x))}, fun.max = function(x){mean(x) + 2*sd(x)/sqrt(length(x))}) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Extinction rate (mu)',
       title = '(b) Extinction rate')

p_speciation <- filter(pars, str_detect(param, 'lambda')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  ggplot(aes(samp_frac, val)) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*sd(x)/sqrt(length(x))}, fun.max = function(x){mean(x) + 2*sd(x)/sqrt(length(x))}) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Speciation rate (lambda)',
       title = '(a) Speciation rate') 

p_speciation

p_transition <- filter(pars, str_detect(param, 'q')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  filter(is.na(parse_number(param))) %>%
  ggplot(aes(samp_frac, val)) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*sd(x)/sqrt(length(x))}, fun.max = function(x){mean(x) + 2*sd(x)/sqrt(length(x))}) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_dodge(width = 0.5), show.legend = FALSE) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Hidden transition rates',
       title = '(c) Transition rates')

p_speciation + p_extinction + p_transition + plot_layout(guides = 'collect')

# save as image
ggsave('plots/manuscript_plots/secsse_params.png', width = 12, height = 3.5)

#------------------------------------------------#
# read in models for the bootstrap replicates ####
#------------------------------------------------#

# read in all the files
files <- list.files('data/sequencing_rpoB/processed/secsse/results/bootstraps', full.names = TRUE)

# read in all files
fits <- files %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         aic = round(aic, 2),
         # extract text between the first and second underscore
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         # extract text between the second and third underscore
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name),
         # extract number after the last underscore
         boot = sub(".*_(.*?)$", "\\1", model_name))

# check the number of models per samp frac and model that have the same aic
check <- fits %>%
  group_by(samp_frac, boot, model, aic) %>%
  tally()

# how often was the lowest AIC score found in each combination
check2 <- slice_min(check, n = 1, order_by = aic) %>%
  mutate(., best = 'yes')
# some were only found once

# filter the files to only keep the lowest aic score for each model
best_models <- fits %>%
  left_join(check2, by = c('samp_frac', 'boot', 'model', 'aic')) %>%
  filter(best == 'yes') %>%
  group_by(boot, samp_frac, model) %>%
  slice_min(order_by = aic, n = 1, with_ties = FALSE) %>%
  ungroup()

# a few of the best models were only found one time in terms of the start values, but not sure what more can be done about that, this may indicate local optima, but we have followed the best practice approaches

# find the files with the lowest aic per group and keep them to read back in
files_best <- files[str_detect(files, paste(best_models$model_name, collapse = '|'))]

# read in best files only
fits <- files_best %>%
  map(., get_loglik) %>%
  list_rbind() %>%
  mutate(., aic = (2*n_params) - 2*loglik,
         aic = round(aic, 2),
         # extract text between the first and second underscore
         model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         # extract text between the second and third underscore
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = parse_number(samp_frac),
         boot = sub(".*_(.*?)$", "\\1", model_name),
         boot = paste('ASV boot', parse_number(boot), sep = ' ')) %>%
  group_by(boot, samp_frac) %>%
  mutate(., weights = MuMIn::Weights(aic) %>% round(2)) %>%
  arrange(aic) %>%
  mutate(order = 1:n()) %>%
  ungroup()

# make table for sampling fraction 0.5
model_table <- fits %>%
  group_by(samp_frac, boot) %>%
  slice_max(weights, n = 1) %>%
  arrange(boot) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

# make table
table <- select(model_table, boot, samp_frac, model, n_params, loglik, aic, weights) %>%
  mutate(across(where(is.numeric), \(x) round(x,2))) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(model = "Model",
                    samp_frac = 'Sampled fraction',
                    boot = 'Bootstrap',
                    n_params = 'Number of estimated parameters',
                    loglik = 'Log Likelihood',
                    aic = 'AIC',
                    weights = "AIC weight") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  merge_v(~boot) %>%
  hline(i = c(5, 10, 15, 20, 25, 30, 35, 40), border = fp_border(color = "grey"), part = 'body') %>%
  valign(valign = 'top', j = 1, part = 'body') %>%
  fix_border_issues(part = 'all') %>%
  autofit() 

table

save_as_image(table, 'plots/manuscript_plots/secsse_table_boots.png')
save_as_docx(table, path ='plots/manuscript_plots/secsse_table.docx')

# read in parameter values of all the best models
pars <- files_best %>%
  map(., get_pars) %>%
  list_rbind() %>%
  mutate(., model = sub("^[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = sub("^[^_]*_[^_]*_(.*?)_.*$", "\\1", model_name),
         samp_frac = parse_number(samp_frac),
         boot = sub(".*_(.*?)$", "\\1", model_name),
         boot = paste('ASV boot', parse_number(boot), sep = ' ')) %>%
  mutate(model = case_when(model == 'muctd2' ~ 'CTD2',
                           model == 'musse' ~ 'MuSSE',
                           model == 'muhisseSSonly' ~ 'MuHiSSE',
                           model == 'muctd3' ~ 'CTD3',
                           model == 'muctd4' ~ 'CTD4'))

p_extinction <- filter(pars, str_detect(param, 'mu')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  ggplot(aes(samp_frac, val)) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05)) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*sd(x)/sqrt(length(x))}, fun.max = function(x){mean(x) + 2*sd(x)/sqrt(length(x))}) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Extinction rate (mu)',
       title = '(b) Extinction rate')

p_speciation <- filter(pars, str_detect(param, 'lambda')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  ggplot(aes(samp_frac, val)) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05)) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*sd(x)/sqrt(length(x))}, fun.max = function(x){mean(x) + 2*sd(x)/sqrt(length(x))}) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Speciation rate (lambda)',
       title = '(a) Speciation rate') 

p_speciation

p_transition <- filter(pars, str_detect(param, 'q')) %>%
  mutate(samp_frac = as.character(samp_frac)) %>%
  filter(is.na(parse_number(param))) %>%
  ggplot(aes(samp_frac, val)) +
  geom_point(aes(group = model, col = model), shape = 21, fill = 'white', position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05), show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "point", col = "black", size = 4) +
  stat_summary(fun = "mean", fun.min = function(x){mean(x) - 2*(sd(x)/sqrt(length(x)))}, fun.max = function(x){mean(x) + 2*(sd(x)/sqrt(length(x)))}) +
  theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(x = 'Sampled fraction',
       y = 'Hidden transition rates',
       title = '(c) Transition rates')

p_speciation + p_extinction + p_transition + plot_layout(guides = 'collect')

# save as image
ggsave('plots/manuscript_plots/secsse_params_boots.png', width = 12, height = 3.5)
