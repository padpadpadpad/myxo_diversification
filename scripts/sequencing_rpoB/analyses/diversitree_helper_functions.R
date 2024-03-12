# write function to plot diversitree model

# function for getting a data frame from a diversitree object
get_diversitree_df <- function(div_obj, trait_vec, replace_vec){
  
  if(is.null(div_obj$par.full)){div_obj$par.full <- div_obj$par}
  
  temp <- tibble(param = names(div_obj$par.full)) %>%
    mutate(state_1_num = as.numeric(substr(param, 2,2)),
           state_2_num = as.numeric(substr(param, 3,3)),
           transition_rate = unlist(div_obj$par.full),
           state_1 = stringi::stri_replace_all_regex(state_1_num, pattern = trait_vec, replacement = replace_vec, vectorize=FALSE),
           state_2 = stringi::stri_replace_all_regex(state_2_num, pattern = trait_vec, replacement = replace_vec, vectorize=FALSE)) %>%
    select(param, state_1, state_2, state_1_num, state_2_num, transition_rate) %>%
    mutate(free_param = ifelse(param %in% names(div_obj$par), 'yes', 'no'),
           num_params = length(div_obj$par))
  
  return(temp)
  
}

# function for creating a simple plot and extracting a legend from it

make_legend <- function(d_habitat, col_vec){
  
  p <- ggplot(d_habitat, aes(habitat_preference, prop, col = habitat_preference)) +
    geom_point() +
    scale_color_manual('Biome preference', values = col_vec, labels = c('freshwater + land generalist', 'freshwater specialist', 'marine generalist', 'marine specialist', 'land specialist')) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=5)))
  p <- cowplot::get_legend(p)
  return(p)
}

# function for plotting diversitree dataset
make_network_diversitree <- function(diversitree_df, d_habpref_summary, cols_hab){
  
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
                      edge_width = transition_rate,
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
    scale_size(range = c(5,20)) +
    scale_edge_width(range = c(0.5, 2)) +
    scale_color_manual('Biome preference', values = cols_hab) +
    scale_fill_manual('Biome preference', values = cols_hab)
  
  # grab data for points
  point_data <- p$data %>%
    select(x, y, name) %>%
    mutate(nudge_x = ifelse(x < 0, -0.2, 0.2),
           nudge_y = ifelse(y < 0, -0.2, 0.2))
  
  p <- p +
    xlim(c(min(point_data$x) + min(point_data$nudge_x)), max(point_data$x) + max(point_data$nudge_x)) +
    ylim(c(min(point_data$y) + min(point_data$nudge_y)), max(point_data$y) + max(point_data$nudge_y)) +
    theme(legend.position = 'none',
          panel.background = element_rect(fill = 'white', colour = 'white'))
  
  return(p)
  
}

# function for making a source sink table

make_source_sink_table <- function(diversitree_df){
  
  temp <- select(diversitree_df, away = state_1, into = state_2, transition_rate) %>%
    pivot_longer(cols = c(away, into), names_to = 'direction', values_to = 'habitat_preference') %>%
    group_by(habitat_preference, direction) %>%
    summarise(total_rate = sum(transition_rate), .groups = 'drop') %>%
    pivot_wider(names_from = direction, values_from = total_rate) %>%
    mutate(source_sink1 = away / into,
           source_sink2 = into - away) %>%
    select(habitat_preference, away, into, source_sink1) %>%
    mutate(across(away:source_sink1, ~round(.x, 2)),
           habitat_preference = gsub(':', ' + ', habitat_preference),
           habitat_preference = gsub('_', ' ', habitat_preference)) %>%
    arrange(desc(source_sink1)) %>%
    flextable(.) %>%
    set_header_labels(habitat_preference = 'biome preference',
                      source_sink1 = 'source sink ratio') %>%
    align(align = 'center', part = 'all') %>%
    align(align = 'left', part = 'body', j = 1) %>%
    font(fontname = 'Times', part = 'all') %>%
    fontsize(size = 16, part = 'all')
  
  return(temp)
}

# function to make ggplot2 friendly ltt plot
plot_ltt <- function(tree, log = 'N'){
  
  # create lineage through time plot
  d_ltt <- ape::ltt.plot.coords(tree) %>%
    data.frame() %>%
    mutate(time2 = time + 1)
  
  # create lineage through time plot
  # plot N or log(N) based on the argument log = Y or N
  if(log == 'Y'){
    p2 <- ggplot(d_ltt, aes(time, log(N))) +
      geom_line() +
      theme_bw(base_size = 12) +
      labs(x = 'Relative time from present',
           y = 'Log number of lineages')
  } else {
    p2 <- ggplot(d_ltt, aes(time2, N)) +
    geom_line() +
    theme_bw(base_size = 12) +
    labs(x = 'Relative time',
         y = 'Number of lineages')
  }
  
  return(p2)
}

# function for wrapping text in a label
label_wrap2 <- function(x, width){
  unlist(lapply(strwrap(x, width = width, simplify = FALSE), 
                paste0, collapse = "\n"))
}


