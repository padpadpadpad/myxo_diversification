# visualise old vs new LTT plot

librarian::shelf(ape, tidyverse, patchwork)

# function to make ggplot2 friendly ltt plot
plot_ltt <- function(tree, log = 'N'){
  
  # create lineage through time plot
  d_ltt <- ape::ltt.plot.coords(tree) %>%
    data.frame() %>%
    mutate(time2 = time + 1)
  
  # create lineage through time plot
  # plot N or log(N) based on the argument log = Y or N
  if(log == 'Y'){
    p2 <- ggplot(d_ltt, aes(time2, log(N))) +
      geom_line() +
      theme_bw(base_size = 12) +
      labs(x = 'Relative time',
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

# load in trees
tree1 <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_chronopl10.tre')
tree2 <- read.tree('data/sequencing_rpoB/raxml/trees/myxo_asv/myxo_asv_treepl_cv_node_labels.tre')

p1 <- plot_ltt(tree1, log = 'Y')
p2 <- plot_ltt(tree2, log = 'Y')

p1 + labs(title = '(a) LTT from chronopl method') + p2 + labs(title = '(b) LTT from treepl method')

ggsave('plots/manuscript_plots/phylogeny_checks/ltt_comparison.png', width = 10, height = 4)
