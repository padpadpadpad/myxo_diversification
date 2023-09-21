# analyse age of transitions

# load in packages
library(phytools)
library(progress)
library(here)
library(tidyverse)
library(castor)
library(tidytree)

# load in actual tree
tree <- read.tree(here('data/sequencing_rpoB/bamm/rerooted-pruned-chronopl10.tre'))

# load in age of transition data
res <- readRDS(here('data/sequencing_rpoB/processed/transition_rates/time_of_transition.rds'))

# edit tree labels
d_labels <- data.frame(tip_label = tree$tip.label) %>%
  separate(., tip_label, c('part1', 'part2', 'part3'), sep = '_', remove = FALSE) %>%
  unite('tip_label_new', c(part1, part2), sep = '_')

tree$tip.label <- d_labels$tip_label_new

str(tree$edge)

head(tree$edge)

root_node <- castor::find_root(tree)

# calculate distances of all tips and nodes to the root
distances <- as_tibble(tree) %>%
  # think this returns the right order
  mutate(dist_to_root = castor::get_all_distances_to_root(tree))

df <- as_tibble(tree)

# results
res2 <- separate(res, mapped_edge, c('node1', 'node2'), sep = ',', remove = FALSE) %>%
  mutate(across(starts_with('node'), as.numeric))

# combine the distances to the root to the original data frame
res2 <- left_join(res2, select(distances, node1 = node, dist_to_root1 = dist_to_root)) %>%
  left_join(., select(distances, node2 = node, dist_to_root2 = dist_to_root))

head(res2)

# as the transition does not happen at a node, but between two nodes, calculate the distance to root as the node 1 + 0.5*(dist_to_root2 - dist_to_root)
res2 <- mutate(res2, dist_to_root = dist_to_root1 + 0.5*(dist_to_root2 - dist_to_root1))

# add column to define the transition
res2 <- mutate(res2, transition = paste(a, b, sep = ' -> '))

# look at when transitions occur in the evolutionary history of myxococcus
ggplot(res2, aes(dist_to_root)) +
  geom_histogram(fill = 'white', col = 'black') +
  facet_wrap(~transition, scales = 'free_y')

ggplot(res2, aes(dist_to_root1)) +
  geom_density() +
  facet_wrap(~transition, scales = 'free_y')
# there do not seem to be obvious differences between any of the transitions here

# calculate summary of average age of transition
res2_summary <- group_by(res2, transition) %>%
  tidybayes::median_qi(dist_to_root)

# look at common nodes for transitions
res2_node_summary <- group_by(res2, mapped_edge, node1, node2, a, b, transition, dist_to_root) %>%
  tally() %>%
  mutate(prop = n/1000)

res2_node_summary2 <- group_by(res2, mapped_edge, node1, node2, dist_to_root) %>%
  tally()

# plot distance to root and how common the transition is
ggplot(res2_node_summary, aes(dist_to_root, prop)) +
  geom_point(shape = 21, fill = 'white') +
  facet_wrap(~transition) +
  theme_bw() +
  theme(strip.text = element_text(size = 8))

ggsave('plots/sequencing_rpoB/analyses/common_transition_nodes_per_transition.png', last_plot(), height = 8, width = 12)

# plot distance to root and how common the transition is
ggplot(res2_node_summary2, aes(dist_to_root, n)) +
  geom_point(shape = 21, fill = 'white') +
  theme_bw() +
  theme(strip.text = element_text(size = 8))

ggsave('plots/sequencing_rpoB/analyses/common_transition_nodes.png', last_plot(), height = 8, width = 12)