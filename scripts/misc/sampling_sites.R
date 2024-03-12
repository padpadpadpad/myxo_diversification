#---------------------------------#
# create map of sampling sites ####
#---------------------------------#

# load in packages
library(tidyverse)
library(ggmap)
library(patchwork)
library(flextable)

# load in data
d_locations <- readxl::read_excel('data/Supplemental Table 1 - Sampling.xlsx') %>%
  select(`sample nr`, location, coordinates, site) %>%
  mutate(id = gsub('\\.0', '', `sample nr`)) %>%
  select(-`sample nr`)

d_samples <- readxl::read_excel('data/Supplemental Table 1 - Sampling.xlsx') %>%
  select(`sample type`) %>%
  group_by_all() %>%
  tally()

# load in metadata
d_meta <- read.csv('data/metadata.csv') %>%
  select(-location)

d <- left_join(d_meta, d_locations) 

# make Table S1
d_table <- separate(d, coordinates, c('lat', 'lon'), sep = ',') %>%
  mutate(across(c(lat, lon), as.numeric),
         predefined_habitat = case_when(habitat_group_16s == 'woodland_oak' ~ 'soil underneath oak',
                                        habitat_group_16s == 'woodland_pine' ~ 'soil underneath monterey pine',
                                        habitat_group_16s == 'field_wheat' ~ 'wheat field soil',
                                        habitat_group_16s == 'river' ~ 'riverbed sediment',
                                        habitat_group_16s == 'reservoir' ~ 'reservoir/lake sediment',
                                        habitat_group_16s == 'pasture' ~ 'pasture soil',
                                        habitat_group_16s == 'marine mud_full saline' ~ 'marine sediment',
                                        habitat_group_16s == 'beach_seaweed' ~ 'beachcast seaweed',
                                        habitat_group_16s == 'rock_samphire' ~ 'rock samphire rhizosphere',
                                        habitat_group_16s == 'woodland_oak' ~ 'low subtidal sand',
                                        habitat_group_16s == 'estuarine mud_full saline' ~ 'estuarine sediment (close to full salinity)',
                                        habitat_group_16s == 'estuarine mud_low polyhaline' ~ 'estuarine sediment (polyhaline/mesohaline)',
                                        habitat_group_16s == 'estuarine mud_oligohaline' ~ 'estuarine sediment (oligohaline)',
                                        habitat_group_16s == 'beach_subtidal' ~ 'low subtidal beach sand',
                                        habitat_group_16s == 'beach_supratidal' ~ 'high supratidal beach sand',
                                        habitat_group_16s == 'thrift_rhizosphere' ~ 'thrift rhizosphere'),
         sequenced_16s = 'yes',
         sequenced_rpob = ifelse(habitat_group_16s %in% c('beach_supratidal', 'thrift_rhizosphere', 'estuarine mud_low polyhaline'), 'no', 'yes'))

# read in cluster assignments
clusters <- read.csv('data/sequencing_16S/sample_cluster_assignments.csv') %>%
  select(X = sample, biome = medoid_nbclust) %>%
  mutate(biome = case_when(biome == '1' ~ 'land',
                           biome == '2' ~ 'freshwater',
                           biome == '3' ~ 'marine'),
         biome = ifelse(X == 's46', 'marine', biome))

d_table <- left_join(d_table, clusters)

# grab info to paste into ena
d_table %>%
  dplyr::arrange(id) %>%
  select(location, site, lat, lon, predefined_habitat, biome) %>%
  clipr::write_clip()

d_table %>%
  dplyr::arrange(id) %>%
  filter(!is.na(biome)) %>%
  mutate(id = paste('myxo_rpoB_', id, sep = '')) %>%
  select(location, site, lat, lon, predefined_habitat, biome, id) %>%
  clipr::write_clip()

table_1 <- select(d_table, id, site, location, lat, lon, predefined_habitat, sequenced_16s, sequenced_rpob) %>%
  arrange(predefined_habitat, site) %>%
  flextable() %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(id = "sample number",
                    lat = 'latitude',
                    lon = 'longitude',
                    predefined_habitat = 'predefined habitat',
                    sequenced_16s = "16s sequencing",
                    sequenced_rpob = 'rpoB sequencing') %>%
  hline(i = c(4, 9, 15, 20, 22, 27, 31, 37, 42, 48, 53, 58, 64, 67), border = fp_border_default()) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 10, part = 'all') %>%
  autofit() 

save_as_image(table_1, 'plots/manuscript_plots/table_s1.png')

# remove sites that are currently not present
d <- filter(d, !is.na(coordinates))

# split lat and long
d <- separate(d, coordinates, c('lat', 'lon'), sep = ',') %>%
  mutate(across(c(lat, lon), as.numeric))

# get average coordinates for each location
d_location <- group_by(d, site) %>%
  summarise(ave_lat = mean(lat),
            ave_lon = mean(lon),
            .groups = 'drop')
d_location

# make bounding box
box <- make_bbox(ave_lon, ave_lat, data = d_location, f = .10)

box[1] <- c(-5.8)
box[2] <- 49.9
box[3] <- -4
box[4] <- 50.57

# set up colours
cols <- tibble(group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_supratidal", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'),
               hab_order = c(1.1, 2.1, 1.2, 3.1, 3.2, 2.2, 2.3, 2.4, 1.3, 2.5, 2.6, 2.7, 1.4, 1.5, 2.8))
cols <- filter(cols, group %in% d$habitat_group_16s)
cols <- mutate(cols, habitat_group_16s = group) %>% arrange(hab_order)
cols <- left_join(cols, select(d_table, group = habitat_group_16s, predefined_habitat) %>% distinct()) %>%
  arrange(habitat_group_16s)

# change colours names
cols <- mutate(cols, group2 = gsub('_', ' ', group))

map1_base <- get_stamenmap(box, zoom = 11, maptype = "terrain") %>% 
  ggmap() +
  geom_point(aes(ave_lon, ave_lat), d_location, size = 3) +
  geom_point(aes(lon, lat, col = habitat_group_16s), d) +
  ggrepel::geom_label_repel(aes(ave_lon, ave_lat, label = site), d_location, size = MicrobioUoE::pts(14), box.padding = 0.7) +
  theme_bw(base_size = 14) +
  labs(x = 'Longitude',
       y = 'Latitude',
       title = '(a)') +
  scale_color_manual('Predefined habitat', values = setNames(cols$col, cols$habitat_group_16s), labels = cols$predefined_habitat)

map1 <- map1_base +
  ggforce::theme_no_axes() +
  annotate(y = box[2], x = box[3], label = 'Map made using ggmap, map tiles by Stamen Design, under CC BY 3.0', geom = 'label', size = MicrobioUoE::pts(8), vjust = -0.03, hjust = 1) +
  NULL

map1

map1 + 
  guides(colour = 'none')

ggsave('sampling_map.pdf', last_plot(), height = 6, width = 8, bg = 'transparent')

# create a map of the UK with Cornwall as a point

box2 <- box
box2[1] <- -10.8
box2[2] <- 49.6
box2[3] <- 1.9
box2[4] <- 59.7

cornwall <- tibble(x = -5.06971, y = 50.29189, label = 'Cornwall')

map2 <- get_stamenmap(box2, zoom = 7, maptype = "terrain") %>% 
  ggmap() +
  geom_point(aes(-5.06971, 50.29189)) +
  ggrepel::geom_label_repel(aes(x, y, label = label), cornwall, size = MicrobioUoE::pts(10), box.padding = 0.5) +
  ggforce::theme_no_axes() +
  theme(plot.margin = unit(c(0,0,-1,-1), 'mm'))
  
map3 <- map1 + inset_element(map2, 0.05, 0.6, 0.25, 0.93, ignore_tag = TRUE)

# save this plot 
saveRDS(map3, file.path('plots', 'map.rds'))



# zoom in on a single site
d_camel <- filter(d, location == 'Camel') %>%
  filter(id != '48')

# set up colours based on habitats
cols <- tibble(habitat_group = c("woodland_oak", "estuarine mud_low polyhaline", "woodland_pine", "reservoir", "river", "estuarine mud_oligohaline", "estuarine mud_full saline", "beach_supratidal", "pasture", "beach_subtidal","thrift_rhizosphere","beach_seaweed","field_wheat","rock_samphire","marine mud_full saline"),
               col = c('#089a2d', '#995a08', '#106c12', '#1170bd', '#9dcdf4', '#663c05', '#b6966b', '#f5e279', '#61dd1e', '#f2f426', '#c5f8ae', '#a6ab52', '#9ff121', '#5e8128', '#714a03'),
               hab_order = c(1.1, 2.1, 1.2, 3.1, 3.2, 2.2, 2.3, 2.4, 1.3, 2.5, 2.6, 2.7, 1.4, 1.5, 2.8))
cols <- filter(cols, habitat_group != 'thrift_rhizosphere') %>% arrange(hab_order)

# make bounding box
box <- make_bbox(lon, lat, data = d_camel, f = .10)

map2_base <- get_stamenmap(box, zoom = 13, maptype = "terrain") %>% 
  ggmap() +
  theme_bw(base_size = 14) +
  labs(x = 'Longitude',
       y = 'Latitude')

map2 <- map2_base +
  geom_point(aes(lon, lat, fill = habitat_group), col = 'black', shape = 21, d_camel, size = 6, show.legend = FALSE) +
  ggforce::theme_no_axes() +
  annotate(y = box[2], x = box[3], label = 'Map made using ggmap, map tiles by Stamen Design, under CC BY 3.0', geom = 'label', size = MicrobioUoE::pts(8), vjust = -0.03, hjust = 1) +
  scale_fill_manual('Predefined habitat', values = setNames(cols$col, cols$habitat_group)) +
  NULL

map2

ggsave('camel_map.pdf', last_plot(), height = 6, width = 8, bg = 'transparent')

d_meta_summary <- group_by(d_meta, habitat_group) %>%
  tally()

# look at the rpoB sequencing
readRDS('data/sequencing_rpoB/phyloseq/myxococcus/ps_phyloseq_myxo.rds') %>%
  phyloseq::sample_data() %>%
  data.frame() %>%
  pull(habitat_group_16s) %>%
  unique()
