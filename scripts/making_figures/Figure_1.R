# assemble Figure 1

# load in packages
librarian::shelf(tidyverse, patchwork)

# load in plots
p1 <- readRDS('plots/map.rds') &
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16)) 
p2 <- readRDS('plots/sequencing_16s/ordination_16S_clean.rds') +
  theme(legend.position = 'bottom') +
  guides(col = 'none',
         shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  ggtitle('(b)')
p3 <- readRDS('plots/sequencing_16s/cluster_pcoa.rds') +
  theme(legend.position = 'bottom') +
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
  ggtitle('(c)')

# grab legend from p1
legend <- cowplot::get_legend(p1)
cowplot::plot_grid(legend)

p4 <- p1 &
  guides(colour = 'none')

# create patchwork layout
layout <- '
AAD
BCD
'
p4 + p2 + p3 + legend + plot_layout(design = layout, heights = c(0.7, 0.3))

ggsave('plots/manuscript_plots/Figure_1.png', last_plot(), height = 10, width = 14)
