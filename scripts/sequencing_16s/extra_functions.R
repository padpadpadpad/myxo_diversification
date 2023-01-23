# extra functions used in amplicon 16S analysis.R

# convert distance matrix to dataframe
dist_2_df <- function(dist_ob){
  m <- as.matrix(dist_ob) # coerce dist object to a matrix
  xy <- t(combn(colnames(m), 2))
  return(data.frame(xy, dist=m[xy], stringsAsFactors = FALSE))
}

# code stolen from phyloseq website
get_top_taxa <- function(ps, tax_rank, to_keep){
  temp <- tapply(phyloseq::taxa_sums(ps), phyloseq::tax_table(ps)[, tax_rank], sum, na.rm = TRUE)
  temp2 <-  names(sort(temp, TRUE))[1:to_keep]
  return(temp2)
}

# get betadisper dataframes

# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get distance from 00
dist_between_points <- function(x1, x2, y1, y2){
  return(abs(sqrt((x1 - x2)^2+(y1-y2)^2)))
}

# calculate pairwise permanovas
calc_pairwise_permanovas <- function (dm, metadata_map, compare_header, extra_header = NULL, n_perm) 
{
  comp_var = metadata_map[, compare_header]
  comp_pairs = combn(levels(comp_var), 2)
  pval = c()
  R2 = c()
  for (i in 1:ncol(comp_pairs)) {
    pair = comp_pairs[, i]
    dm_w_map = list(dm_loaded = dm, map_loaded = metadata_map)
    dm_w_map$map_loaded$in_pair = comp_var %in% pair
    dm_w_map_filt = mctoolsr::filter_dm(dm_w_map, filter_cat = "in_pair", 
                                        keep_vals = TRUE)
    if (!missing(n_perm)) {
      if(!is.null(extra_header)){
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, compare_header] + dm_w_map_filt$map_loaded[, extra_header], permutations = n_perm)
      }
      m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, compare_header], permutations = n_perm)
    }
    else {
      if(!is.null(extra_header)){
        m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, compare_header] + dm_w_map_filt$map_loaded[, extra_header])
      }
      m = vegan::adonis2(dm_w_map_filt$dm_loaded ~ dm_w_map_filt$map_loaded[, compare_header])
    }
    pval = c(pval, m$`Pr(>F)`[1])
    R2 = c(R2, m$R2[1])
  }
  results = data.frame(t(comp_pairs), R2, pval)
  results$pvalBon = p.adjust(pval, 'bonferroni')
  results$pvalFDR = p.adjust(pval, 'fdr')
  results
}

# custom rarefaction curve
calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=FALSE) {
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- reshape2::melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

