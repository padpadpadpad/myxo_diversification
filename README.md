
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Processed data and analysis for the project exploring the macroevolution of Myxobacteria in natural environments.

## Summary

This repository contains the processed data and scripts used to
investigate the macroevolution of Myxobacteria in natural environments.
The repository has evolved through the project and now contains the
processed data and scripts used to recreate the analyses and figures in
the manuscript (as it is now). In some instances where the code to run
an analysis (e.g. for running BAMM) is not available here, we have tried
our best to highlight that fact.

**This repository contains a lot of previous code and scripts that are
no longer used, but that we did not want to remove completely as they
might be of use to somebody. They live in a variety of `previous`
folders in their respective analysis folders.**

## To use

The GitHub repo should be super easy to use on everyone’s machine using
RStudio projects. RStudio projects automatically assigns the root
directory to the directory in which the project resides. Consequently
all of the analyses should be runnable without altering paths. Simply
double-click on the `.Rproj` folder that is in the base folder of this
repository and it should open the project in RStudio.

I have avoided uploading overly large files to GitHub as there is a file
size limit of 50MB. The processed sequencing data is stored as
[phyloseq](https://joey711.github.io/phyloseq/) objects and should be
easy to read in. Raw sequencing files will be uploaded upon publishing
of the manuscript.

## Folder structure

### data

**data** includes the processed data from each part of the project.

- `sequencing_16s` contains the data for the 16s sequencing
- `sequencing_rpoB` contains the data for the rpoB sequencing

### scripts

**scripts** includes the scripts to analyse the processed datasets.

- `sequencing_16s` contains the scripts used to analyse the 16s data.
  1.  `prevalence_filtering.R` - removes low prevalence ASVs
  2.  `rarefy_data.R` - rarefies data and sets root to tree
  3.  `first_clustering.R` - runs initial PCoA plot to look at
      differences between habitats
  4.  `second_clustering.R` - runs cleaned PCoA and clustering analysis
      to identify habitat clusters in the data.
- `sequencing_rpoB` contains the scripts used to analyse the rpoB data.
- `misc` contains some miscellaneous scripts that are not used much in
  the final analysis.
  - **sampling_sites.R** makes **Table_S1** and saves a map out that is
    used in **Figure_1.R** to make Figure 1.
  - all other scripts are not used in the final analyses and are not of
    particular importance.

### plots

**plots** contains figures and tables made during the project.

- `manuscript_plots` contains the figures and tables used in the
  manuscript

- `sequencing_16s` contains figures and tables produced during the
  analysis of the 16s data

- `sequencing_rpoB` contains figures and tables produced during the
  analysis of the rpoB data

- **metadata.csv** - completed metadata for each sample

- **sample_cluster_assignments.csv** - cluster assignments for each
  sample based on 16S clustering.

- **sequencing 16s** - folder containing data from the 16s sequencing

1.  `ps_16s_complete.rds` - the original phyloseq object from the
    processing pipeline
2.  `ps_16s_prev_filtered.rds` - low abundance and low prevalence ASVs
    removed
3.  `ps_16s_low_depth_removed.rds` - low depth samples removed **no
    samples were too low here**
4.  `ps_16s_rarefied.rds` - rarefied samples

- **sequencing_rpoB** - folder containing data from the rpoB sequencing

1.  `alignments` - folder containing the alignments of myxobacteria
    reads after clustering at different levels of OTU similarity and
    prevalence filtering
2.  `raxml` - the trees of each level of OTU similarity. To see the
    meaning of the output files see the [raxml-ng help
    page](https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings).
    Rooted trees (myxo_asv.raxml.bestTreeRooted) were made in R using
    `ape::root` and has only been done for the asv tree so far.
    `_LogLik.txt` files show the log likelihood of all the trees
    outputted by the raxml command.

### scripts

**scripts** includes the scripts from each part of the project so far.
`sequencing_16s` contains the scripts for the 16s sequencing etc. At the
top of each script there is a summary of what has been done in each
script. Order of scripts to be run is the same as numbered here.

- **sequencing_16s**

1.  `prevalence_filtering.R` - removes low prevalence ASVs
2.  `rarefy_data.R` - rarefies data and sets root to tree
3.  `first_clustering.R` - runs initial PCoA plot to look at differences
    between habitats
4.  `second_clustering.R` - runs cleaned PCoA and clustering analysis to
    identify habitat clusters in the data.

- **sequencing_rpoB**

### plots

**plots** includes the plots from each part of the project so far.
`sequencing_16s` contains the plots for the 16s sequencing etc.
