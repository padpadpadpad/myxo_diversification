
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Processed data and analysis for the project exploring the macroevolution of Myxobacteria in natural environments.

## Summary

This repository contains the processed data and scripts used to
investigate the macroevolution of Myxobacteria in natural environments.
The repository has evolved through the project and now contains the
processed data and scripts used to recreate the analyses and figures in
the manuscript (as it is now). In some instances where the code to run
an analysis (e.g. for running BAMM) is not available here, we have tried
our best to highlight where that is the case.

## To use

The GitHub repo should be easy to use on everyone’s machine using
RStudio projects. RStudio projects automatically assigns the root
directory to the directory in which the project resides. Consequently
all of the analyses should be runnable without altering paths. Simply
double-click on the `.Rproj` folder that is in the base folder of this
repository and it should open the project in RStudio. All the packages
used in the script should be at the top of each script, and RStudio
usually requests you install them if they are not installed already.

I have avoided uploading overly large files to GitHub as there is a file
size limit of 50MB. The processed sequencing data is stored as
[phyloseq](https://joey711.github.io/phyloseq/) objects and should be
easy to read in. Raw sequencing files will be uploaded upon publishing
of the manuscript to the European Nucleotide Archive.

## Folder structure

### data

**data** includes the processed data from each part of the project.

- `sequencing_16s` contains the data for the 16s sequencing
  - The first dataset to be read in is **ps_16s_complete.rds**, this is
    the original phyloseq object from the processing pipeline. All other
    data files are produced from the scripts outlined below.
- `sequencing_rpoB` contains the data for the rpoB sequencing
  - The first dataset to be read in is
    **phyloseq/myxococcus/ps_phyloseq_myxo.rds**. This is the dataset of
    the phyloseq object after filtering has been done to only retain
    ASVs that are assigned to the Myxbacteria. All other data files are
    produced from this first dataset and the scripts outlined below.
- in the base folder there are some data files that are used across
  analyses
  - **metadata.csv** - completed metadata for each sample
  - **Supplemental Table 1 - Sampling.xlsx** - information on sampling
    sites used in **Table S1**.

### scripts

**scripts** includes the scripts to analyse the processed datasets.

- `sequencing_16s` contains the scripts used to analyse the 16s data.
  - The order of running the scripts is:
    1.  **prevalence_filtering.R** - removes low prevalence ASVs
    2.  **rarefy_data.R** - rarefies data and sets the root of the tree
        for the 16s data
    3.  **first_clustering.R** - runs initial PCoA plot to look at
        differences between pre-defined habitats using permutational
        ANOVAs
    4.  **second_clustering.R** - runs PCoA on data with three
        pre-defined habitats removed and does clustering analysis to
        identify biome clusters. Saves out plots that are used in
        **Figure_1.R** to make **Figure 1**.
  - **extra_functions.R** contains extra functions used in the other
    scripts.
- `sequencing_rpoB` contains the scripts used to analyse the rpoB data.
  - The scripts in `processing` and `tree_building` can be run, then
    followed by the scripts in `analysis`.
    - The order of the processing scripts is as follows:
      1.  **remove_mmseqs_nonmyxo.R** - removes ASVs from the rpoB
          dataset that are not assigned to Myxobacteria in both MMSeqs2
          and the naive bayesian classifier used by dada2.
      2.  **asvs_to_otus.R** - does prevalence filtering on the rpoB
          Myxobacteria dataset, creates a sequence alignment of the ASV
          dataset and then creates new datasets at different levels of
          OTU similarity.
      3.  **assign_habitat_preference_all.R** - assigns biome preference
          to each ASV in each OTU similarity cut-off dataset.
      4.  **align_myxo.R** - creates a sequence alignment for each
          Myobacteria dataset across each OTU similarity and creates a
          backbone tree for that alignment using
          **tree_processing/csv2constraint.pl**.
      5.  **rarefaction_curves.R** - creates rarefaction curves from the
          Myxobacteria dataset created at the end of
          **remove_mmseqs_nonmyxo.R**. Recreates **Figure S1**.
      6.  **check_habitat_proportions.R** looks at how many ASVs are in
          each biome preference across different OTU similarity
          cut-offs. Recreates **Figure S2**.
      7.  **check_myxo_clustering.R** runs a PCoA on the prevalence
          filtered Mxyobacteria ASV dataset, looking for clustering of
          communities based on their different biomes. Recreates
          **Figure S3**.
    - The order of the tree building scripts is as follows:
      1.  **run_raxml_many_screens.sh** - runs raxml on all the sequence
          alignments.
      2.  **append_tiplabels.R** - reads the tree and taxonomic
          information into R and runs **extend_labels.pl**. Each tree
          after this step was rooted manually using FigTree. A key step
          for **chronos()** to work in the next script is to include
          annotations as Nexus when exporting the tree from FigTree.
      3.  **make_tree_ultrametric.R** - makes each tree ultrametric and
          remove family labels from tree.
      4.  **csv2constraint.pl** creates a constraint tree from a csv
          specifiying groupings and a list of ASVs and their grouping in
          relation to the constraint tree.
      5.  **extend_labels.pl** adds the family name onto the tip label
          for the tree to make rooting in FigTree easier.
    - The order of the analysis scripts is as follows.

    3.  
    4.  
    5.  
    6.  
    7.  
  - 
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
- `sequencing_16s` contains extra figures and tables produced during the
  analysis of the 16s data
- `sequencing_rpoB` contains extra figures and tables produced during
  the analysis of the rpoB data
- in the base folder there are some files that do not belong to either
  16s or rpoB analyses
  - **map_legend.rds** - saved legend used to create **Figure 1**.
  - **map.rds** - saved plot used to create **Figure 1**.
  - **rough_manuscript_plots.ppt** - rough slides planning out the
    figures for the manuscript.
- **sequencing 16s** - folder containing data from the 16s sequencing

1.  
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