
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Processed data and analysis for the project exploring the macroevolution of Myxobacteria in natural environments.

## Summary

This repository contains the processed data and scripts used to
investigate the macroevolution of Myxobacteria in natural environments.
The repository has evolved through the project and now contains the
processed data and scripts used to recreate the analyses and figures in
the manuscript (as it is now, for the revised manuscript). We have tried
our best to highlight the instances where the code to run an analysis
(e.g. for running BAMM) is not available in this repository.

It will be cleaned up of intermediate and surplus files upon acceptance
of the manuscript.

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
  - `phyloseq` contains the phyloseq objects from which all analyses are
    done, as well as habitat preference assignments for each ASV at each
    level of OTU similarity used.
  - `raxml` contains the input files to raxml for the tree building and
    the output files from raxml and treePL.
  - `bamm` contains the input and output files for the BAMM analyses.
  - `processed` contains saved intermediate data files and models from
    the discrete character evolution models and models from *secsse*.
  - `trees` contains the final trees used in the analyses. We have
    copied them over from their original folders and renamed them to be
    simpler than the names of the files used in the other scripts.
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
      6.  **visualise_biome_preference.R** - visualises the biome
          preference of several ASV across to demonstrate how the
          bootstrap method worked. Recreates **Figure S2**.
      7.  **check_habitat_proportions.R** looks at how many ASVs are in
          each biome preference across different OTU similarity
          cut-offs. Recreates **Figure S3**.
      8.  **check_myxo_clustering.R** runs a PCoA on the prevalence
          filtered Mxyobacteria ASV dataset, looking for clustering of
          communities based on their different biomes. Recreates
          **Figure S4**.
    - The order of the tree building scripts is as follows:
      1.  **csv2constraint.pl** creates a constraint tree from a csv
          specifiying groupings and a list of ASVs and their grouping in
          relation to the constraint tree.
      2.  **run_raxml_many_screens.sh** - runs *raxml* on all the
          sequence alignments.
      3.  **extend_labels.pl** adds the family name onto the tip label
          for the tree to make rooting in FigTree easier.
      4.  **bootstrap_raxml.sh** - runs *raxml* bootstrapping on the
          best tree for the ASV tree, 95% and 97.7% trees.
      5.  **append_tiplabels.R** - reads the tree and taxonomic
          information into R and runs **extend_labels.pl**. Each tree
          after this step was rooted manually using FigTree, and the
          rooted tree was made ultrametric using *treePL*.
      6.  **pick_boostrap_trees.R** - picks 9 random bootstrapped trees
          from the ASV bootstrap replicates.
      7.  **treepl_bootstraps.sh** shows how *treePL* was ran on the ASV
          bootstrap replicate trees.
      8.  **match_treepl_with_bootstrapped_trees.R** - matches the
          output of the best tree from *treePL* with the non-ultrametric
          tree from *raxml* bootstrapping for the ASV tree, 95% and
          97.7% trees.
      9.  **bootstrap_force_ultrametric.R** - makes the output of
          *treePL* ultrametric using **phytools::force.ultrametric()**.
      10. **compare_asv_trees.R** - compares the best ASV tree to the 9
          bootstrap replicate trees.
      11. **ltt_plot_new_old.R** - compares the lineage through time
          plot between *ape::chronos()* and *treePL*.
    - The analysis scripts are split into folders for `asv`,
      `otu_level`, and `bootstraps` when essentially the same scripts
      are replicated for the ASV tree, the 95% and 97.7% trees, and
      bootstrap replicate trees respectively. The base folder contains
      scripts that aggregate the results from the other folders, as well
      as helper scripts. I will summarise the scripts from the `asv`
      folder, with their counterparts in the other folders being
      similar.
      1.  **visualise_phylogeny.R** in the `asv` folder plots the
          phylogeny and groups some biome preferences together.
          Recreates **Figure 2**. It also plots the bootstrap values of
          the trees (**Figures S5** to **S13**).
      2.  **diversitree_mk_models.R** in the `asv` folder fits the
          Markov models looking at biome preference transitions.
          Recreates **Figure 3**, **Table 1**, **Figure S14**, **Figure
          S15**, **Table S2**. This script has its equivalents in the
          `otu_level` and `bootstraps` folders.
      3.  **bamm_analysis.R** in the `asv` folder analyses the output of
          BAMM. Recreates **Figure 4**.
      4.  **fit_secsse\_\*\_start_vals** finds the maximum likelihood of
          potential start values for the musse, muhisse, ctd2, ctd3, and
          ctd4 models fit using secsse. They also create tables
          demonstrating the model parameters of each model (**Tables
          S7** to **S11**). These scripts have their equivalents in the
          `otu_level` and `bootstraps` folders.
      5.  **fit_secsse\_\*\_v2.R** fits the musse, muhisse, ctd2, ctd3,
          and ctd4 models using secsse. These scripts have their
          equivalents in the `otu_level` and `bootstraps` folders.
      6.  **check_fixed_transition_rates.R** looks at the impact on
          parameter estimates of fixing some transition rates in the
          fitting of models using secsse.
      7.  **diversitree_helper_functions.R** contains helper functions
          used to summarise and visualise the output of the diversitree
          models.
      8.  **bamm_analysis_all.R** aggregates the output of BAMM for the
          ASV, 95% and 97.7% trees. Recreates **Figures S28** to
          **S31**.
      9.  **compare_transition_rates.R** compares the transition rates
          between between ASV tree, bootstrap trees, 95% and 97.7%
          trees. Recreates **Figure S27** and **Table S3**.
      10. **diversitree_musse_models.R** fits MuSSE models using
          diversitree to the ASV, 95% and 97.7% trees.
      11. **run_bamm.sh** runs BAMM on the ASV, 95% and 97.7% trees.
      12. **seccse_model_compare.R** compares the secsse models at the
          ASV level, bootstraps, and 95% and 97.7% trees. Recreates
          **Figure S32**, **Table S4**, **Table S5** and **Table S6**.
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

### Software versions

- These scripts were ran on R version 4.3.1 and MacOS Ventura 13.4.1.
